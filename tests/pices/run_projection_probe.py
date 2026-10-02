"""Frozen P5 initialization: no advection, heat solve, source, TA update or Stokes.
Native double snapshots are independently checked with NumPy/SciPy.
The consistent-matrix projection is an offline candidate, not a solver change.
"""
import argparse,gc,gzip,hashlib,json,os,re,shutil,subprocess,time
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix,diags
from scipy.sparse.linalg import cg
from scipy.optimize import minimize

def audit(d,profile):
 blocks=[np.fromfile(d/f'probe.{rank}.nodes.bin',dtype='=f8').reshape(-1,8) for rank in range(12)]
 xyz=np.concatenate([b[:,:3] for b in blocks]);_,mapping=np.unique(np.round(xyz,10),axis=0,return_inverse=True);N=int(mapping.max()+1)
 allnodes=np.concatenate(blocks);count=np.bincount(mapping,minlength=N)
 initial=np.bincount(mapping,weights=allnodes[:,3],minlength=N)/count
 fixed=(np.bincount(mapping,weights=allnodes[:,6],minlength=N)>0)
 mass=np.bincount(mapping,weights=allnodes[:,7],minlength=N)
 assert np.max(abs(allnodes[:,3]-initial[mapping]))<2e-12
 assert np.all((mass>0)&np.isfinite(mass))
 ids=[];ww=[];tp=[];offset=0
 for rank,b in enumerate(blocks):
  a=np.fromfile(d/f'probe.{rank}.particles.bin',dtype='=f8').reshape(-1,17)
  local=a[:,1:17:2].astype(int)-1;assert local.min()>=0 and local.max()<len(b)
  ids.append(mapping[offset:offset+len(b)][local]);ww.append(a[:,2:17:2]);tp.append(a[:,0]);offset+=len(b)
 ids=np.concatenate(ids);ww=np.concatenate(ww);tp=np.concatenate(tp);P=len(tp)
 assert np.isfinite(ww).all() and ww.min()>=0 and np.max(abs(ww.sum(axis=1)-1))<2e-12
 W=csr_matrix((ww.ravel(),(np.repeat(np.arange(P),8),ids.ravel())),shape=(P,N))
 q_error=float(np.max(abs(W@initial-tp)));assert q_error<2e-12
 den=np.asarray(W.sum(axis=0)).ravel();assert (den>0).all()
 raw=np.asarray(W.T@tp).ravel()/den;clamped=raw.copy();clamped[fixed]=initial[fixed]
 p_error=float(np.max(abs(raw[mapping]-allnodes[:,4])));bc_error=float(np.max(abs(clamped[mapping]-allnodes[:,5])))
 assert p_error<2e-12 and bc_error<2e-12
 energies=np.array([mass@initial,mass@raw,mass@clamped]);actual=np.loadtxt(d/'probe.energy.txt')
 assert np.max(abs(energies-actual[:3]))<3e-12
 defect=float(np.max(abs(clamped-initial)));assert abs(defect-actual[3])<2e-12
 if profile=='constant':assert defect<2e-12
 # Consistent least-squares candidate: enforce only the same fixed grid nodes.
 A=(W.T@W).tocsr();free=np.flatnonzero(~fixed);bound=np.flatnonzero(fixed)
 Aff=A[free][:,free].tocsr();Afb=A[free][:,bound].tocsr();M=diags(1/Aff.diagonal())
 def solve(values,boundary):
  rhs=np.asarray(W.T@values).ravel()[free]-Afb@boundary[bound];iterations=[0]
  def tick(_):iterations[0]+=1
  t=time.monotonic();z,info=cg(Aff,rhs,M=M,rtol=1e-12,atol=1e-15,maxiter=2000,callback=tick)
  assert info==0,info
  field=boundary.copy();field[free]=z
  return field,iterations[0],time.monotonic()-t
 consistent,it,seconds=solve(tp,initial)
 cerror=float(np.max(abs(consistent-initial)));assert cerror<1e-8,cerror
 result=dict(global_nodes=N,particles=P,Q_reference_max_error=q_error,P_reference_max_error=p_error,boundary_reference_max_error=bc_error,initial_energy=float(energies[0]),raw_projection_energy=float(energies[1]),clamped_projection_energy=float(energies[2]),raw_change_percent=float(100*(energies[1]-energies[0])/energies[0]),boundary_correction_percent=float(100*(energies[2]-energies[1])/energies[0]),net_change_percent=float(100*(energies[2]-energies[0])/energies[0]),grid_max_error_K=3400*defect,grid_volume_weighted_rms_K=float(3400*np.sqrt(mass@((clamped-initial)**2)/mass.sum())),consistent_candidate=dict(max_error_K=3400*cerror,energy_change_percent=float(100*(mass@consistent-energies[0])/energies[0]),iterations=it,solve_seconds=seconds))
 if profile=='constant':
  # Non-FE particle step profile bounded in [0,1]. Tests a candidate's monotonicity,
  # not temporal accuracy: no physical evolution or exact solution is claimed.
  radii=np.linalg.norm(xyz,axis=1);rnodal=np.bincount(mapping,weights=radii,minlength=N)/count
  jump=(W@rnodal<.775).astype(float);boundary=(rnodal<.775).astype(float)
  candidate,its,_=solve(jump,boundary);lumped=np.asarray(W.T@jump).ravel()/den;lumped[fixed]=boundary[fixed]
  # A bounded QP is only an offline feasibility probe; it supplies no physical
  # particle masses and makes no moving-particle energy-conservation claim.
  target=np.asarray(W.T@jump).ravel()[free]-Afb@boundary[bound]
  def objective(x):
   ax=Aff@x
   return .5*float(x@ax)-float(target@x),ax-target
  opt=minimize(objective,lumped[free],jac=True,method='L-BFGS-B',bounds=[(0.,1.)]*len(free),options=dict(ftol=1e-15,gtol=1e-9,maxiter=1000,maxls=40))
  assert opt.success,opt.message
  bounded=boundary.copy();bounded[free]=opt.x;gradient=Aff@opt.x-target
  projected_gradient=opt.x-np.clip(opt.x-gradient,0.,1.)
  kkt=float(np.max(abs(projected_gradient)));assert kkt<1e-5,kkt
  assert bounded.min()>=0 and bounded.max()<=1
  result['bounded_sharp_candidate']=dict(min=float(bounded.min()),max=float(bounded.max()),min_K=float(300+3400*bounded.min()),max_K=float(300+3400*bounded.max()),iterations=int(opt.nit),projected_gradient_inf=kkt,FE_energy=float(mass@bounded),lumped_FE_energy=float(mass@lumped),unbounded_FE_energy=float(mass@candidate),physical_energy_conservation_proven=False)
  result['sharp_particle_profile']=dict(particle_min=float(jump.min()),particle_max=float(jump.max()),lumped_min=float(lumped.min()),lumped_max=float(lumped.max()),consistent_min=float(candidate.min()),consistent_max=float(candidate.max()),consistent_min_K=float(300+3400*candidate.min()),consistent_max_K=float(300+3400*candidate.max()),consistent_iterations=its)
 return result

def control(build,root,results):
 d=root/'assim_full_step_control';d.mkdir();(d/'DATA').mkdir();source=root/'assim_base'
 cfg=(source/'case.cfg').read_text()
 for key,value in [('maxstep','1'),('maxtotstep','2'),('storage_spacing','1'),('CBF_frequency','1')]:
  cfg=re.sub(r'^'+key+r'=.*$',key+'='+value,cfg,flags=re.M)
 (d/'case.cfg').write_text(cfg);shutil.copy(source/'refstate.txt',d/'refstate.txt');shutil.copytree(source/'forcing',d/'forcing')
 with (d/'stdout').open('w') as out,(d/'stderr').open('w') as err:
  run=subprocess.run([os.environ.get('MPIEXEC','/usr/local/bin/mpiexec'),'--oversubscribe','-n','12',str(build.resolve()/'CitcomSFull'),'case.cfg'],cwd=d,stdout=out,stderr=err,timeout=600)
 assert run.returncode in (0,8),(run.returncode,(d/'stderr').read_text()[-1200:])
 return analyze_control(root,results,run.returncode)

def analyze_control(root,results,exit_code):
 d=root/'assim_full_step_control';source=root/'assim_base'
 max_particle_error=0.
 for rank in range(12):
  expected=np.fromfile(source/f'probe.{rank}.particles.bin',dtype='=f8').reshape(-1,17)[:,0]
  with gzip.open(d/f'DATA/{rank}/0/tracer.{rank}.0.gz','rt') as f:
   header=f.readline();actual=np.array([float(line.split()[3]) for line in f])
  assert actual.shape==expected.shape
  # Output_gzdir writes tracer values with %.5e. Compare like precision.
  rounded=np.array([float(format(x,'.5e')) for x in expected]);assert np.array_equal(actual,rounded)
  max_particle_error=max(max_particle_error,float(np.max(abs(actual-expected))))
 ls=(d/'DATA/0/log').read_text().splitlines();steps=[dict(t.split('=',1) for t in l.split()[1:]) for l in ls if l.startswith('PICES_STEP ')]
 assert len(steps)==1 and steps[0]['step']=='1'
 metric=[dict(t.split('=',1) for t in l.split()[1:]) for l in ls if l.startswith('P5_METRIC step=0 ')][0]
 initial=float(metric['energy_proxy']);assert abs(initial-results['assim_base']['initial_energy'])<3e-12
 full=100*float(steps[0]['remap_energy'])/initial;frozen=results['assim_base']['net_change_percent']
 return dict(mpi_exit=exit_code,initial_particle_T_matches_printed_precision=True,initial_particle_T_max_rounding_difference=max_particle_error,full_first_remap_percent=full,frozen_remap_percent=frozen,difference_percentage_points=full-frozen,frozen_to_full_magnitude_ratio=abs(frozen/full),note='Same-build one-step P5 assimilation control; initial particle temperatures matching at ASCII output precision. Frozen probe itself executes no temporal operators.')

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('build',type=Path);p.add_argument('runs',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True)
 cases=[]
 for scenario in ['cold','hot','assim']:
  for variant in ['grid_coarse','base','grid_fine','np_low','np_high']:cases.append((f'{scenario}_{variant}',f'{scenario}_pices_{variant}','keep'))
 cases += [(mode,'cold_pices_base',mode) for mode in ['constant','cartesian_linear']]
 src=a.runs/'pices_p5';manifest={x['name']:x for x in json.loads((src/'matrix.json').read_text())['cases']};results={}
 for name,input_name,profile in cases:
  row=manifest[input_name];d=a.output/name;d.mkdir();(d/'DATA').mkdir()
  shutil.copy(src/'cases'/f'{input_name}.cfg',d/'case.cfg');kind='assim' if row['scenario']=='assim' else 'transport';shutil.copy(src/f'refstate_{kind}_{row["nodes"]}.txt',d/'refstate.txt');shutil.copytree(src/f'forcing_{row["nodes"]}',d/'forcing',ignore=shutil.ignore_patterns('._*'))
  with (d/'stdout').open('w') as out,(d/'stderr').open('w') as err:
   run=subprocess.run([os.environ.get('MPIEXEC','/usr/local/bin/mpiexec'),'--oversubscribe','-n','12',str(a.build.resolve()/'PICESProjectionProbe'),'case.cfg'],cwd=d,env=dict(os.environ,PICES_PROBE_PROFILE=profile),stdout=out,stderr=err,timeout=180)
  assert run.returncode==0,(name,run.returncode,(d/'stderr').read_text()[-1600:])
  for rank in range(12):
   log=(d/f'DATA/{rank}/log').read_text();assert 'PICES_STEP ' not in log and 'PICES_TA ' not in log and 'PICES_STOKES ' not in log
  result=audit(d,profile);assert result['particles']==12*(row['nodes']-1)**3*row['particles']
  result['input_case']=input_name;result['profile']=profile;result['mpi_exit']=run.returncode;results[name]=result
  (a.output/'partial.json').write_text(json.dumps(results,indent=2)+'\n');print(name,'net energy %',result['net_change_percent'],'max error K',result['grid_max_error_K'],flush=True);gc.collect()
 comparison=control(a.build,a.output,results)
 report=dict(full_step_control=comparison,status='PASS',scope='Frozen production P/Q diagnostic plus offline consistent-projection candidate',production_algorithm_changed=False,production_approved=False,cases=results,binary_sha256=hashlib.sha256((a.build/'PICESProjectionProbe').read_bytes()).hexdigest(),numpy=np.__version__,source_sha256={str(f.relative_to(Path(__file__).resolve().parents[2])):hashlib.sha256(f.read_bytes()).hexdigest() for f in [Path(__file__).resolve(),Path(__file__).with_name('projection_probe.inc'),Path(__file__).with_name('build_projection_probe.py'),Path(__file__).resolve().parents[2]/'lib/Pices.c']})
 (a.output/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
