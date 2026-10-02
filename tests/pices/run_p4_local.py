"""P4 full-kernel tests: subcycling, exact TA/Q transfer, mapping convergence, CBF."""
import argparse,json,math,os,re,subprocess
from pathlib import Path
import numpy as np
from verify_p2 import fields
from verify_p4 import boundary
p=argparse.ArgumentParser();p.add_argument('build',type=Path);p.add_argument('cfg',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True)
base=a.cfg.read_text();results={};consistent='pices_projection=bounded_consistent' in base
def run(name,nz=5,ta=True,source=250000,resolved=False):
 d=a.output/name;d.mkdir();(d/'DATA').mkdir();s=base
 opts=dict(maxstep=1,maxtotstep=2,tracers_per_element=32,pices_checkpoint='off',Q0=source,dissipation_number=0,lith_age=int(ta),lith_age_asml=int(ta),lith_age_asml_tau_Ma=.1,
 phase_delta_rho='0,0,0',phase_delta_s='0,0,0',fixed_timestep='1e-7',kT_exponent=0,kd_upper_prefactor=4,kd_lower_prefactor=4,
 nodex=nz,nodey=nz,nodez=nz,mgunitx=nz-1,mgunity=nz-1,mgunitz=nz-1,max_plate_age_Ma=20000 if resolved else 70)
 for k,v in opts.items():s=re.sub(r'^'+k+'=.*$',k+'='+str(v),s,flags=re.M)
 (d/'case.cfg').write_text(s);ref=re.search(r'^refstate_file=(.+)$',s,re.M)[1];(d/ref).write_text(('1 1 .6 1 1\n')*nz)
 forcing=d/'pices_p4_forcing';forcing.mkdir()
 for age in range(4):
  (forcing/f'trench.{age}.xyz').write_text('')
  for cap in range(12):(forcing/f'age.{age}.{cap}').write_text((str(10000 if resolved else 30+10*age)+'\n')*nz*nz)
 with (d/'stdout').open('w') as out,(d/'stderr').open('w') as err:
  r=subprocess.run([os.environ.get('MPIEXEC','/usr/local/bin/mpiexec'),'--oversubscribe','-n','12',str(a.build.resolve()/'PICESManufactured'),'case.cfg'],cwd=d,env=dict(os.environ,PICES_TEST_MODE='p4_linear' if resolved else 'constant'),stdout=out,stderr=err,timeout=180)
 assert r.returncode in (0,8),(name,r.returncode,(d/'stderr').read_text()[-1000:])
 assert all('PICES_STEP step=1 ' in (d/f'DATA/{rank}/log').read_text() for rank in range(12)),name
 print(name,'completed',flush=True);return d

heat=run('heat_only',ta=False);ta=run('heat_and_TA');err=0.;qerr=0.;mapping_errors=[];local_mapping=[]
for rank in range(12):
 before=np.loadtxt(heat/f'nodes.{rank}.1.txt');after=np.loadtxt(ta/f'nodes.{rank}.1.txt');data=np.loadtxt(ta/f'ta.{rank}.1.txt')
 h={k:float(v) for k,v in fields((ta/f'ta.{rank}.1.txt').read_text().splitlines()[0]).items()}
 expected=np.zeros(len(data));delta=after[:,4]-before[:,4]
 for n,(_,inc,depth,tref,age,flag) in enumerate(data):
  if before[n,5] or depth>=h['depth']:continue
  age=min(age,h['cap']/h['scalet']);target=tref-h['surface']*math.erfc(.5*depth/math.sqrt(age))
  w=math.expm1(-h['exp']*depth/h['depth'])/math.expm1(-h['exp'])
  # The legacy TA time-scale product is evaluated in float before w.
  physical_dt=float(np.float32(np.float32(h['dt'])*np.float32(h['scalet'])))
  alpha=-math.expm1(-physical_dt*(1-w)/h['tau'])
  expected[n]=alpha*(target-before[n,4])
 err=max(err,float(np.max(abs(delta-expected))),float(np.max(abs(data[:,1]-delta))))
 oldp=np.loadtxt(heat/f'particles.{rank}.1.txt');newp=np.loadtxt(ta/f'particles.{rank}.1.txt')
 assert np.array_equal(oldp[:,1:],newp[:,1:]),'TA moved particles'
 ids=newp[:,1:17:2].astype(int)-1;weights=newp[:,2:17:2];q=(delta[ids]*weights).sum(axis=1)
 qerr=max(qerr,float(np.max(abs(newp[:,0]-oldp[:,0]-q))))
 num=np.zeros(len(delta));den=np.zeros(len(delta))
 for j in range(8):
  np.add.at(num,ids[:,j],weights[:,j]*q);np.add.at(den,ids[:,j],weights[:,j])
 local_mapping.append((after[:,1:4],delta,num,den))
 log=(ta/f'DATA/{rank}/log').read_text();assert log.count('PICES_TA ')==1
 assert int(fields(next(l for l in log.splitlines() if l.startswith('PICES_STEP ')))['substeps'])>=3
 # Cached residual must be independent of TA, including byte-identical native output.
 for side in ['surf','botm']:assert (heat/f'DATA/{rank}/q.{side}.{rank}.1').read_bytes()==(ta/f'DATA/{rank}/q.{side}.{rank}.1').read_bytes()
assert err<2e-14 and qerr<2e-14,(err,qerr)
# Sum all cap contributions at common Cartesian nodes before normalization.
assembled={}
for xyz,delta,num,den in local_mapping:
 for x,t,u,v in zip(xyz,delta,num,den):
  key=tuple(np.round(x,10));row=assembled.setdefault(key,[t,0.,0.]);assert abs(row[0]-t)<1e-12
  row[1]+=u;row[2]+=v
mapping=max(abs(u/v-t) for t,u,v in assembled.values())
reported=float(fields(next(l for l in (ta/'DATA/0/log').read_text().splitlines() if l.startswith('PICES_TA ')))['mapping_error'])
if consistent:
 assert reported<1e-11,reported
 mapping=reported # Independent consistent operator is checked by run_consistent_projection.py.
else:assert abs(mapping-reported)<1e-13,(mapping,reported)
rows=(ta/'DATA/0/log').read_text().splitlines();ledger=fields(next(x for x in rows if x.startswith('PICES_EBA ')))
reaction=float(ledger['boundary_reaction'])/h['dt']*4*3400*6371000
flux=boundary(ta,1,'botm')-boundary(ta,1,'surf');balance=abs(flux-reaction)/max(1,abs(flux)+abs(reaction));assert balance<1e-12,balance
results['multisubstep_TA_and_CBF']=dict(projection='bounded_consistent' if consistent else 'lumped',mapping_error=reported,status='PASS',TA_grid_error=err,particle_Q_error=qerr,mapping_diagnostic_error=None if consistent else abs(mapping-reported),CBF_before_after_TA_equal=True,relative_CBF_boundary_residual=balance)
# Mapping is intentionally not fed back to T. Verify P(Q(delta))-delta shrinks
# under radial/horizontal refinement with the same physical TA law.
# A synthetic 10000 Ma age resolves the HSC layer on these coarse meshes;
# a continuous linear initial T gives a continuous increment at the top BC.
grids=[5,9] if consistent else [5,9,17,25]
for nz in grids:
 d=run('mapping_'+str(nz),nz,source=0,resolved=True)
 log=(d/'DATA/0/log').read_text().splitlines();row=fields(next(x for x in log if x.startswith('PICES_TA ')))
 mapping_errors.append(float(row['mapping_error']))
if consistent:assert max(mapping_errors)<1e-11,mapping_errors
else:assert mapping_errors[2]<min(mapping_errors[:2]) and mapping_errors[3]<.85*mapping_errors[2],mapping_errors
results['mapping_refinement']=dict(status='PASS',nodes=grids,max_errors=mapping_errors,coarse_to_fine_ratio=mapping_errors[0]/max(1e-30,mapping_errors[-1]))
(a.output/'summary.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(results,indent=2))
