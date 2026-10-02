"""Independent matrix and signed-field checks for the distributed bounded projection."""
import argparse,gc,hashlib,json,os,shutil,subprocess
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix

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
 raw=np.bincount(mapping,weights=allnodes[:,4],minlength=N)/count
 bounded=np.bincount(mapping,weights=allnodes[:,5],minlength=N)/count
 den=np.asarray(W.sum(axis=0)).ravel();rhs=np.asarray(W.T@tp).ravel()
 residual=float(np.max(abs((rhs-W.T@(W@raw))/den)))
 assert residual<3e-12,residual
 lo=min(tp.min(),initial.min());hi=max(tp.max(),initial.max())
 gradient=(rhs-W.T@(W@bounded))/den
 projected=bounded-np.clip(bounded+gradient,lo,hi);projected[fixed]=0
 kkt=float(np.max(abs(projected)));assert kkt<3e-12,kkt
 assert np.max(abs(bounded[fixed]-initial[fixed]))<2e-12
 assert bounded.min()>=lo-2e-12 and bounded.max()<=hi+2e-12
 signed=np.concatenate([np.fromfile(d/f'probe.{rank}.signed.bin',dtype='=f8').reshape(-1,2) for rank in range(12)])
 signed_error=float(np.max(abs(signed[:,0]-signed[:,1])));assert signed_error<1e-8,signed_error
 assert signed[:,1].min()<0
 error=float(np.max(abs(bounded-initial)))
 if profile!='particle_jump':
  assert np.max(abs(W@initial-tp))<2e-12
  assert error<1e-8,error
 else:
  assert raw.min()<lo or raw.max()>hi
 energies=np.array([mass@initial,mass@raw,mass@bounded]);actual=np.loadtxt(d/'probe.energy.txt')
 assert np.max(abs(energies-actual[:3]))<3e-12
 return dict(global_nodes=N,particles=P,raw_residual=residual,bounded_KKT=kkt,signed_error=signed_error,grid_max_error_K=3400*error,net_change_percent=float(100*(energies[2]-energies[0])/energies[0]),minimum=float(bounded.min()),maximum=float(bounded.max()))

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('build',type=Path);p.add_argument('runs',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True)
 cases=[]
 for scenario in ['cold','hot','assim']:
  for variant in ['grid_coarse','base','grid_fine','np_low','np_high']:cases.append((f'{scenario}_{variant}',f'{scenario}_pices_{variant}','keep'))
 cases += [(mode,'cold_pices_base',mode) for mode in ['constant','cartesian_linear','particle_jump']]
 src=a.runs/'pices_p5';manifest={x['name']:x for x in json.loads((src/'matrix.json').read_text())['cases']};results={}
 for name,input_name,profile in cases:
  row=manifest[input_name];d=a.output/name;d.mkdir();(d/'DATA').mkdir()
  (d/'case.cfg').write_text((src/'cases'/f'{input_name}.cfg').read_text()+'\npices_projection=bounded_consistent\n');kind='assim' if row['scenario']=='assim' else 'transport';shutil.copy(src/f'refstate_{kind}_{row["nodes"]}.txt',d/'refstate.txt');shutil.copytree(src/f'forcing_{row["nodes"]}',d/'forcing',ignore=shutil.ignore_patterns('._*'))
  with (d/'stdout').open('w') as out,(d/'stderr').open('w') as err:
   run=subprocess.run([os.environ.get('MPIEXEC','/usr/local/bin/mpiexec'),'--oversubscribe','-n','12',str(a.build.resolve()/'PICESProjectionProbe'),'case.cfg'],cwd=d,env=dict(os.environ,PICES_PROBE_PROFILE=profile),stdout=out,stderr=err,timeout=180)
  assert run.returncode==0,(name,run.returncode,(d/'stderr').read_text()[-1600:])
  for rank in range(12):
   log=(d/f'DATA/{rank}/log').read_text();assert 'PICES_STEP ' not in log and 'PICES_TA ' not in log and 'PICES_STOKES ' not in log
  result=audit(d,profile);assert result['particles']==12*(row['nodes']-1)**3*row['particles']
  result['input_case']=input_name;result['profile']=profile;result['mpi_exit']=run.returncode;results[name]=result
  (a.output/'partial.json').write_text(json.dumps(results,indent=2)+'\n');print(name,'net energy %',result['net_change_percent'],'max error K',result['grid_max_error_K'],flush=True);gc.collect()
 report=dict(status='PASS',production_approved=False,cases=results,binary_sha256=hashlib.sha256((a.build/'PICESProjectionProbe').read_bytes()).hexdigest())
 (a.output/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
