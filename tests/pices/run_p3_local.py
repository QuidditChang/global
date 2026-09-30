"""P3 source solution, nonlinear steady convergence and fail-closed integration."""
import argparse,json,os,re,subprocess
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('build',type=Path);p.add_argument('cfg',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True)
base=a.cfg.read_text();results={}
def run(name,changes,mode='constant',expect=None,variable=False):
 d=a.output/name;d.mkdir();(d/'DATA').mkdir();s=base
 opts=dict(pices_eba='on',maxstep='1',maxtotstep='2',tracers_per_element='32',pices_checkpoint='off',Q0='0',dissipation_number='0',fixed_timestep='.001');opts.update(changes)
 for k,v in opts.items():
  v=str(v);s=re.sub(r'^'+k+'=.*$',k+'='+v,s,flags=re.M) if re.search(r'^'+k+'=',s,re.M) else s+'\n'+k+'='+v+'\n'
 (d/'case.cfg').write_text(s);ref=re.search(r'^refstate_file=(.+)$',s,re.M)[1];nz=int(opts.get('nodez',5))
 (d/ref).write_text('\n'.join(f'{1+.1*(nz-1-i)/(nz-1) if variable else 1} 1 {1-i/(nz-1)} 1 {1+.2*(nz-1-i)/(nz-1) if variable else 1}' for i in range(nz))+'\n')
 with (d/'stdout').open('w') as out,(d/'stderr').open('w') as err:
  r=subprocess.run([os.environ.get('MPIEXEC','/usr/local/bin/mpiexec'),'--oversubscribe','-n','12',str(a.build.resolve()/'PICESManufactured'),'case.cfg'],cwd=d,env=dict(os.environ,PICES_TEST_MODE=mode),stdout=out,stderr=err,timeout=180)
 if expect:assert r.returncode==72 and expect in (d/'stderr').read_text(),(name,r.returncode,(d/'stderr').read_text()[-2000:])
 else:
  assert r.returncode in (0,8),(name,r.returncode,(d/'stderr').read_text()[-2000:])
  assert all('PICES_STEP step=1 ' in (d/f'DATA/{rank}/log').read_text() for rank in range(12))
 results[name]={'status':'PASS'};print(name,'PASS',flush=True);return d
# A uniform specific heat source with no diffusion has dT/dt=Q/Cp=1.
d=run('source_exact',{'Q0':'1','fixed_timestep':'.025','pices_test_no_diffusion':'on'})
err=0
for rank in range(12):
 rows=np.loadtxt(d/f'nodes.{rank}.1.txt');free=rows[:,5]==0
 err=max(err,float(np.max(abs(rows[free,4]-(.5+float(np.float32(.025)))))))
 assert int(re.search(r'PICES_STEP .*substeps=(\d+)',(d/f'DATA/{rank}/log').read_text())[1])>=3
assert err<1e-13,err
results['source_exact']['max_temperature_error']=err
# The exact spherical steady solution integrates k(T)dT = const*dr/r^2.
# Check the full remap + heat update with independently known continuum T(r).
errors=[]
for nz in [5,9]:
 d=run('nonlinear_steady_'+str(nz),{'nodex':nz,'nodey':nz,'nodez':nz,'mgunitx':nz-1,'mgunity':nz-1,'mgunitz':nz-1,'kT_exponent':'.3'},'nonlinear_steady',variable=True)
 sq=[]
 for rank in range(12):
  initial=np.loadtxt(d/f'nodes.{rank}.0.txt');final=np.loadtxt(d/f'nodes.{rank}.1.txt');free=final[:,5]==0
  sq.extend((final[free,4]-initial[free,4])**2)
 errors.append(float(np.sqrt(np.mean(sq))))
assert errors[1]<errors[0]/1.5,errors
results['nonlinear_steady_convergence']={'status':'PASS','rms_errors':errors,'ratio':errors[0]/errors[1]}
run('reject_capacity',{'phase_delta_s':'100,100,100','phase_clapeyron':'1,1,1','phase_width':'1,1,1'},expect='effective phase capacity')
run('reject_source_limit',{'Q0':'1','fixed_timestep':'.025','pices_test_no_diffusion':'on','pices_max_substeps':'1'},expect='heat substep limit')
run('reject_composition',{'kC_ratio':'2'},expect='P3 requires EBA')
(a.output/'summary.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(results,indent=2))
