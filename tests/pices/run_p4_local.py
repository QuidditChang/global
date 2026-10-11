"""P4 full-kernel tests: subcycling, particle reaction, separate ledgers, CBF."""
import argparse,json,math,os,re,subprocess
from pathlib import Path
import numpy as np
from verify_p2 import fields
from verify_p4 import boundary
p=argparse.ArgumentParser();p.add_argument('build',type=Path);p.add_argument('cfg',type=Path);p.add_argument('--output',type=Path,required=True);p.add_argument('--cold',action='store_true');a=p.parse_args();a.output.mkdir(parents=True)
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
  r=subprocess.run([os.environ.get('MPIEXEC','/usr/local/bin/mpiexec'),'--oversubscribe','-n','12',str(a.build.resolve()/'PICESManufactured'),'case.cfg'],cwd=d,env=dict(os.environ,PICES_TEST_MODE='p4_linear' if resolved else ('p4_cold' if a.cold else 'constant')),stdout=out,stderr=err,timeout=180)
 assert r.returncode in (0,8),(name,r.returncode,(d/'stderr').read_text()[-1000:])
 assert all('PICES_STEP step=1 ' in (d/f'DATA/{rank}/log').read_text() for rank in range(12)),name
 print(name,'completed',flush=True);return d

heat=run('heat_only',ta=False);ta=run('heat_and_TA');particle_error=0.;energy_delta={};energy_ta={}
for rank in range(12):
 before=np.loadtxt(heat/f'nodes.{rank}.1.txt');after=np.loadtxt(ta/f'nodes.{rank}.1.txt');data=np.loadtxt(ta/f'ta.{rank}.1.txt')
 h={k:float(v) for k,v in fields((ta/f'ta.{rank}.1.txt').read_text().splitlines()[0]).items()}
 alpha=np.zeros(len(data));target=np.zeros(len(data));delta=after[:,4]-before[:,4]
 for n,(_,inc,depth,tref,age,flag) in enumerate(data):
  if before[n,5] or depth>=h['depth']:continue
  age=min(age,h['cap']/h['scalet']);target[n]=tref-h['surface']*math.erfc(.5*depth/math.sqrt(age))
  w=math.expm1(-h['exp']*depth/h['depth'])/math.expm1(-h['exp'])
  physical_dt=float(np.float32(np.float32(h['dt'])*np.float32(h['scalet'])))
  alpha[n]=-math.expm1(-physical_dt*(1-w)/h['tau'])
 oldp=np.loadtxt(heat/f'particles.{rank}.1.txt');newp=np.loadtxt(ta/f'particles.{rank}.1.txt')
 assert np.array_equal(oldp[:,1:],newp[:,1:]),'TA moved particles'
 ids=newp[:,1:17:2].astype(int)-1;weights=newp[:,2:17:2]
 expected=(1-(weights*alpha[ids]).sum(axis=1))*oldp[:,0]+(weights*alpha[ids]*target[ids]).sum(axis=1)
 particle_error=max(particle_error,float(np.max(abs(newp[:,0]-expected))))
 assert np.min(300+3400*newp[:,0])>=0
 # Node masses already include MPI contributions: count shared nodes once.
 for xyz,d,t,mass in zip(after[:,1:4],delta,data[:,1],after[:,6]):
  key=tuple(np.round(xyz,10));energy_delta[key]=d*mass;energy_ta[key]=t*mass
 log=(ta/f'DATA/{rank}/log').read_text();assert log.count('PICES_TA ')==1
 assert int(fields(next(l for l in log.splitlines() if l.startswith('PICES_STEP ')))['substeps'])>=3
 for side in ['surf','botm']:assert (heat/f'DATA/{rank}/q.{side}.{rank}.1').read_bytes()==(ta/f'DATA/{rank}/q.{side}.{rank}.1').read_bytes()
assert particle_error<2e-14,particle_error
rows=(ta/'DATA/0/log').read_text().splitlines()
ta_ledger=fields(next(x for x in rows if x.startswith('PICES_TA ')))
assert ta_ledger['operator']=='particle_relaxation_v1'
assert abs(sum(energy_ta.values())-float(ta_ledger['storage']))<1e-12
assert abs(sum(energy_delta.values())-float(ta_ledger['total_grid_storage']))<1e-12
ledger=fields(next(x for x in rows if x.startswith('PICES_EBA ')))
reaction=float(ledger['boundary_reaction'])/h['dt']*4*3400*6371000
flux=boundary(ta,1,'botm')-boundary(ta,1,'surf');balance=abs(flux-reaction)/max(1,abs(flux)+abs(reaction));assert balance<1e-12,balance
results['multisubstep_TA_and_CBF']=dict(status='PASS',particle_reaction_error=particle_error,
 TA_storage_error=abs(sum(energy_ta.values())-float(ta_ledger['storage'])),
 total_storage_error=abs(sum(energy_delta.values())-float(ta_ledger['total_grid_storage'])),
 CBF_before_after_TA_equal=True,relative_CBF_boundary_residual=balance)
(a.output/'summary.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(results,indent=2))
