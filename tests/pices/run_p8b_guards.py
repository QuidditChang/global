"""Rejection guards, frozen-viscosity plates and an active particle-CFL limit."""
import argparse,gzip,json,math,os,shutil,subprocess
from pathlib import Path
from prepare_p8a import edit
from verify_p2 import fields
p=argparse.ArgumentParser();p.add_argument('build',type=Path);p.add_argument('runs',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True)
src=a.runs/'pices_p8b';results={}
def run(name,changes,plate_scale=1):
 d=a.output/name;d.mkdir();(d/'DATA').mkdir()
 (d/'case.cfg').write_text(edit((src/'coupled.cfg').read_text(),changes))
 shutil.copy(src/'refstate.txt',d);shutil.copytree(src/'forcing',d/'forcing',ignore=shutil.ignore_patterns('._*'))
 if plate_scale!=1:
  for f in (d/'forcing').glob('bvel.*'):
   f.write_text(''.join(' '.join(str(float(v)*plate_scale) for v in line.split())+'\n' for line in f.read_text().splitlines()))
 with (d/'stdout').open('w') as out,(d/'stderr').open('w') as err:
  r=subprocess.run([os.environ.get('MPIEXEC','/usr/local/bin/mpiexec'),'--oversubscribe','-n','12',str(a.build.resolve()/'CitcomSFull'),'case.cfg'],cwd=d,stdout=out,stderr=err,timeout=900)
 return d,r.returncode
for name,changes,message in [('checkpoint',dict(pices_checkpoint='on',p5_case='off'),'requires P8c'),('rotation',dict(remove_rigid_rotation=1),'remove_rigid_rotation=off'),('invalid_limit',dict(pices_max_timestep_Ma=0),'must be finite and positive')]:
 d,code=run(name,changes)
 assert code not in (0,8) and message in (d/'stderr').read_text(),name
 assert not any('PICES_STEP ' in f.read_text() for f in (d/'DATA').glob('*/log')),name
 results[name]='PASS'
scalev=6371000/(4/(3300*1250))/(100*365.25*86400)
scalet=6371000**2/(4/(3300*1250))/(1e6*365.25*86400)
for name,changes,plate_scale in [('frozen_viscosity',dict(VISC_UPDATE='off'),1),('particle_cfl',{},10000)]:
 d,code=run(name,dict(changes,maxstep=1,maxtotstep=2),plate_scale);assert code in (0,8),name
 for rank in range(12):
  folder=d/f'DATA/{rank}';log=(folder/'log').read_text().splitlines()
  steps=[fields(s) for s in log if s.startswith('PICES_STEP ')];assert len(steps)==1 and steps[0]['step']=='1',name
  age=2.13-float(steps[0]['time'])*scalet
  data=gzip.open(folder/f'1/velo.{rank}.1.gz','rt').read().splitlines()[2:]
  assert all(math.isclose(float(line.split()[1]),plate_scale*.01*(1+age)*scalev,rel_tol=2e-6) for line in data[8::9]),name
  if name=='particle_cfl':
   dt=[fields(s) for s in log if s.startswith('PICES_TIMESTEP ')][0]
   assert math.isclose(float(dt['dt']),float(dt['particle_limit']),rel_tol=2e-7),dt
   assert .239999<float(steps[0]['cfl'])<=.240001,steps[0]
 results[name]='PASS'
(a.output/'summary.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(results))
