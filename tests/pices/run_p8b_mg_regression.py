"""Compare PG CG/MG initial fields to cover the shared V-cycle correction."""
import argparse,gzip,json,math,os,shutil,subprocess
from pathlib import Path
from prepare_p8a import edit
p=argparse.ArgumentParser();p.add_argument('build',type=Path);p.add_argument('runs',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True)
src=a.runs/'pices_p8b';values={}
for name in ('coupled','multigrid'):
 d=a.output/name;d.mkdir();(d/'DATA').mkdir()
 (d/'case.cfg').write_text(edit((src/(name+'.cfg')).read_text(),dict(energy_solver='pg',pices_projection='lumped',pices_eba='off',pices_p4='off',CBF_use_advection='on',maxstep=1,maxtotstep=2,fixed_timestep=1e-9)))
 shutil.copy(src/'refstate.txt',d);shutil.copytree(src/'forcing',d/'forcing',ignore=shutil.ignore_patterns('._*'))
 with (d/'stdout').open('w') as out,(d/'stderr').open('w') as err:
  r=subprocess.run([os.environ.get('MPIEXEC','/usr/local/bin/mpiexec'),'--oversubscribe','-n','12',str(a.build.resolve()/'CitcomSFull'),'case.cfg'],cwd=d,stdout=out,stderr=err,timeout=3600)
 assert r.returncode in (0,8),name
 v=[]
 for rank in range(12):
  folder=d/f'DATA/{rank}'
  assert (folder/f'1/velo.{rank}.1.gz').exists(),name
  lines=gzip.open(folder/f'0/velo.{rank}.0.gz','rt').read().splitlines()[2:]
  v.extend(float(x) for line in lines for x in line.split()[:3])
 assert all(map(math.isfinite,v)),name
 values[name]=v
rel=math.sqrt(sum((x-y)**2 for x,y in zip(values['coupled'],values['multigrid']))/sum(x*x for x in values['coupled']))
assert rel<1e-4,rel
result=dict(status='PASS',PG_initial_velocity_CG_MG_relative_L2=rel,tolerance=1e-4)
(a.output/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
