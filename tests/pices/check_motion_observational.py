"""24-rank real solver regression: motion audit must not change temperatures/paths."""
import os,re,json,shutil,subprocess,argparse
from pathlib import Path
from check_motion_run import audit
p=argparse.ArgumentParser();p.add_argument('build',type=Path);p.add_argument('runs',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
out=a.output.resolve();out.mkdir();source=a.runs/'pices_p8c';exe=a.build.resolve()/'PICESManufactured'
for case in ('off','on'):
 d=out/case;d.mkdir();(d/'DATA').mkdir()
 s=(source/'continuous.cfg').read_text()
 for key,val in dict(nprocz=2,mgunitz=1,maxstep=2,maxtotstep=3,tracers_per_element=4,pices_projection='lumped',pices_ta_before_heat='on',fixed_timestep=1e-8).items():
  if re.search('^'+key+'=',s,re.M):s=re.sub('^'+key+'=.*$',key+'='+str(val),s,flags=re.M)
  else:s+='\n'+key+'='+str(val)+'\n'
 (d/'case.cfg').write_text(s);shutil.copy(source/'refstate.txt',d);shutil.copytree(source/'forcing',d/'forcing',ignore=shutil.ignore_patterns('._*'))
 env=dict(os.environ,PICES_TEST_MODE='rotation',PICES_MOTION_AUDIT='1' if case=='on' else '0')
 with (d/'stdout').open('w') as stdout,(d/'stderr').open('w') as stderr:
  rc=subprocess.run(['/usr/local/bin/mpiexec','--oversubscribe','-n','24',str(exe),'case.cfg'],cwd=d,env=env,stdout=stdout,stderr=stderr,timeout=600).returncode
 assert rc in (0,8),(case,rc,(d/'stderr').read_text()[-1500:])
 assert 'PICES_STEP step=2 ' in (d/'DATA/0/log').read_text(),case
result=audit(out/'on',2,24,False)
files=list((out/'off').glob('particles.*.2.txt'))+list((out/'off').glob('nodes.*.2.txt'))
assert len(files)==48
for f in files:assert f.read_bytes()==(out/'on'/f.name).read_bytes(),f
result['physical_snapshots_byte_identical']=len(files)
(out/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
