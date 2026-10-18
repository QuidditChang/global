"""Real MPI regression for Pyre unset values passed to the standalone replay.
Requires the checkpoint/velocity fixture used by check_motion_replay.py.
"""
import argparse,os,re,shutil,subprocess
from pathlib import Path
from motion_window import standalone_runtime,standalone_mesh,restore_input_precision,parameters,replace
p=argparse.ArgumentParser(description=__doc__);p.add_argument('build',type=Path);p.add_argument('single_case',type=Path);p.add_argument('velocity_root',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
source=a.single_case.resolve();out=a.output.resolve();out.mkdir()
original=(source/'case.cfg').read_text()
# Match production input.cfg array whitespace, absent from the old fixture.
for key,value in parameters(original).items():
 if ',' in value:
  try:[float(v) for v in value.split(',')]
  except ValueError:continue
  original=replace(original,key,', '.join(v.strip() for v in value.split(',')))
raw=re.sub(r'(?m)^mgunit[xyz]=.*$', '',original)+'\ninputdiffusivity=-1\nreference_conductivity=-1\n'
fixed,removed=standalone_runtime(raw);assert set(removed)=={'inputdiffusivity','reference_conductivity'}
fixed,restored=restore_input_precision(fixed,original)
missing_mesh=fixed
fixed,units,mesh=standalone_mesh(fixed)
env=dict(os.environ,PICES_MOTION_AUDIT='1',PICES_REPLAY_ROOT=str(a.velocity_root.resolve()),PICES_REPLAY_SUBSTEPS='1')
for key in ('PICES_REPLAY_IDS','PICES_REPLAY_END_STEP'):env.pop(key,None)
for name,text in (('serialized_pyre',raw),('missing_mesh',missing_mesh),('standalone',fixed)):
 d=out/name;d.mkdir();(d/'DATA').mkdir();shutil.copy(source/'refstate.txt',d);shutil.copytree(source/'forcing',d/'forcing');(d/'case.cfg').write_text(text)
 with (d/'stdout').open('w') as stdout,(d/'stderr').open('w') as stderr:
  rc=subprocess.run(['/usr/local/bin/mpiexec','--oversubscribe','-n','24',str(a.build.resolve()/'PICESMotionReplay'),'case.cfg'],cwd=d,env=env,stdout=stdout,stderr=stderr,timeout=600).returncode
 if name=='serialized_pyre':
  assert rc!=0 and 'remove the deprecated independent value' in (d/'stderr').read_text(),rc
  assert not list(d.glob('DATA/*/*.replay.*'))
 elif name=='missing_mesh':
  assert rc!=0 and 'accepted velocity length mismatch' in (d/'stderr').read_text(),rc
  assert not list(d.glob('DATA/*/*.replay.*'))
 else:
  assert rc==0,(rc,(d/'stderr').read_text()[-1500:])
  files=list(d.glob('DATA/*/*.replay.*'));assert len(files)==24
  for f in files:assert f.read_bytes()==(source/f.relative_to(d)).read_bytes(),f
print('PASS: reproduced conductivity and checkpoint-length rejections; translated cfg restarted and produced 24 byte-identical particle endpoint files')
