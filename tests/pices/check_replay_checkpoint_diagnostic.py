"""Real MPI positive/negative checks against a freshly generated checkpoint.
The build must contain PICESMotionReplay, PICESMotionReplay_plain and
PICESMotionReplay_legacy_probe. Only /tmp test source checkpoint is read.
"""
import argparse,json,os,re,shutil,struct,subprocess
from pathlib import Path
from audit_replay_checkpoint import audit
from motion_window import parameters,replace

def main():
 p=argparse.ArgumentParser(description=__doc__)
 for key in ('build','source','output'):p.add_argument(key,type=Path)
 a=p.parse_args();source=a.source.resolve();out=a.output.resolve();out.mkdir()
 for rank in range(24):
  d=source/'DATA'/str(rank);meta=json.loads((d/('PICES_P8c.chkpt.%d.0.pices.json'%rank)).read_text());n=1
  for v in meta['local_mesh']:n*=v
  raw=(d/('PICES_P8c.chkpt.%d.0.pices.state'%rank)).read_bytes()[:12*n]
  dt=struct.unpack('=f',struct.pack('=f',1e-7))[0]
  (d/('global.motion_velocity.%d.1'%rank)).write_bytes(struct.pack('=4i2d',0x504d5631,rank,1,n,meta['time'],dt)+raw)
 cases=[('plain','PICESMotionReplay_plain'),('exact','PICESMotionReplay'),('physics_mismatch','PICESMotionReplay'),('unused_fields','PICESMotionReplay_legacy_probe'),('checkpoint_only','PICESMotionReplay')]
 for name,exe in cases:
  d=out/name;d.mkdir();(d/'DATA').mkdir();shutil.copy(source/'refstate.txt',d);shutil.copytree(source/'forcing',d/'forcing')
  text=(source/'case.cfg').read_text()
  updates=dict(restart='on',solution_cycles_init=0,datadir_old=str(source/'DATA/%RANK'),datafile_old='PICES_P8c',datafile='global',datadir='DATA/%RANK')
  if name=='physics_mismatch':updates['qvis_friction_angle_rad']='0.086'
  for key,value in updates.items():text=replace(text,key,value) if key in parameters(text) else text+'\n'+key+'='+str(value)+'\n'
  (d/'case.cfg').write_text(text)
  env=dict(os.environ,PICES_MOTION_AUDIT='1',PICES_REPLAY_ROOT=str(source),PICES_REPLAY_SUBSTEPS='1')
  for key in ('PICES_REPLAY_END_STEP','PICES_REPLAY_IDS','PICES_REPLAY_CHECKPOINT_ONLY','PICES_TEST_MODE','PICES_RESTART_SOURCE_COMMIT'):env.pop(key,None)
  if name=='checkpoint_only':env['PICES_REPLAY_CHECKPOINT_ONLY']='1'
  with (d/'stdout').open('w') as stdout,(d/'stderr').open('w') as stderr:
   rc=subprocess.run(['/usr/local/bin/mpiexec','--oversubscribe','-n','24',str(a.build.resolve()/exe),'case.cfg'],cwd=d,env=env,stdout=stdout,stderr=stderr,timeout=300).returncode
  if name=='plain':assert rc==0,(rc,(d/'stderr').read_text()[-1000:]);continue
  result=audit(d,24);(d/'audit.json').write_text(json.dumps(result,indent=2)+'\n');assert not result['missing_ranks'],result
  files=list(d.glob('DATA/*/*.replay.*'))
  if name in ('physics_mismatch','unused_fields'):
   assert rc!=0 and not files and result['field_counts']=={'physics_mesh_sha256':24},(rc,result)
   if name=='unused_fields':assert result['zero_unused_candidate_matches_saved_physics']==list(range(24)),result
   else:assert result['zero_unused_candidate_matches_saved_physics']==[],result
  else:
   assert rc==0 and result['status']=='EXACT',(rc,result)
   if name=='exact':
    assert len(files)==24
    for f in files:assert f.read_bytes()==(out/'plain'/f.relative_to(d)).read_bytes(),f
   else:assert not files and 'motion_steps=0' in (d/'stderr').read_text()
  print(name,'PASS',flush=True)
 print('PASS: exact metadata and unchanged endpoints; physical and obsolete-field differences rejected with all-rank reports; checkpoint-only stops without motion')
if __name__=='__main__':main()
