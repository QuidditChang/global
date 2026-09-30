"""P2 migration/restart failure tests using real production thermal objects."""
import argparse,gzip,hashlib,json,os,re,shutil,subprocess
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('build',type=Path);p.add_argument('cfg',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True)
base=a.cfg.read_text();results={}
def run(name,changes,mode='rotation',restart=None,expect=None):
 d=a.output/name;d.mkdir();(d/'DATA').mkdir();s=base
 opts=dict(datafile='PICES_P2',maxstep='32',maxtotstep='33',fixed_timestep='.02',tracers_per_element='32',checkpointFrequency='16',pices_checkpoint='on',pices_test_no_diffusion='on');opts.update(changes)
 if restart:opts.update(restart='on',solution_cycles_init='16',datafile_old='PICES_P2',datadir_old=str(restart.resolve()/'DATA/%RANK'))
 for k,v in opts.items():s=re.sub(r'^'+k+r'=.*$',k+'='+v,s,flags=re.M) if re.search(r'^'+k+'=',s,re.M) else s+'\n'+k+'='+v+'\n'
 (d/'case.cfg').write_text(s);ref=re.search(r'^refstate_file=(.+)$',s,re.M)[1];shutil.copy(a.cfg.parent/ref,d/ref)
 with (d/'stdout').open('w') as out,(d/'stderr').open('w') as err:
  r=subprocess.run([os.environ.get('MPIEXEC','/usr/local/bin/mpiexec'),'--oversubscribe','-n','12',str(a.build.resolve()/'PICESManufactured'),'case.cfg'],cwd=d,env=dict(os.environ,PICES_TEST_MODE=mode,PICES_TEST_IDS='1'),stdout=out,stderr=err,timeout=180)
 text=(d/'stderr').read_text()
 if expect:
  assert r.returncode==72 and expect in text,(name,r.returncode,text[-1500:])
  assert not any('PICES_STEP' in f.read_text() for f in (d/'DATA').glob('*/log'))
 else:
  assert r.returncode in (0,8),(name,r.returncode,text[-1500:])
  assert all('PICES_STEP step=32 ' in (d/'DATA'/str(i)/'log').read_text() for i in range(12))
 results[name]={'status':'PASS'};print(name,'PASS',flush=True);return d

def particles(d,step):
 rows=[]
 for rank in range(12):
  z=np.loadtxt(d/f'particles.{rank}.{step}.txt');rows.append(np.column_stack((z[:,0],z[:,-3:],np.full(len(z),rank))))
 v=np.concatenate(rows);return v[np.argsort(v[:,0])]

for mode in ['rotation','rotation_x','rotation_y']:
 d=run(mode,{ },mode)
 first,last=particles(d,0),particles(d,32)
 assert np.array_equal(first[:,0],last[:,0]),'Tp/ID loss or corruption'
 migrated=first[:,-1]!=last[:,-1];assert migrated.sum()>100 and len(set(first[migrated,-1]))==12
 # Compare to exact solid-body rotation as a diagnostic; coarse-grid interpolation
 # has spatial error, so this is not a time-order assertion.
 axis={'rotation':2,'rotation_x':0,'rotation_y':1}[mode];xyz=first[:,1:4];u=(axis+1)%3;v=(axis+2)%3;exact=xyz.copy();angle=32*np.float32(.02)
 exact[:,u]=xyz[:,u]*np.cos(angle)-xyz[:,v]*np.sin(angle);exact[:,v]=xyz[:,u]*np.sin(angle)+xyz[:,v]*np.cos(angle)
 err=np.linalg.norm(last[:,1:4]-exact,axis=1);assert np.isfinite(err).all() and err.max()<.1,err.max()
 results[mode].update(migrated=int(migrated.sum()),origin_ranks=12,max_position_error=float(err.max()))
 if mode=='rotation':source=d
r=run('rotation_restart',{},restart=source)
assert np.array_equal(particles(source,32),particles(r,32)),'restart changed Tp or coordinates'
results['rotation_restart']['bitwise_particle_snapshot_equal']=True
# Declared physics mismatch must reject even when the checkpoint is intact.
run('reject_physics',{'fixed_timestep':'.01'},restart=source,expect='metadata/physics')
run('reject_flavor',{'tracer_flavors':'1'},restart=source,expect='tracer schema')
for name,suffix,mode,token in [('missing_metadata','.pices.json','delete','missing'),('missing_manifest','.pices.manifest','delete','missing'),('corrupt_payload','','flip','checksum'),('corrupt_velocity','.pices.state','flip','checksum'),('wrong_slot','.pices.json','slot','metadata/physics')]:
 copy=a.output/(name+'_input');shutil.copytree(source/'DATA',copy/'DATA')
 path=copy/'DATA/0'/('PICES_P2.chkpt.0.16'+suffix)
 if mode=='delete':path.unlink()
 elif mode=='flip':v=bytearray(path.read_bytes());v[-1]^=1;path.write_bytes(v)
 else:
  path.write_text(path.read_text().replace('"slot":0','"slot":1'))
  # Re-sign manifests so slot validation, not merely the collective hash, is tested.
  hashes=b''.join(hashlib.sha256((copy/f'DATA/{i}/PICES_P2.chkpt.{i}.16.pices.json').read_bytes()).hexdigest().encode()+b'\0' for i in range(12))
  digest=hashlib.sha256(hashes).hexdigest()
  for i in range(12):(copy/f'DATA/{i}/PICES_P2.chkpt.{i}.16.pices.manifest').write_text('{"magic":"CITCOMS_EBA_PICES_COMPLETE","schema":1,"mpi_size":12,"metadata_set_sha256":"'+digest+'"}\n')
 run('reject_'+name,{},restart=copy,expect=token)
(a.output/'summary.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(results,indent=2))
