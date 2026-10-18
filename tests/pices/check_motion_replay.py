"""Exercise real restart, frozen-field subdivisions and cross-rank migration.
Uses the observational test's checkpoint zero and an amplified prescribed
rotation field solely to make transfer measurable on a tiny local mesh.
"""
import argparse,os,re,json,struct,shutil,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('build',type=Path);p.add_argument('source',type=Path);p.add_argument('output',type=Path);p.add_argument('--steps',type=int,default=1);p.add_argument('--single-reference',type=Path);a=p.parse_args();source=a.source.resolve();out=a.output.resolve();out.mkdir()
for rank in range(24):
 d=source/'DATA'/str(rank);meta=json.loads((d/f'PICES_P8c.chkpt.{rank}.0.pices.json').read_text());n=1
 for v in meta['local_mesh']:n*=v
 raw=(d/f'PICES_P8c.chkpt.{rank}.0.pices.state').read_bytes()[:3*n*4]
 values=struct.unpack('='+str(3*n)+'f',raw)
 # Native-endian fixture. The production reader rejects inconsistent clocks.
 clock=meta['time'];dt=struct.unpack('=f',struct.pack('=f',1e-7))[0]
 for step in range(1,a.steps+1):
  (d/f'global.motion_velocity.{rank}.{step}').write_bytes(struct.pack('=4i2d',0x504d5631,rank,step,n,clock,dt)+struct.pack('='+str(3*n)+'f',*[v*10000 for v in values]))
  clock=struct.unpack('=f',struct.pack('=f',clock+dt))[0]
positions=[];migrations=[]
for factor in (1,2,4):
 d=out/str(factor);d.mkdir();(d/'DATA').mkdir();shutil.copy(source/'refstate.txt',d);shutil.copytree(source/'forcing',d/'forcing')
 s=(source/'case.cfg').read_text()
 for k,v in dict(restart='on',solution_cycles_init=0,datadir_old=str(source/'DATA/%RANK'),datafile_old='PICES_P8c').items():
  if re.search('^'+k+'=',s,re.M):s=re.sub('^'+k+'=.*$',k+'='+str(v),s,flags=re.M)
  else:s+='\n'+k+'='+str(v)+'\n'
 (d/'case.cfg').write_text(s)
 env=dict(os.environ,PICES_MOTION_AUDIT='1',PICES_REPLAY_ROOT=str(source),PICES_REPLAY_SUBSTEPS=str(factor))
 if a.steps>1:env['PICES_REPLAY_END_STEP']=str(a.steps)
 with (d/'stdout').open('w') as stdout,(d/'stderr').open('w') as stderr:
  rc=subprocess.run(['/usr/local/bin/mpiexec','--oversubscribe','-n','24',str(a.build.resolve()/'PICESMotionReplay'),'case.cfg'],cwd=d,env=env,stdout=stdout,stderr=stderr,timeout=600).returncode
 assert rc==0,(factor,rc,(d/'stderr').read_text()[-1500:])
 state={};moved=0
 for f in d.glob('DATA/*/*.replay.*'+('.'+str(a.steps) if a.steps>1 else '')):
  rank=int(f.parent.name)
  for line in f.read_text().splitlines():
   if line.startswith('#'):continue
   uid,x,y,z,host=line.split();uid=int(uid);assert uid not in state;state[uid]=tuple(map(float,(x,y,z)));moved+=uid//2**32!=rank
 assert len(state)==24576,len(state)
 positions.append(state);migrations.append(moved)
assert positions[0].keys()==positions[1].keys()==positions[2].keys()
errors=[max(sum((positions[i][k][d]-positions[2][k][d])**2 for d in range(3))**.5 for k in positions[2]) for i in (0,1)]
assert all(n>0 for n in migrations),migrations
assert errors[1]<errors[0],errors
if a.single_reference:
 for factor in (1,2,4):
  for f in (out/str(factor)).glob('DATA/*/*.replay.*.1'):
   original=a.single_reference/str(factor)/f.relative_to(out/str(factor)).with_name(f.name.rsplit('.',1)[0])
   assert f.read_bytes()==original.read_bytes(),(f,original)
result=dict(status='PASS',steps=a.steps,particles=24576,subdivisions=[1,2,4],cross_rank_endpoints=migrations,max_endpoint_difference_to_four=errors)
(out/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
