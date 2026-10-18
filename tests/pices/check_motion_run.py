"""Audit motion signatures, replay inventory and sparse path pairing; no physics changes."""
import argparse,json,re,struct
from pathlib import Path

def audit(root,steps=300,ranks=384,replay=True):
 root=Path(root);text=(root/'DATA/0/global.log').read_text() if (root/'DATA/0/global.log').exists() else (root/'DATA/0/log').read_text()
 rows=[dict(re.findall(r'(\w+)=([^ ]+)',line)) for line in text.splitlines() if line.startswith('PICES_MOTION ')]
 assert len(rows)==2*steps,('motion rows',len(rows),steps)
 assert [(int(r['step']),int(r['stage'])) for r in rows]==[(s,k) for s in range(1,steps+1) for k in (1,2)]
 assert len({(r['count'],r['hash_sum'],r['hash_xor']) for r in rows})==1,'ID signatures changed'
 files=sorted(root.glob('DATA/*/*.motion.*'))
 assert len(files)==ranks,('motion log count',len(files))
 # Disk-backed path join avoids retaining millions of trajectory IDs in RAM.
 import sqlite3,tempfile
 totals={'paths':0,'cells':0,'valid_gradients':0,'gradient_unavailable':0}
 with tempfile.TemporaryDirectory() as d:
  db=sqlite3.connect(str(Path(d)/'paths.db'));db.execute('create table paths(step int,id text,stage int,primary key(step,id,stage))')
  for f in files:
   for line in f.open():
    v=line.split()
    if not v:continue
    if v[0]=='CELL':totals['cells']+=1
    if v[0]=='HOST':assert int(v[5])==1,('invalid host',str(f),line)
    if v[0]!='PATH':continue
    step,stage=int(v[1]),int(v[2]);uid=v[3]
    db.execute('insert into paths values(?,?,?)',(step,uid,stage));totals['paths']+=1
    if stage==1:totals['valid_gradients' if int(v[18]) else 'gradient_unavailable']+=1
   db.commit()
  assert db.execute('select count(*) from (select step,id from paths group by step,id having count(*)!=2)').fetchone()[0]==0,'unpaired source/destination paths'
 inventory={}
 if replay:
  for rank in range(ranks):
   d=root/'DATA'/str(rank);mesh=d/f'global.motion_mesh.{rank}'
   with mesh.open('rb') as f:
    magic,rr,cap,nno,nel,elx,elz=struct.unpack('=7i',f.read(28))
   assert magic==0x504d4d31 and rr==rank
   assert mesh.stat().st_size==28+40+3*nno*8+nel*8*4
   for step in range(81,106):
    v=d/f'global.motion_velocity.{rank}.{step}'
    with v.open('rb') as f:head=struct.unpack('=4i',f.read(16))
    assert head==(0x504d5631,rank,step,nno)
    assert v.stat().st_size==32+3*nno*4
   for step in (80,90,100):
    stem=d/f'global.chkpt.{rank}.{step}'
    for suffix in ('','.pices.json','.pices.state','.pices.manifest'):assert Path(str(stem)+suffix).stat().st_size>0
    metadata=json.loads(Path(str(stem)+'.pices.json').read_text())
    assert {'motion_id','motion_selected'}<={q['name'] for q in metadata['extraq']}
  inventory={'ranks':ranks,'mesh_files':ranks,'velocity_files':25*ranks,'checkpoint_steps':[80,90,100]}
 return dict(status='PASS',steps=steps,identity_signatures='constant_count_and_two_hashes_not_exact_uniqueness_proof',
             clipped=sum(int(x['clipped']) for x in rows),limited=sum(int(x['limited']) for x in rows),
             sampled_div_max_nd=max(float(x['sampled_div_max_nd']) for x in rows),sparse=totals,replay_inventory=inventory,
             interpretation='observational capture validated; frozen-field convergence and alternative MPI decomposition still require replay')
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path);p.add_argument('--steps',type=int,default=300);p.add_argument('--ranks',type=int,default=384);p.add_argument('--no-replay-inventory',action='store_true');a=p.parse_args()
 result=audit(a.root,a.steps,a.ranks,not a.no_replay_inventory);(a.root/'motion_audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
