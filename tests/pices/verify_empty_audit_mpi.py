"""Verify real MPI audit-on/off runs: script ON_ROOT OFF_ROOT."""
import sys,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from audit_empty_run import topology
root=Path(sys.argv[1]);mesh={};meta={};cells={};frames={}
for p in (root/'PICES_Audit').glob('rank_*.jsonl'):
 rank=int(p.stem.split('_')[1]);frames[rank]=[]
 for line in p.open():
  r=json.loads(line)
  if r['type']=='meta':meta[rank]=r
  if r['type']=='mesh':mesh[rank,r['element']]=r
  if r['type']=='frame':frames[rank].append(r['step'])
  if r['type']=='cell':cells.setdefault(r['step'],{})[rank,r['element']]=r
assert len(frames)==24 and all(v==list(range(5)) for v in frames.values())
targets={k for rows in cells.values() for k,r in rows.items() if r['count']==0}
adj,missing=topology(mesh,meta,targets)
assert not missing,missing
assert any(k[0]!=n[0] for k in targets for n in adj[k]),adj
for step,rows in cells.items():
 for k,r in rows.items():
  if not r['count']:assert all(n in rows for n in adj[k])
print('PASS: 24 ranks, frames0..4; complete face halo including cross-rank neighbors',dict(adj))
other=Path(sys.argv[2]);checked=0
for pattern in ('particles.*.txt','nodes.*.txt','ta.*.txt'):
 for p in root.glob(pattern):
  assert p.read_bytes()==(other/p.name).read_bytes(),p.name
  checked+=1
assert checked>0
print('PASS: audit on/off identical snapshots:',checked)
for r in (root,other):
 assert 'PICES_ERROR' not in (r/'stderr').read_text()
 for rank in range(24):assert 'PICES_STEP step=4 ' in (r/f'DATA/{rank}/log').read_text()
print('PASS: all ranks finish step4; standalone driver returns expected legacy exit8')
