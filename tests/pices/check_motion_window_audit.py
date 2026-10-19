"""Exercise the audit on real multi-step MPI replay outputs, then inject bad evidence."""
import argparse,contextlib,io,json,os,tempfile
from pathlib import Path
from motion_window import audit,fields,read_positions
p=argparse.ArgumentParser();p.add_argument('replay',type=Path);a=p.parse_args()
source=a.replay.resolve();root=Path(tempfile.mkdtemp(prefix='pices-window-audit-'));steps=range(1,4);ranks=24
for factor in (1,2,4):
 case=root/('sub%d'%factor);case.mkdir();(case/'returncode.txt').write_text('0\n')
 for rank in range(ranks):
  d=case/'DATA'/str(rank);d.mkdir(parents=True)
  for f in (source/str(factor)/'DATA'/str(rank)).glob('*.replay*'):
   (d/f.name.replace('PICES_P8c','global')).symlink_to(f)
  if rank==0:(d/'global.log').symlink_to(source/str(factor)/'DATA/0/log')
lines=(root/'sub1/DATA/0/global.log').read_text().splitlines()
signature=next([f[k] for k in ('count','hash_sum','hash_xor')] for l in lines if l.startswith('PICES_MOTION ') for f in [fields(l)])
hist={int(f['step']):list(map(int,f['histogram'].split(','))) for l in lines if l.startswith('PICES_REPLAY ') for f in [fields(l)]}
low={s:[] for s in steps};references=[];ids=None
for step in steps:
 state=read_positions(root/'sub1',step,ranks);ids=set(state) if ids is None else ids;assert ids==state.keys()
 for uid,(xyz,rank,host) in state.items():references.append('%d %d %.17g %.17g %.17g %d %d\n'%((step,uid)+xyz+(rank,host)))
 for rank in range(ranks):
  for line in (root/'sub1/DATA'/str(rank)/('global.replay_cells.%d.%d'%(rank,step))).open():
   e,n=map(int,line.split())
   if n<=4:low[step].append([rank,e,n])
(root/'reference.tsv').write_text(''.join(references));(root/'tracked_ids.txt').write_text(''.join(str(i)+'\n' for i in sorted(ids)))
(root/'replay_manifest.json').write_text(json.dumps(dict(first=1,last=3,ranks=ranks,signature=signature,histograms=hist,low_cells=low,scope='manufactured integration regression')))
with contextlib.redirect_stdout(io.StringIO()):audit(root)
assert json.loads((root/'replay_audit.json').read_text())['status']=='PASS_REVIEW_REQUIRED'
# The baseline gate must work before either comparison run exists.
for factor in (2,4):(root/('sub%d'%factor)).rename(root/('pending%d'%factor))
with contextlib.redirect_stdout(io.StringIO()):audit(root,baseline_only=True)
assert json.loads((root/'baseline_audit.json').read_text())['status']=='PASS_REVIEW_REQUIRED'
original=(root/'reference.tsv').read_text();v=original.splitlines();first=v[0].split();first[2]=str(float(first[2])+1e-5);v[0]=' '.join(first);(root/'reference.tsv').write_text('\n'.join(v)+'\n')
try:
 with contextlib.redirect_stdout(io.StringIO()):audit(root,baseline_only=True)
except ValueError:pass
else:raise AssertionError('altered reference position was accepted')
assert json.loads((root/'baseline_audit.json').read_text())['status']=='FAIL_BASELINE'
(root/'reference.tsv').write_text(original)
m=json.loads((root/'replay_manifest.json').read_text());m['low_cells']['1'][0][2]+=1;(root/'replay_manifest.json').write_text(json.dumps(m))
try:
 with contextlib.redirect_stdout(io.StringIO()):audit(root,baseline_only=True)
except ValueError:pass
else:raise AssertionError('altered source cell census was accepted')
print('PASS: real MPI outputs accepted; altered endpoint and low-cell map rejected')
print(root)
