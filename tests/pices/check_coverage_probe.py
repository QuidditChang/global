"""Validate prescribed-rotation probe on/off runs (same cfg, binary and MPI layout)."""
from pathlib import Path
import argparse,re
p=argparse.ArgumentParser();p.add_argument('on',type=Path);p.add_argument('off',type=Path);a=p.parse_args()
for d in (a.on,a.off):
 log=(d/'DATA/0/log').read_text()
 assert 'PICES_STEP step=4 ' in log, str(d)
 assert 'PICES_ERROR' not in (d/'stderr').read_text(),str(d)
on=(a.on/'DATA/0/PICES_P8c.tracer_log.0').read_text()
assert 'PICES_COVERAGE step=0 moved=0' in on
assert 'PICES_PATH step=4 ' in on
assert 'stored_valid=0' not in on
support=[s for s in on.splitlines() if s.startswith('PICES_SUPPORT ')]
assert len(support)==8
for s in support:
 assert float(re.search(r'weight_sum=([^ ]+)',s)[1])>0
 assert 'global_invalid_weights=0' in s
# Compare physical accepted-step ledger, excluding walltime and other logs.
for rank in range(24):
 x=a.on/f'DATA/{rank}';y=a.off/f'DATA/{rank}'
 def physics(d):
  return [s for s in (d/'log').read_text().splitlines() if s.startswith(('PICES_STEP ','PICES_EBA ','PICES_TA ','PICES_DIAGNOSTICS '))]
 assert physics(x)==physics(y),rank
files=list(a.on.glob('DATA/*/*.pices.state'))
assert files,'missing particle checkpoints'
for x in files:
 y=a.off/x.relative_to(a.on)
 assert x.read_bytes()==y.read_bytes(),str(x)
print('PASS: 24-rank accepted diagnostics identical; %d particle states byte-identical; positive MPI node support; paths valid'%len(files))
