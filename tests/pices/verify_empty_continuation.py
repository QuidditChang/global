"""Check prescribed-motion empty-element and normal MPI regression runs."""
from pathlib import Path
import argparse,re
p=argparse.ArgumentParser();p.add_argument('empty',type=Path);p.add_argument('normal',type=Path);a=p.parse_args()
worst=0
for case in (a.empty,a.normal):
 assert 'PICES_ERROR' not in (case/'stderr').read_text()
 for rank in range(24):
  log=(case/f'DATA/{rank}/log').read_text();assert 'PICES_STEP step=4 ' in log,(case,rank)
 log=(case/'DATA/0/log').read_text()
 for s in log.splitlines():
  if s.startswith('PICES_DIAGNOSTICS '):
   d=dict(re.findall(r'(\w+)=([^\s]+)',s));worst=max(worst,abs(float(d['closure_residual_nd'])))
   assert d['particle_count_global']=='393216'
   assert float(d['particle_Tmin_K'])>=0
  if s.startswith('PICES_PARTICLE_BALANCE '):assert 'delta=0 invalid_hosts=0 count_mismatch_elements=0' in s
 assert 'PICES_EMPTY_MONITOR step=4 empty='+('1' if case==a.empty else '0') in log
assert worst<1e-12,worst
assert 'event=begin' in (a.empty/'DATA/0/PICES_P8c.tracer_log.0').read_text()
print('PASS: empty and normal 24-rank runs reach step4; global particle conservation; max energy closure residual',worst)
