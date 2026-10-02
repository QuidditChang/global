"""Copy the reviewed P5 matrix and P4 restart inputs, opting PIC into P6."""
import argparse,json,shutil
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('runs',type=Path);a=p.parse_args();r=a.runs;d=r/'pices_p6'
shutil.copytree(r/'pices_p5',d,dirs_exist_ok=True,ignore=shutil.ignore_patterns('._*'))
for f in (d/'cases').glob('*pices*.cfg'):
 if f.name.startswith('._'):continue
 s=f.read_text()
 if 'pices_projection=' in s:raise ValueError('P5 input already selects a projection')
 f.write_text(s+'\n# P6 bounded consistent absolute projection; linear signed increments.\npices_projection=bounded_consistent\n')
m=json.loads((d/'matrix.json').read_text());m['projection']='bounded_consistent_v1';(d/'matrix.json').write_text(json.dumps(m,indent=2)+'\n')
q=d/'restart';q.mkdir(exist_ok=True)
for suffix in ['', '_split','_restart']:
 f=r/f'cmbhf_EBA_PICES_P4{suffix}.cfg';(q/f.name).write_text(f.read_text()+'\npices_projection=bounded_consistent\n')
shutil.copy(r/'refstate_EBA_PICES_P4.txt',q)
shutil.copytree(r/'pices_p4_forcing',q/'pices_p4_forcing',dirs_exist_ok=True,ignore=shutil.ignore_patterns('._*'))
print('P6: 46 matrix cases + continuous/split/restart')
