"""Restart rejects changed composition physics before advancing any step."""
import argparse,json,os,re,shutil,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('build',type=Path);p.add_argument('suite',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True)
results={}
for name,key,value in [('wrong_buoyancy','buoyancy_ratio',','.join(['0']*24)),('wrong_conductivity','kC_ratio','0.9')]:
 d=a.output/name;d.mkdir();(d/'DATA').mkdir();s=(a.suite/'input/restart.cfg').read_text();s=re.sub(r'^'+key+r'=.*$',key+'='+value,s,flags=re.M);s=re.sub(r'^datadir_old=.*$','datadir_old='+str(a.suite.resolve()/'split/DATA/%RANK'),s,flags=re.M)
 (d/'case.cfg').write_text(s);shutil.copy(a.suite/'input/refstate_EBA_PICES_P4.txt',d);shutil.copytree(a.suite/'input/pices_p4_forcing',d/'pices_p4_forcing',ignore=shutil.ignore_patterns('._*'))
 with (d/'stdout').open('w') as out,(d/'stderr').open('w') as err:r=subprocess.run([os.environ.get('MPIEXEC','/usr/local/bin/mpiexec'),'--oversubscribe','-n','12',str(a.build.resolve()/'CitcomSFull'),'case.cfg'],cwd=d,stdout=out,stderr=err,timeout=90)
 assert r.returncode not in (0,8) and 'checkpoint metadata/physics/normalization/commit mismatch' in (d/'stderr').read_text(),name
 assert not any('PICES_STEP ' in f.read_text() for f in (d/'DATA').glob('*/log')),name
 results[name]='PASS'
(a.output/'summary.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(results))
