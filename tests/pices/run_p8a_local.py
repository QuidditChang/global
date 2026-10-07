"""Execute the exact P8a inputs locally; audit with verify_p8a.py."""
import argparse,json,os,shutil,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('build',type=Path);p.add_argument('runs',type=Path);p.add_argument('--output',type=Path,required=True);p.add_argument('--cases',nargs='*');a=p.parse_args();a.output.mkdir(parents=True)
src=a.runs/'pices_p8a';shutil.copytree(src,a.output/'input',ignore=shutil.ignore_patterns('._*'))
for row in json.loads((src/'matrix.json').read_text())['cases']:
 name=row['name']
 if a.cases and name not in a.cases:continue
 d=a.output/name;d.mkdir();(d/'DATA').mkdir();shutil.copy(src/(name+'.cfg'),d/'case.cfg');shutil.copy(src/'refstate_EBA_PICES_P4.txt',d);shutil.copytree(src/'pices_p4_forcing',d/'pices_p4_forcing',ignore=shutil.ignore_patterns('._*'))
 with (d/'solver.stdout').open('w') as out,(d/'solver.stderr').open('w') as err:
  result=subprocess.run([os.environ.get('MPIEXEC','/usr/local/bin/mpiexec'),'--oversubscribe','-n','12',str(a.build.resolve()/'CitcomSFull'),'case.cfg'],cwd=d,stdout=out,stderr=err,timeout=900)
 (d/'mpi_exit_code.txt').write_text(str(result.returncode)+'\n')
 assert result.returncode in (0,8),(name,result.returncode,(d/'solver.stderr').read_text()[-1500:])
 assert all(f"PICES_STEP step={row['last']} " in (d/f'DATA/{rank}/log').read_text() for rank in range(12)),(name,'premature exit')
 print(name,'complete',flush=True)
