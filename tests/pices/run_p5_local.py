"""Run selected or all exact P5 matrix inputs with the standalone binary."""
import argparse,os,shutil,subprocess,time,json,hashlib,platform
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('build',type=Path);p.add_argument('runs',type=Path);p.add_argument('--output',type=Path,required=True);p.add_argument('--cases',nargs='*');p.add_argument('--timeout',type=float,default=3600);a=p.parse_args();assert a.timeout>0,'timeout must be positive';a.output.mkdir(parents=True)
src=a.runs/'pices_p5';shutil.copytree(src,a.output/'input');rows=json.loads((src/'matrix.json').read_text())['cases']
assert not a.cases or set(a.cases)<={r['name'] for r in rows},'unknown case'
(a.output/'local_provenance.json').write_text(json.dumps(dict(binary=str(a.build.resolve()/'CitcomSFull'),binary_sha256=hashlib.sha256((a.build/'CitcomSFull').read_bytes()).hexdigest(),platform=platform.platform(),timeout_seconds=a.timeout,selected_cases=a.cases or [r['name'] for r in rows]),indent=2)+'\n')
for row in rows:
 if a.cases and row['name'] not in a.cases:continue
 d=a.output/row['name'];d.mkdir();(d/'DATA').mkdir();shutil.copy(src/'cases'/f"{row['name']}.cfg",d/'case.cfg');kind='assim' if row['scenario']=='assim' else 'transport';shutil.copy(src/f"refstate_{kind}_{row['nodes']}.txt",d/'refstate.txt');shutil.copytree(src/f"forcing_{row['nodes']}",d/'forcing')
 started=time.monotonic()
 with (d/'solver.stdout').open('w') as out,(d/'solver.stderr').open('w') as err:
  r=subprocess.run([os.environ.get('MPIEXEC','/usr/local/bin/mpiexec'),'--oversubscribe','-n','12',str(a.build.resolve()/'CitcomSFull'),'case.cfg'],cwd=d,stdout=out,stderr=err,timeout=a.timeout)
 (d/'mpi_exit_code.txt').write_text(str(r.returncode)+'\n');(d/'wall_seconds.txt').write_text(str(time.monotonic()-started)+'\n')
 assert r.returncode in (0,8),(row['name'],r.returncode,(d/'solver.stderr').read_text()[-1200:])
 assert f"P5_METRIC step={row['steps']} " in (d/'DATA/0/log').read_text(),row['name']
 print(row['name'],'completed',flush=True)
