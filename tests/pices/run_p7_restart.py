"""Run the P7 endurance configurations with local MPI, then audit."""
import argparse,os,shutil,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('build',type=Path);p.add_argument('runs',type=Path);p.add_argument('--stage',choices=['P4'],default='P4');p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True)
a.runs=a.runs/'pices_p7/restart'
for name,suffix,first,last in [('continuous','',1,64),('split','_split',1,32),('restart','_restart',33,64)]:
 d=a.output/name;d.mkdir();(d/'DATA').mkdir();shutil.copy(a.runs/('cmbhf_EBA_PICES_'+a.stage+suffix+'.cfg'),d/'case.cfg');shutil.copy(a.runs/('refstate_EBA_PICES_'+a.stage+'.txt'),d)
 if a.stage=='P4':shutil.copytree(a.runs/'pices_p4_forcing',d/'pices_p4_forcing')
 with (d/'solver.stdout').open('w') as out,(d/'solver.stderr').open('w') as err:
  r=subprocess.run([os.environ.get('MPIEXEC','/usr/local/bin/mpiexec'),'--oversubscribe','-n','12',str(a.build.resolve()/'CitcomSFull'),'case.cfg'],cwd=d,stdout=out,stderr=err,timeout=1800)
 (d/'mpi_exit_code.txt').write_text(str(r.returncode)+'\n');assert r.returncode in (0,8),(name,r.returncode,(d/'solver.stderr').read_text()[-2000:])

 assert all(f'PICES_STEP step={last} ' in (d/f'DATA/{rank}/log').read_text() for rank in range(12)),(name,'incomplete run despite exit code 8')
 print(name,'completed',flush=True)
from verify_p2 import verify
import json
r=verify(a.output,True,'P4',64,32);(a.output/'summary.json').write_text(json.dumps(r,indent=2)+'\n');print(r['status'])
