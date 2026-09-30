"""Run the exact three P2 HPC configurations with local MPI, then audit."""
import argparse,os,shutil,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('build',type=Path);p.add_argument('runs',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True)
lsf=(a.runs/'cmbhf_EBA_PICES_P2.lsf').read_text();gate=lsf[lsf.index('awk -v first='):lsf.index('    cd "$RUN_ROOT"\ndone')]
for name,suffix,first,last in [('continuous','',1,4),('split','_split',1,2),('restart','_restart',3,4)]:
 d=a.output/name;d.mkdir();(d/'DATA').mkdir();shutil.copy(a.runs/('cmbhf_EBA_PICES_P2'+suffix+'.cfg'),d/'case.cfg');shutil.copy(a.runs/'refstate_EBA_PICES_P2.txt',d)
 with (d/'solver.stdout').open('w') as out,(d/'solver.stderr').open('w') as err:
  r=subprocess.run([os.environ.get('MPIEXEC','/usr/local/bin/mpiexec'),'--oversubscribe','-n','12',str(a.build.resolve()/'CitcomSFull'),'case.cfg'],cwd=d,stdout=out,stderr=err,timeout=180)
 (d/'mpi_exit_code.txt').write_text(str(r.returncode)+'\n');assert r.returncode in (0,8),(name,r.returncode,(d/'solver.stderr').read_text()[-2000:])
 subprocess.run(['bash','-c',gate],cwd=d,env=dict(os.environ,first_step=str(first),final_step=str(last)),check=True)
 print(name,'completed',flush=True)
subprocess.run(['python3',str(Path(__file__).with_name('verify_p2.py').resolve()),str(a.output),'--local','--summary',str(a.output/'summary.json')],check=True)
