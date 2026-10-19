"""Build an isolated frozen-field replay; never relink the production executable.
--compile checks lib/bin against the checkpoint commit before compiling them.
The helper driver has separate provenance; checkpoint physics checks stay intact.
"""
import argparse, concurrent.futures, hashlib, json, os, subprocess
from pathlib import Path
from replay_checkpoint_instrument import instrument,VERSION
ROOT=Path(__file__).resolve().parents[2]
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('build',type=Path)
p.add_argument('--compile',action='store_true',help='compile isolated solver objects')
p.add_argument('--solver-commit',help='checkpoint commit; lib/bin must match exactly')
a=p.parse_args();b=a.build.resolve();b.mkdir(parents=True,exist_ok=True)
cc=os.environ.get('MPICC','mpicc')
commit=subprocess.check_output(['git','-C',str(ROOT),'rev-parse',a.solver_commit or 'HEAD'],text=True).strip()
flags=['-std=gnu99','-O2','-w','-Wno-error=implicit-function-declaration','-Wno-error=implicit-int','-Wno-error=int-conversion','-DUSE_GZDIR','-DPICES_SOLVER_COMMIT="'+commit+'"','-I'+str(ROOT/'lib')]
sources=(ROOT/'lib/Makefile.am').read_text().split('sources =',1)[1].split('EXTRA_DIST',1)[0]
sources=[ROOT/'lib'/v for v in sources.replace('\\','').split() if v.endswith('.c')]
objects=[str(b/(v.stem+'.o')) for v in sources]
diag_header=Path(__file__).with_name('replay_checkpoint_diagnostic.h').read_bytes()
(b/'replay_checkpoint_diagnostic.h').write_bytes(diag_header)
diag_source=b/'Pices_checkpoint.c'
diag_source.write_text(instrument((ROOT/'lib/Pices_checkpoint.c').read_text()))
if a.compile:
 subprocess.run(['git','-C',str(ROOT),'diff','--exit-code',commit,'--','lib','bin',':(exclude,glob)**/Makefile.in'],check=True)
 def compile_one(source):
  compile_source=diag_source if source.name=='Pices_checkpoint.c' else source
  result=subprocess.run([cc,*flags,'-c',str(compile_source),'-o',str(b/(source.stem+'.o'))],capture_output=True,text=True)
  if result.returncode:raise RuntimeError(str(source)+'\n'+result.stderr)
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
  list(pool.map(compile_one,sources+[ROOT/'bin/CitcomSFull.c']))
else:
 raise SystemExit('Use --compile to ensure all objects and diagnostic provenance agree')
helper=Path(__file__).with_name('motion_replay_driver.h').read_text()
s=(ROOT/'bin/Citcom.c').read_text().replace('int main(argc,argv)',helper+'\nint main(argc,argv)')
needle='      read_checkpoint(E);';assert s.count(needle)==1
s=s.replace(needle,needle+'\n      motion_replay(E);\n      MPI_Finalize();\n      return 0;')
f=b/'PICESMotionReplay.c';f.write_text(s)
subprocess.run([cc,*flags,'-c',str(f),'-o',str(b/'PICESMotionReplay.o')],check=True)
exe=b/'PICESMotionReplay'
subprocess.run([cc,*objects,str(b/'CitcomSFull.o'),str(b/'PICESMotionReplay.o'),'-lz','-lm','-o',str(exe)],check=True)
(b/'replay_build.json').write_text(json.dumps(dict(solver_commit=commit,isolated_compile=a.compile,
 driver_commit=subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip(),
 checkpoint_diagnostic_version=VERSION,checkpoint_diagnostic_sha256=hashlib.sha256(diag_source.read_bytes()+diag_header).hexdigest(),
 helper_sha256=hashlib.sha256(helper.encode()).hexdigest(),binary_sha256=hashlib.sha256(exe.read_bytes()).hexdigest()),indent=2)+'\n')
print(exe)
