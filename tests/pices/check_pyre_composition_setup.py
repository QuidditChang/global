"""Reproduce Pyre's unset derived capability at the common C setup boundary.

Uses real C objects and the existing coupled MPI gate. This tests the entry
boundary, not the unavailable local Python 2/Pythia runtime.
"""
import argparse,json,os,shutil,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
p=argparse.ArgumentParser();p.add_argument('build',type=Path);p.add_argument('runs',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output=a.output.resolve();a.output.mkdir(parents=True)
source=(ROOT/'lib/Instructions.c').read_text()
entry='void initial_mesh_solver_setup(struct All_variables *E)\n{'
assert source.count(entry)==1
source=source.replace(entry,entry+'\n    E->composition.on=0;\n    E->composition.icompositional_rheology=0;',1)
objects=[str(f.resolve()) for f in a.build.glob('*.o') if f.name not in ('Instructions.o','Citcom_manufactured.o')]
assert len(objects)>20
for variant in ('legacy','fixed'):
    b=a.output/variant;b.mkdir()
    text=source.replace('    composition_set_capabilities(E);','',1) if variant=='legacy' else source
    c=b/'Instructions.c';c.write_text(text);obj=b/'Instructions.o'
    compiler=os.environ.get('MPICC','mpicc')
    subprocess.run([compiler,'-std=gnu99','-w','-Wno-error=implicit-function-declaration','-Wno-error=implicit-int','-Wno-error=int-conversion','-O2','-DUSE_GZDIR','-I'+str(ROOT/'lib'),'-c',str(c),'-o',str(obj)],check=True)
    subprocess.run([compiler,*objects,str(obj),'-lz','-lm','-o',str(b/'CitcomSFull')],check=True)
src=a.runs/'pices_p8c';d=a.output/'legacy_case';d.mkdir();(d/'DATA').mkdir()
shutil.copy(src/'continuous.cfg',d/'case.cfg');shutil.copy(src/'refstate.txt',d);shutil.copytree(src/'forcing',d/'forcing',ignore=shutil.ignore_patterns('._*'))
with (d/'stdout').open('w') as out,(d/'stderr').open('w') as err:
    r=subprocess.run([os.environ.get('MPIEXEC','/usr/local/bin/mpiexec'),'--oversubscribe','-n','12',str(a.output/'legacy/CitcomSFull'),'case.cfg'],cwd=d,stdout=out,stderr=err,timeout=180)
assert r.returncode==72 and 'P8a reclassification requires 25 flavors, primordial 24 and P4 age forcing' in (d/'stderr').read_text()
print('PASS: reproduced legacy Pyre capability abort',flush=True)
subprocess.run([sys.executable,str(ROOT/'tests/pices/run_p8c_local.py'),str(a.output/'fixed'),str(a.runs.resolve()),'--output',str(a.output/'coupled')],check=True)
subprocess.run([sys.executable,str(ROOT/'tests/pices/verify_p8c.py'),str(a.output/'coupled'),'--local','--summary',str(a.output/'summary.json')],check=True)
