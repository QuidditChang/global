"""Link test-only prescribed-field driver against the real P1 build objects."""
from pathlib import Path
import argparse, os, subprocess
ROOT=Path(__file__).resolve().parents[2]
p=argparse.ArgumentParser();p.add_argument('build',type=Path);a=p.parse_args();build=a.build.resolve()
s=(ROOT/'bin/Citcom.c').read_text().replace('int main(argc,argv)', '#include "p1_manufactured.h"\nint main(argc,argv)')
s=s.replace('      initial_conditions(E);','      initial_conditions(E);\n      p1_test_initial(E);\n      p1_test_snapshot(E);')
s=s.replace('(E->next_buoyancy_field)(E);','(E->next_buoyancy_field)(E);\n    p1_test_snapshot(E);')
s=s.replace('general_stokes_solver(E);','p1_test_velocity(E);')
f=build/'PICESManufactured.c';f.write_text(s)
cc=os.environ.get('MPICC','/usr/local/bin/mpicc')
subprocess.run([cc,'-std=gnu99','-w','-Wno-error=implicit-function-declaration','-Wno-error=implicit-int', '-DUSE_GZDIR','-I'+str(ROOT/'lib'),'-I'+str(ROOT/'tests/pices'),'-c',str(f),'-o',str(build/'PICESManufactured.o')],check=True)
objs=[str(x) for x in build.glob('*.o') if x.name not in ('Citcom.o','Citcom_manufactured.o','PICESManufactured.o')]
subprocess.run([cc,*objs,str(build/'PICESManufactured.o'),'-lz','-lm','-o',str(build/'PICESManufactured')],check=True)
print(build/'PICESManufactured')
