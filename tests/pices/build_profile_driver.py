"""Test-only full Stokes driver writing profiles at each normal output."""
from pathlib import Path
import argparse,os,subprocess
root=Path(__file__).resolve().parents[2]
p=argparse.ArgumentParser();p.add_argument('build',type=Path);a=p.parse_args();b=a.build.resolve()
s=(root/'bin/Citcom.c').read_text()
s='struct All_variables;\nvoid output_profiles(struct All_variables *,int);\n'+s
s=s.replace('output_time(E, E->monitor.solution_cycles);','output_time(E, E->monitor.solution_cycles);\n  output_profiles(E,E->monitor.solution_cycles);')
f=b/'PICESProfiles.c';f.write_text(s)
cc=os.environ.get('MPICC','/usr/local/bin/mpicc')
subprocess.run([cc,'-std=gnu99','-w','-Wno-error=implicit-function-declaration','-Wno-error=implicit-int','-DUSE_GZDIR','-I'+str(root/'lib'),'-c',str(f),'-o',str(b/'PICESProfiles.o')],check=True)
ss=(root/'lib/Makefile.am').read_text().split('sources =',1)[1].split('EXTRA_DIST',1)[0]
objs=[str(b/(Path(v).stem+'.o')) for v in ss.replace('\\','').split() if v.endswith('.c')]+[str(b/'CitcomSFull.o'),str(b/'PICESProfiles.o')]
subprocess.run([cc,*objs,'-lz','-lm','-o',str(b/'PICESProfiles')],check=True)
