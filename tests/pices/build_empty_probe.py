"""Test-only prescribed-motion driver creating one empty element without losing particles."""
from pathlib import Path
import argparse,subprocess,os
root=Path(__file__).resolve().parents[2]
p=argparse.ArgumentParser();p.add_argument('build',type=Path);a=p.parse_args();b=a.build.resolve()
subprocess.run(['python3',str(root/'tests/pices/build_manufactured.py'),str(b)],check=True)
s=(b/'PICESManufactured.c').read_text()
inject=r'''
      /* Only the test driver moves these points. Production never resamples. */
      if(E->parallel.me==0) {
        int donor=0,p,k;double c[25];
        for(p=1;p<=E->trace.ntracers[1];p++)if(E->trace.ielement[1][p]==101){donor=p;break;}
        if(!donor)pices_fail(E,"empty test missing donor");
        for(k=0;k<E->composition.ncomp;k++)c[k]=E->composition.comp_el[1][k][100];
        for(p=1;p<=E->trace.ntracers[1];p++)if(E->trace.ielement[1][p]==100) {
          for(k=0;k<6;k++)E->trace.basicq[1][k][p]=E->trace.basicq[1][k][donor];
          E->trace.ielement[1][p]=101;
        }
        /* Collective calls below must also be entered by other ranks. */
        recount_tracers_of_flavors(E);
        fill_composition(E);
        for(k=0;k<E->composition.ncomp;k++)if(c[k]!=E->composition.comp_el[1][k][100])pices_fail(E,"empty test changed retained composition");
      } else {
        recount_tracers_of_flavors(E);fill_composition(E);
      }
'''
s=s.replace('      p1_test_initial(E);','      p1_test_initial(E);'+inject)
f=b/'PICESEmptyProbe.c';f.write_text(s)
cc=os.environ.get('MPICC','/usr/local/bin/mpicc')
subprocess.run([cc,'-std=gnu99','-w','-Wno-error=implicit-function-declaration','-Wno-error=implicit-int','-DUSE_GZDIR','-I'+str(root/'lib'),'-I'+str(root/'tests/pices'),'-c',str(f),'-o',str(b/'PICESEmptyProbe.o')],check=True)
ss=(root/'lib/Makefile.am').read_text().split('sources =',1)[1].split('EXTRA_DIST',1)[0]
objs=[str(b/(Path(v).stem+'.o')) for v in ss.replace('\\','').split() if v.endswith('.c')]+[str(b/'CitcomSFull.o'),str(b/'PICESEmptyProbe.o')]
subprocess.run([cc,*objs,'-lz','-lm','-o',str(b/'PICESEmptyProbe')],check=True)
