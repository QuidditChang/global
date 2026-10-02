"""Link a test-only probe of real P/Q; requires a build_validation.py build."""
import argparse,os,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
def main():
 p=argparse.ArgumentParser();p.add_argument('build',type=Path);a=p.parse_args();b=a.build.resolve()
 source=(ROOT/'lib/Pices.c').read_text()+'\n'+(ROOT/'tests/pices/projection_probe.inc').read_text()
 (b/'PicesProjectionProbe.c').write_text(source)
 driver=(ROOT/'bin/Citcom.c').read_text();anchor='      initial_conditions(E);'
 assert driver.count(anchor)==1
 driver=driver.replace(anchor,anchor+'\n      pices_projection_probe(E);\n      MPI_Finalize();\n      return 0;')
 driver=driver.replace('int main(argc,argv)','void pices_projection_probe(struct All_variables *);\nint main(argc,argv)')
 (b/'ProjectionMain.c').write_text(driver)
 cc=os.environ.get('MPICC','/usr/local/bin/mpicc');flags=['-std=gnu99','-O2','-w','-Wno-error=implicit-function-declaration','-Wno-error=implicit-int','-Wno-error=int-conversion','-DUSE_GZDIR','-I'+str(ROOT/'lib')]
 for name in ['PicesProjectionProbe','ProjectionMain']:
  subprocess.run([cc,*flags,'-c',str(b/(name+'.c')),'-o',str(b/(name+'.o'))],check=True)
 sources=(ROOT/'lib/Makefile.am').read_text().split('sources =',1)[1].split('EXTRA_DIST',1)[0]
 objects=[b/(Path(v).stem+'.o') for v in sources.replace('\\','').split() if v.endswith('.c') and v!='Pices.c']
 subprocess.run([cc,*map(str,objects),str(b/'PicesProjectionProbe.o'),str(b/'ProjectionMain.o'),str(b/'CitcomSFull.o'),'-lz','-lm','-o',str(b/'PICESProjectionProbe')],check=True)
 print(b/'PICESProjectionProbe')
if __name__=='__main__':main()
