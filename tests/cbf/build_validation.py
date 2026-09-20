"""Build standalone validation drivers without modifying generated build files.

This is a local test build, not the production Pyre installation recipe.
All sources, including the parameter parser, are compiled without patches.
"""
import argparse
import concurrent.futures
from pathlib import Path
import os
import shlex
import subprocess

ROOT=Path(__file__).resolve().parents[2]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build-dir',type=Path,required=True)
    args=parser.parse_args();build=args.build_dir.resolve();build.mkdir(parents=True,exist_ok=True)
    source_list=(ROOT/'lib/Makefile.am').read_text().split('sources =',1)[1].split('EXTRA_DIST',1)[0]
    sources=[ROOT/'lib'/v for v in source_list.replace('\\','').split() if v.endswith('.c')]
    (build/'VALIDATION_BUILD.txt').write_text('parser_workaround=False\nsource=%s\n'%ROOT)
    flags=[]
    libs=[]
    compiler=os.environ.get('MPICC','mpicc')
    command=[compiler,'-std=gnu99','-w','-Wno-error=implicit-function-declaration',
             '-Wno-error=implicit-int','-Wno-error=int-conversion','-O0','-g',
             '-DUSE_GZDIR','-I'+str(ROOT/'lib'),
             '-I'+str(ROOT/'tests/cbf'),*flags]
    sources += [ROOT/'bin/Citcom.c',ROOT/'bin/CitcomSFull.c']
    driver=(ROOT/'bin/Citcom.c').read_text().replace('int main(argc,argv)',
                        '#include "manufactured_state.h"\nint main(argc,argv)')
    driver=driver.replace('      initial_conditions(E);',
            '      initial_conditions(E);\n      CBF_manufactured_state(E);\n'
            '      save_CBF_if_due(E);\n      MPI_Finalize();\n      return 0;',1)
    manufactured=build/'Citcom_manufactured.c';manufactured.write_text(driver)
    sources.append(manufactured)
    def compile_one(source):
        obj=build/(source.stem+'.o')
        result=subprocess.run(command+['-c',str(source),'-o',str(obj)],capture_output=True,text=True)
        return source,obj,result
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        results=list(pool.map(compile_one,sources))
    for source,obj,result in results:
        if result.returncode:raise RuntimeError(str(source)+'\n'+result.stderr)
    for name,exclude in [('CitcomSFull','Citcom_manufactured.o'),('CBFManufactured','Citcom.o')]:
        objects=[str(obj) for source,obj,result in results if obj.name!=exclude]
        subprocess.run([compiler,*objects,*libs,'-lz','-lm','-o',str(build/name)],check=True)
        print(build/name)

if __name__=='__main__':main()
