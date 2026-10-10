"""Check real ASCII output against gzip and coupled split/restart state."""
import argparse
import gzip
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
from verify_p2 import checkpoint

HERE = Path(__file__).resolve().parent

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('build',type=Path);p.add_argument('runs',type=Path)
    p.add_argument('--resume',action='store_true');p.add_argument('--output',type=Path,required=True);p.add_argument('--assembler',type=Path,required=True)
    a=p.parse_args();root=a.output.resolve();root.mkdir(parents=True,exist_ok=a.resume)
    for fmt in ('ascii-gz','ascii'):
        inputs=root/(fmt+'-inputs')
        if not inputs.exists():shutil.copytree(a.runs/'pices_p8c',inputs/'pices_p8c',ignore=shutil.ignore_patterns('._*'))
        for f in (inputs/'pices_p8c').glob('*.cfg'):
            f.write_text(f.read_text().replace('output_format=ascii-gz','output_format='+fmt))
        cmd=[sys.executable,str(HERE/'run_p8c_local.py'),str(a.build.resolve()),str(inputs),'--output',str(root/fmt)]
        cases=['continuous'] if fmt=='ascii-gz' else ['continuous','split','restart']
        pending=[]
        for case in cases:
            d=root/fmt/case
            if not d.exists():pending.append(case);continue
            assert a.resume and (d/'mpi_exit_code.txt').read_text().strip() in ('0','8')
            assert (d/'case.cfg').read_bytes()==(inputs/'pices_p8c'/(case+'.cfg')).read_bytes()
            last=2 if case=='split' else 4
            log_name='log' if fmt=='ascii-gz' else 'PICES_P8c.log'
            assert all(f'PICES_STEP step={last} ' in (d/f'DATA/{rank}'/log_name).read_text() for rank in range(12))
        if pending:subprocess.run(cmd+['--cases',*pending],check=True)
    raw=root/'ascii';packed=root/'ascii-gz';comparisons=0
    for rank in range(12):
        def folder(base,case):return base/case/'DATA'/str(rank)
        for step in range(5):
            for field in ('velo','visc','comp_nd'):
                ascii_file=folder(raw,'continuous')/f'PICES_P8c.{field}.{rank}.{step}'
                gz_file=folder(packed,'continuous')/str(step)/f'{field}.{rank}.{step}.gz'
                ascii_lines=ascii_file.read_text().splitlines();gz_lines=gzip.decompress(gz_file.read_bytes()).decode().splitlines()
                # ASCII composition has two metadata rows; gzip has one.
                ah=1 if field=='visc' else 2; gh=2 if field=='velo' else 1
                x=list(map(float,' '.join(ascii_lines[ah:]).split()));y=list(map(float,' '.join(gz_lines[gh:]).split()))
                assert len(x)==len(y) and all(math.isfinite(v) for v in x)
                assert all(math.isclose(u,v,rel_tol=1e-5,abs_tol=1e-7) for u,v in zip(x,y)),ascii_file
                comparisons+=1
            if step:
                for side in ('surf','botm'):
                    f=f'q.{side}.{rank}.{step}'
                    assert (folder(raw,'continuous')/f).read_bytes()==(folder(packed,'continuous')/f).read_bytes()
        for step in (0,2,4):
            f=f'PICES_P8c.chkpt.{rank}.{step}'
            assert checkpoint(folder(raw,'continuous')/f)[1]==checkpoint(folder(packed,'continuous')/f)[1]
        for case,steps in [('split',range(3)),('restart',range(2,5))]:
            for step in steps:
                for field in ('velo','visc','comp_nd'):
                    f=f'PICES_P8c.{field}.{rank}.{step}'
                    assert (folder(raw,'continuous')/f).read_bytes()==(folder(raw,case)/f).read_bytes()
                if step%2==0:
                    f=f'PICES_P8c.chkpt.{rank}.{step}'
                    assert checkpoint(folder(raw,'continuous')/f)[1]==checkpoint(folder(raw,case)/f)[1]
    # Restore the accepted state and exercise production k/CBF output without
    # another Stokes solve. Rebuild Output.o so this probe uses current sources.
    source_root=HERE.parents[1]
    probe=root/'probe-build';probe.mkdir(exist_ok=True)
    source=(source_root/'bin/Citcom.c').read_text()
    source=source.replace('  output_checkpoint(E);', '  output_checkpoint(E);\n  MPI_Finalize(); return 0;',1)
    (probe/'OutputProbe.c').write_text(source)
    cc=os.environ.get('MPICC','mpicc')
    flags=[cc,'-std=gnu99','-w','-Wno-error=implicit-function-declaration','-Wno-error=implicit-int','-Wno-error=int-conversion','-O2','-DUSE_GZDIR','-I'+str(source_root/'lib')]
    for src in (source_root/'lib/Output.c',probe/'OutputProbe.c'):
        subprocess.run(flags+['-c',str(src),'-o',str(probe/(src.stem+'.o'))],check=True)
    objects=[str(f) for f in a.build.resolve().glob('*.o') if f.name not in ('Citcom.o','Citcom_manufactured.o','Output.o')]
    subprocess.run([cc,*objects,str(probe/'Output.o'),str(probe/'OutputProbe.o'),'-lz','-lm','-o',str(probe/'OutputProbe')],check=True)
    for fmt in ('ascii','ascii-gz'):
        d=root/('probe-'+fmt);d.mkdir(exist_ok=True);(d/'DATA').mkdir(exist_ok=True)
        cfg_text=(raw/'restart/case.cfg').read_text()
        import re
        for key,val in dict(output_format=fmt,output_optional='comp_nd,k,surf,botm',solution_cycles_init='4',datadir_old=str(raw/'continuous/DATA/%RANK')).items():
            cfg_text=re.sub('^'+key+'=.*$',key+'='+val,cfg_text,flags=re.M)
        (d/'case.cfg').write_text(cfg_text);shutil.copy(raw/'continuous/refstate.txt',d)
        shutil.copytree(raw/'continuous/forcing',d/'forcing',dirs_exist_ok=True)
        with (d/'stdout').open('w') as out,(d/'stderr').open('w') as err:
            subprocess.run([os.environ.get('MPIEXEC','/usr/local/bin/mpiexec'),'--oversubscribe','-n','12',str(probe/'OutputProbe'),'case.cfg'],cwd=d,stdout=out,stderr=err,check=True,timeout=120)
    for rank in range(12):
        d=root/f'probe-ascii/DATA/{rank}'
        shutil.copy(raw/f'continuous/DATA/{rank}/PICES_P8c.coord.{rank}',d)
        x=list(map(float,(d/f'PICES_P8c.k.{rank}.4').read_text().split()))
        y=list(map(float,gzip.decompress((root/f'probe-ascii-gz/DATA/{rank}/4/k.{rank}.4.gz').read_bytes()).split()))
        assert x==y
        for side in ('surf','botm'):
            assert (d/f'q.{side}.{rank}.4').read_bytes()==(root/f'probe-ascii-gz/DATA/{rank}/q.{side}.{rank}.4').read_bytes()
            assert not (d/f'PICES_P8c.{side}.{rank}.4').exists()
    # Exercise the unchanged TA assembler against actual PICES ASCII files.
    cfg=root/'postproc.cfg';cfg.write_text(f'''[CitcomS.solver]
datadir={root}/probe-ascii/DATA/%RANK
datafile=PICES_P8c
[CitcomS.solver.mesher]
nodex=9
nodey=9
nodez=9
nproc_surf=12
nprocx=1
nprocy=1
nprocz=1
[CitcomS.solver.output]
output_format=ascii
output_optional=comp_nd
''')
    subprocess.run([sys.executable,str(a.assembler.resolve()),'--config',str(cfg),'--output-dir',str(root/'caps'),'--cap-fields','coord,velo,visc,k','--workers','1','4'],check=True,stdout=subprocess.DEVNULL)
    assert len(list((root/'caps').glob('PICES_P8c.cap*')))==12
    result=dict(status='PASS',ascii_gzip_fields_compared=comparisons,CBF_identical=True,accepted_binary_state_identical=True,ascii_restart_identical=True,TA_assembler_caps=12,conductivity_fields_identical=True,
                scope='Local standalone coupled fixture; HPC Pyre and full production duration not executed by this test.')
    (root/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':main()
