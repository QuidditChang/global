"""Regression for schema-4 trench history across an integer-age crossing."""
import argparse,json,os,re,shutil,struct,subprocess
from pathlib import Path
p=argparse.ArgumentParser()
p.add_argument('baseline',type=Path);p.add_argument('tb_fixed',type=Path)
p.add_argument('fixed',type=Path);p.add_argument('runs',type=Path)
p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-only',action='store_true')
a=p.parse_args();root=a.output.resolve()
if not a.verify_only:
    root.mkdir();src=a.runs/'pices_p8c'
    for case in ('continuous','restart'):
        d=root/case;d.mkdir();(d/'DATA').mkdir()
        s=(src/'continuous.cfg').read_text()
        values=dict(nprocz=2,mgunitz=1,lith_age_depth=.04,start_age=2.01,
                    accuracy='1e-4',maxstep=4,maxtotstep=5,checkpointFrequency=1,
                    solution_cycles_init=0)
        if case=='restart':values.update(restart='on',solution_cycles_init=3,datadir_old='../continuous/DATA/%RANK',datafile_old='PICES_P8c')
        for k,v in values.items():
            line=k+'='+str(v)
            s=re.sub('^'+k+'=.*$',lambda m:line,s,flags=re.M) if re.search('^'+k+'=',s,re.M) else s+'\n'+line+'\n'
        (d/'case.cfg').write_text(s);shutil.copy(src/'refstate.txt',d)
        shutil.copytree(src/'forcing',d/'forcing',ignore=shutil.ignore_patterns('._*'))
        for age in range(4):
            (d/'forcing'/('trench.%d.xyz'%age)).write_text('> L\n'+''.join('%s %s\n'%(lat,40+age*2) for lat in range(-70,71,10)))
    for case,build,prefix,expect in [('continuous',a.baseline,'',None),('restart',a.tb_fixed,'old.',72),('restart',a.fixed,'fixed.',None)]:
        env=os.environ.copy();env.pop('PICES_RESTART_COMPAT_COMMIT',None)
        if case=='restart':env['PICES_RESTART_COMPAT_COMMIT']='22399749866f452b3c06960e3d00eea5333a37f0'
        d=root/case
        with (d/(prefix+'stdout')).open('w') as out,(d/(prefix+'stderr')).open('w') as err:
            rc=subprocess.run([os.environ.get('MPIEXEC','/usr/local/bin/mpiexec'),'--oversubscribe','-n','24',str(build.resolve()/'CitcomSFull'),'case.cfg'],cwd=d,env=env,stdout=out,stderr=err,timeout=1800).returncode
        assert rc==expect if expect is not None else rc in (0,8),(case,prefix,rc)
        if case=='continuous':shutil.copy(d/'DATA/0/log',d/'DATA/0/PICES_P8c.log')
old=(root/'restart/old.stderr').read_text(errors='replace')
fixed=(root/'restart/fixed.stderr').read_text(errors='replace')
assert 'checkpoint metadata/physics/normalization/commit mismatch' in old
assert 'recovered_legacy_trench_geometry' in fixed
assert 'all_other_metadata=exact' in fixed and 'PICES_ERROR' not in fixed
raw_identical=0
for rank in range(24):
    for suffix in ('','.pices.state'):
        file=f'DATA/{rank}/PICES_P8c.chkpt.{rank}.4'+suffix
        before=bytearray((root/'continuous'/file).read_bytes())
        after=bytearray((root/'restart'/file).read_bytes())
        raw_identical+=before==after
        if not suffix:
            # Checkpoints.c explicitly writes unused P[0]/comp_el[*][0] slots.
            # Exclude only these slots, never any physical element/state.
            nx,ny,nz=struct.unpack_from('3i',before)
            nn=nx*ny*nz;ne=(nx-1)*(ny-1)*(nz-1)
            pressure_zero=44+16+2*(nn+1)*8+24
            before[pressure_zero:pressure_zero+8]=after[pressure_zero:pressure_zero+8]=bytes(8)
            offset=pressure_zero+(ne+1+3*nn)*8
            header=struct.unpack_from('9i',before,offset)
            assert header==struct.unpack_from('9i',after,offset)
            offset+=36+(header[8]+1)*((6+header[5])*8+4)
            nc=struct.unpack_from('i',before,offset+16)[0]
            assert nc==struct.unpack_from('i',after,offset+16)[0]
            offset+=20+16*nc
            assert offset+nc*(ne+1)*8==len(before)==len(after)
            for component in range(nc):
                unused=offset+component*(ne+1)*8
                before[unused:unused+8]=after[unused:unused+8]=bytes(8)
        assert before==after,(rank,suffix)
result=dict(status='PASS',tb_only_restart_rejected=True,original_fingerprint_exact_after_recovery=True,physical_state_files_identical=48,raw_files_identical=raw_identical,excluded_only_unused_pressure_and_comp_el_zero_slots=True,ranks=24,integer_age_crossing=True)
(root/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
