"""Reproduce schema-4 shallow-TA restart failure, then verify exact recovery.

The baseline build must be the original 2239974 executable (or its logging-only
variant). The fixed build must include initial-interface TB recovery.
"""
import argparse,json,os,re,shutil,subprocess
from pathlib import Path
p=argparse.ArgumentParser()
p.add_argument('baseline',type=Path);p.add_argument('fixed',type=Path);p.add_argument('runs',type=Path)
p.add_argument('--output',type=Path,required=True);p.add_argument('--verify-only',action='store_true')
a=p.parse_args();root=a.output.resolve()
if not a.verify_only:
    root.mkdir()
    src=a.runs/'pices_p8c'
    for case in ('continuous','restart'):
        d=root/case;d.mkdir();(d/'DATA').mkdir()
        s=(src/'continuous.cfg').read_text()
        values=dict(nprocz=2,mgunitz=1,lith_age_depth=.04,maxstep=2,maxtotstep=3,checkpointFrequency=1)
        if case=='restart':values.update(restart='on',solution_cycles_init=1,datadir_old='../continuous/DATA/%RANK',datafile_old='PICES_P8c')
        for k,v in values.items():
            line=k+'='+str(v)
            s=re.sub('^'+k+'=.*$',lambda m:line,s,flags=re.M) if re.search('^'+k+'=',s,re.M) else s+'\n'+line+'\n'
        (d/'case.cfg').write_text(s);shutil.copy(src/'refstate.txt',d)
        shutil.copytree(src/'forcing',d/'forcing',ignore=shutil.ignore_patterns('._*'))
    for case,build,prefix,expect in [('continuous',a.baseline,'',None),('restart',a.baseline,'old.',72),('restart',a.fixed,'fixed.',None)]:
        env=os.environ.copy();env.pop('PICES_RESTART_COMPAT_COMMIT',None)
        if prefix=='fixed.':env['PICES_RESTART_COMPAT_COMMIT']='22399749866f452b3c06960e3d00eea5333a37f0'
        d=root/case
        with (d/(prefix+'stdout')).open('w') as out,(d/(prefix+'stderr')).open('w') as err:
            rc=subprocess.run([os.environ.get('MPIEXEC','/usr/local/bin/mpiexec'),'--oversubscribe','-n','24',str(build.resolve()/'CitcomSFull'),'case.cfg'],cwd=d,env=env,stdout=out,stderr=err,timeout=1800).returncode
        assert rc==expect if expect is not None else rc in (0,8),(case,prefix,rc)
old=(root/'restart/old.stderr').read_text(errors='replace')
fixed=(root/'restart/fixed.stderr').read_text(errors='replace')
assert 'checkpoint metadata/physics/normalization/commit mismatch' in old
assert 'recovered_inactive_interface_TB=1944 from_verified_step=0' in fixed
assert 'all_other_metadata=exact' in fixed and 'PICES_ERROR' not in fixed
for rank in range(24):
    for case in ('continuous','restart'):
        assert 'PICES_STEP step=2 ' in (root/case/f'DATA/{rank}/log').read_text()
    for suffix in ('','.pices.state'):
        file=f'DATA/{rank}/PICES_P8c.chkpt.{rank}.2'+suffix
        assert (root/'continuous'/file).read_bytes()==(root/'restart'/file).read_bytes(),(rank,suffix)
result=dict(status='PASS',old_restart_rejected=True,original_fingerprint_exact_after_recovery=True,accepted_payloads_identical=48,ranks=24,nprocz=2,lith_age_depth=.04)
(root/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
