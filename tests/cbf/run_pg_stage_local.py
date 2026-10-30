"""Three-step, 12-rank coupled PG/TA/CBF regression using synthetic forcing.
Build with build_validation.py first. MPIEXEC must match the build's mpicc.
"""
import argparse, os, re, subprocess
from pathlib import Path
from verify_native_outputs import verify
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('build',type=Path);p.add_argument('output',type=Path)
a=p.parse_args();root=a.output.resolve();root.mkdir(parents=True,exist_ok=False)
d={}
for line in Path(__file__).with_name('smoke.cfg').read_text().splitlines():
    if '=' in line and not line.startswith('#'):
        k,v=line.split('=',1);d[k]=v
d.update(maxstep=3,maxtotstep=4,cpu_limits_in_seconds=10000000,
    datadir='DATA/%RANK',lith_age=1,lith_age_asml=1,lith_age_time=1,
    pg_ta_before_heat='on',lith_age_depth=.25,lith_age_asml_tau_Ma=.1,
    lith_age_asml_exp=3,max_plate_age_Ma=70,temperature_bound_adj=0,start_age=2,
    lith_age_file='forcing/age.',flag_depth_file='forcing/trench.',
    flag_depth_new_file='forcing/trench.',tf_file='forcing/trench.',
    fixed_timestep='1e-7',temperature_audit='on',kd_upper_prefactor=4,
    kd_lower_prefactor=4,p5_case='assim',p5_prescribed_velocity='off',
    pices_checkpoint='off')
(root/'case.cfg').write_text(''.join(f'{k}={v}\n' for k,v in d.items()))
(root/'refstate.txt').write_text('1 1 .6 1 1\n'*5)
(root/'forcing').mkdir()
for rank in range(12):(root/'DATA'/str(rank)).mkdir(parents=True)
for age in range(4):
    (root/'forcing'/f'trench.{age}.xyz').write_text('')
    for cap in range(12):(root/'forcing'/f'age.{age}.{cap}').write_text('30\n'*25)
with (root/'stdout').open('w') as out,(root/'stderr').open('w') as err:
    r=subprocess.run([os.environ.get('MPIEXEC','mpiexec'),'--oversubscribe','-n','12',
        str(a.build.resolve()/'CitcomSFull'),'case.cfg'],cwd=root,stdout=out,stderr=err,timeout=120)
# This legacy standalone driver's normal finalizer calls exit(8).
assert r.returncode==8,(r.returncode,(root/'stderr').read_text()[-1500:])
for rank in range(12):
    log=(root/'DATA'/str(rank)/'log').read_text()
    assert 'stage=post_assimilation' not in log
    for step in range(1,4):
        ta=f'stage=pg_ta_before_heat step={step} '
        predictor=f'stage=predictor step={step} '
        assert log.count(ta)==1 and log.index(ta)<log.index(predictor)
        ta_row=next(x for x in log.splitlines() if ta in x)
        pred_row=next(x for x in log.splitlines() if predictor in x)
        def val(row,key):return float(re.search(key+r'=([^ ]+)',row)[1])
        assert abs(val(ta_row,'time_nd')-val(pred_row,'time_nd')-val(ta_row,'dt_nd'))<1e-13
    assert 'nonfinite=0' in log
log=(root/'DATA/0/log').read_text()
assert log.count('PG_CBF_ENDPOINT ')==3
assert 'ta_before_heat=1' in log
assert not list(root.rglob('q.surf.*.0'))
for step in range(1,4):verify(root,step,12)
print('PASS: TA once before predictor, arrival clock restored, three accepted CBF outputs, initial CBF absent')
