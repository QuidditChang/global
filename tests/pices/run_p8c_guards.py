"""Reject changed checkpoint physics, input contents and corrupt accepted state."""
import argparse,json,os,shutil,subprocess
from pathlib import Path
from prepare_p8a import edit
p=argparse.ArgumentParser();p.add_argument('build',type=Path);p.add_argument('baseline',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True)
shutil.copytree(a.baseline/'split',a.output/'split',ignore=shutil.ignore_patterns('._*'))
results={}
for name,changes in [('qvis',dict(qvis_cohesion_pa=2e7)),('cold_scale',dict(cold_scale=.5)),('step_limit',dict(pices_max_timestep_Ma=.04)),('plate_file',{}),('state_corruption',{})]:
    d=a.output/name;d.mkdir();(d/'DATA').mkdir()
    (d/'case.cfg').write_text(edit((a.baseline/'restart/case.cfg').read_text(),changes))
    shutil.copy(a.baseline/'restart/refstate.txt',d);shutil.copytree(a.baseline/'restart/forcing',d/'forcing',ignore=shutil.ignore_patterns('._*'))
    if name=='plate_file':
        f=d/'forcing/bvel.2.0';f.write_text(f.read_text().replace('0.03000000','0.03100000'))
    if name=='state_corruption':
        f=a.output/'split/DATA/0/PICES_P8c.chkpt.0.2.pices.state';b=bytearray(f.read_bytes());b[-1]^=1;f.write_bytes(b)
    with (d/'stdout').open('w') as out,(d/'stderr').open('w') as err:
        r=subprocess.run([os.environ.get('MPIEXEC','/usr/local/bin/mpiexec'),'--oversubscribe','-n','12',str(a.build.resolve()/'CitcomSFull'),'case.cfg'],cwd=d,stdout=out,stderr=err,timeout=180)
    text=(d/'stderr').read_text()
    assert r.returncode not in (0,8) and ('checkpoint metadata/physics/normalization/commit mismatch' in text or 'accepted velocity checksum mismatch' in text),(name,r.returncode,text[-1000:])
    assert not any('PICES_STEP ' in f.read_text() for f in (d/'DATA').glob('*/log')),name
    results[name]='PASS'
(a.output/'summary.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(results))
