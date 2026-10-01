"""P4 checkpoint rejects changed TA law/target/forcing and corrupted CBF cache."""
import argparse,json,os,re,shutil,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('build',type=Path);p.add_argument('run',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True)
source=a.run/'split';base=(a.run/'restart/case.cfg').read_text();results={}
for name,key,value in [('tau','lith_age_asml_tau_Ma','11'),('shape','lith_age_asml_exp','2'),('Tref',None,None),('age',None,None),('CBF_corrupt',None,None)]:
 d=a.output/name;d.mkdir();(d/'DATA').mkdir();cfg=base;old=source
 if name=='CBF_corrupt':
  old=d/'old';shutil.copytree(source/'DATA',old/'DATA')
  state=old/'DATA/0/PICES_P4.chkpt.0.2.pices.state';b=bytearray(state.read_bytes());b[-1]^=1;state.write_bytes(b)
 cfg=re.sub(r'^datadir_old=.*$','datadir_old='+str(old.resolve()/'DATA/%RANK'),cfg,flags=re.M)
 if key:cfg=re.sub('^'+key+'=.*$',key+'='+value,cfg,flags=re.M)
 (d/'case.cfg').write_text(cfg);ref=re.search(r'^refstate_file=(.*)$',cfg,re.M)[1];content=(source/ref).read_text()
 if name=='Tref':
  lines=[line.split() for line in content.splitlines()];lines[2][2]=str(float(lines[2][2])+.01);content='\n'.join(' '.join(x) for x in lines)+'\n'
 (d/ref).write_text(content);shutil.copytree(source/'pices_p4_forcing',d/'pices_p4_forcing')
 if name=='age':
  for f in (d/'pices_p4_forcing').glob('age.*'):f.write_text('61\n'*25)
 with (d/'stdout').open('w') as out,(d/'stderr').open('w') as err:
  r=subprocess.run([os.environ.get('MPIEXEC','/usr/local/bin/mpiexec'),'--oversubscribe','-n','12',str(a.build.resolve()/'CitcomSFull'),'case.cfg'],cwd=d,stdout=out,stderr=err,timeout=90)
 expected='checksum' if name=='CBF_corrupt' else 'metadata/physics'
 assert r.returncode==72 and expected in (d/'stderr').read_text(),(name,r.returncode,(d/'stderr').read_text()[-2000:])
 assert not any('PICES_STEP ' in f.read_text() for f in (d/'DATA').glob('*/log'))
 results[name]='PASS';print(name,'PASS',flush=True)
(a.output/'summary.json').write_text(json.dumps(results,indent=2)+'\n')
