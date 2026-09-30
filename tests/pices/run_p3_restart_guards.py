"""Reject changed P3 physics when restarting otherwise intact checkpoints."""
import argparse,json,os,re,shutil,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('build',type=Path);p.add_argument('run',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True)
source=a.run/'split';base=(a.run/'restart/case.cfg').read_text();results={}
for name,key,value in [('conductivity','kT_exponent','0.4'),('phase','phase_delta_s','-0.03,-0.02,0.02'),('source','Q0','2'),('capacity',None,None)]:
 d=a.output/name;d.mkdir();(d/'DATA').mkdir();cfg=base
 cfg=re.sub(r'^datadir_old=.*$', 'datadir_old='+str(source.resolve()/'DATA/%RANK'),cfg,flags=re.M)
 if key:cfg=re.sub('^'+key+'=.*$',key+'='+value,cfg,flags=re.M)
 (d/'case.cfg').write_text(cfg);ref=re.search(r'^refstate_file=(.*)$',cfg,re.M)[1]
 content=(source/ref).read_text()
 if name=='capacity':
  lines=[line.split() for line in content.splitlines()];lines[2][4]=str(float(lines[2][4])+.1);content='\n'.join(' '.join(x) for x in lines)+'\n'
 (d/ref).write_text(content)
 with (d/'stdout').open('w') as out,(d/'stderr').open('w') as err:
  r=subprocess.run([os.environ.get('MPIEXEC','/usr/local/bin/mpiexec'),'--oversubscribe','-n','12',str(a.build.resolve()/'CitcomSFull'),'case.cfg'],cwd=d,stdout=out,stderr=err,timeout=90)
 assert r.returncode==72 and 'metadata/physics' in (d/'stderr').read_text(),(name,r.returncode,(d/'stderr').read_text()[-2000:])
 assert not any('PICES_STEP ' in f.read_text() for f in (d/'DATA').glob('*/log'))
 results[name]='PASS';print(name,'PASS',flush=True)
(a.output/'summary.json').write_text(json.dumps(results,indent=2)+'\n')
