"""Exercise real P3 Stokes failure paths and verify no checkpoint publication."""
import argparse,os,re,subprocess,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('build',type=Path);p.add_argument('runs',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True)
base=(a.runs/'cmbhf_EBA_PICES_P3.cfg').read_text();results={}
for name,key,token in [('outer','piterations','P3 Stokes outer solve did not converge'),('inner','vlowstep','P3 initial momentum correction did not converge')]:
 d=a.output/name;d.mkdir();(d/'DATA').mkdir();cfg=base
 value='52' if name=='outer' else '1'
 cfg=re.sub('^'+key+'=.*$',key+'='+value,cfg,flags=re.M) if re.search('^'+key+'=',cfg,re.M) else cfg+'\n'+key+'='+value+'\n'
 (d/'case.cfg').write_text(cfg);ref='refstate_EBA_PICES_P3.txt';(d/ref).write_bytes((a.runs/ref).read_bytes())
 with (d/'stdout').open('w') as out,(d/'stderr').open('w') as err:
  r=subprocess.run([os.environ.get('MPIEXEC','/usr/local/bin/mpiexec'),'--oversubscribe','-n','12',str(a.build.resolve()/'CitcomSFull'),'case.cfg'],cwd=d,stdout=out,stderr=err,timeout=90)
 assert r.returncode==72 and token in (d/'stderr').read_text(),(name,r.returncode,(d/'stderr').read_text()[-1500:])
 failed_step=1 if name=='outer' else 0
 assert re.search(r'PICES_ERROR rank=\d+ step='+str(failed_step)+':', (d/'stderr').read_text()),'unexpected failure step'
 manifests=list((d/'DATA').glob('*/*.pices.manifest'))
 assert all(int(f.name.split('.')[-3])<failed_step for f in manifests),'failed solve published a checkpoint'
 assert not any('PICES_STEP step=2 ' in f.read_text() for f in (d/'DATA').glob('*/log')),'advanced after failed solve'
 results[name]={'status':'PASS','exit':r.returncode,'failed_step':failed_step,'checkpoint_published_for_failed_step':False};print(name,'PASS',flush=True)
(a.output/'summary.json').write_text(json.dumps(results,indent=2)+'\n')
