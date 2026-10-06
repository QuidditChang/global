"""P7 targeted preproduction audit. Engineering PASS is not production approval."""
import argparse,json,math,hashlib
from pathlib import Path
from verify_p6 import verify as base
from verify_p2 import require,fields,config

def verify(root,local=False,partial=False):
 root=Path(root);r=base(root,local,partial,'P7')
 target=json.loads((root/'input/production_target.json').read_text())
 require(hashlib.sha256((root/'input/production_target.cfg').read_bytes()).hexdigest()==target['target_sha256'],'production target snapshot hash')
 r['production_target']=target
 sharp=r['projection'].get('sharp_pices_base')
 require(sharp and sharp['max_qp']>0 and sharp['max_active']>0,'sharp test did not activate bounds')
 rows=[fields(l) for l in (root/'sharp_pices_base/DATA/0/log').read_text().splitlines() if l.startswith('PICES_PROJECTION ') and 'kind=absolute' in l]
 require(all(float(x['lower'])==0 and float(x['upper'])==1 for x in rows),'sharp bounds')
 for method in ['pg','pices']:
  d=root/f'assim_{method}_rheol7'
  if partial and not (d/'mpi_exit_code.txt').exists():continue
  c=config(d/'case.cfg');require(c['rheol']=='7' and c['TDEPV']=='on' and c['pices_checkpoint']=='off','rheology input')
 if not local:
  path=(root/'mpi_path.txt').read_text().strip();lib=(root/'binary_ldd.txt').read_text()
  require(path.startswith('/share/intel/') and '/mpi/intel64/bin/' in path,'non-Intel launcher')
  mpi_root=path.rsplit('/bin/',1)[0];require(mpi_root+'/lib/libmpi.so.12' in lib,'launcher/library mismatch')
  require('Intel' in (root/'mpi_provider_check.txt').read_text(),'MPI version provider')
 r['velocity_difference_percent']={}
 for variant in ['base','dt_half','dt_quarter','dt_eighth','rheol7']:
  pg=r['cases'].get('assim_pg_'+variant);pic=r['cases'].get('assim_pices_'+variant)
  if pg and pic:r['velocity_difference_percent'][variant]=100*(pic['final']['vrms']/pg['final']['vrms']-1)
 r['production_decision']='HOLD: rheol7 restart and full target configuration remain unsupported/unvalidated'
 return r
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('--local',action='store_true');p.add_argument('--partial',action='store_true');p.add_argument('--summary',type=Path);a=p.parse_args()
 try:r=verify(a.root,a.local,a.partial)
 except (AssertionError,ValueError,OSError,KeyError) as e:r=dict(status='FAIL',error=str(e))
 if a.summary:a.summary.write_text(json.dumps(r,indent=2)+'\n')
 print(json.dumps({k:v for k,v in r.items() if k in ['status','error','cases_completed','velocity_difference_percent','production_decision']},indent=2));raise SystemExit(r['status']=='FAIL')
