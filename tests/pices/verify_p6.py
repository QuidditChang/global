"""P6 distributed projection diagnostics plus P5 physics and P4 restart audit."""
import argparse,json,math
from pathlib import Path
from verify_p5 import verify as matrix
from verify_p2 import verify as restart,require,config,fields

def verify(root,local=False,partial=False,stage="P6"):
 require(not partial or local,'partial audit allowed only for local validation')
 root=Path(root);report=matrix(root,local,partial,"P7" if stage=="P7" else "P5");checks={}
 cases=[r['name'] for r in json.loads((root/'input/matrix.json').read_text())['cases'] if r['method']=='pices']
 for name in cases+['restart_suite/'+s for s in ['continuous','split','restart']]:
  d=root/name
  if partial and not (d/'mpi_exit_code.txt').exists():continue
  require(config(d/'case.cfg').get('pices_projection')=='bounded_consistent','projection config '+name)
  rows=[fields(l) for l in (d/'DATA/0/log').read_text().splitlines() if l.startswith('PICES_PROJECTION ')]
  require(rows and {r['kind'] for r in rows}=={'absolute','signed'},'missing projection calls '+name)
  steps=[int(fields(l)['step']) for l in (d/'DATA/0/log').read_text().splitlines() if l.startswith('PICES_STEP ')]
  require([int(r['step']) for r in rows if r['kind']=='absolute']==steps,'absolute projection lifecycle '+name)
  require(all(sum(r['kind']=='signed' and int(r['step'])==step for r in rows)>=2 for step in steps),'signed projection lifecycle '+name)
  for r in rows:
   require(r['method']=='bounded_consistent_v1','projection implementation')
   require(all(math.isfinite(float(r[k])) for k in ['residual','tolerance','lower','upper']),'finite projection diagnostic')
   require(0<=float(r['residual'])<=float(r['tolerance']),'projection convergence')
  checks[name]=dict(calls=len(rows),max_cg=max(int(r['cg']) for r in rows),max_qp=max(int(r['qp']) for r in rows),max_active=max(int(r['active']) for r in rows))
 suite=root/'restart_suite';archived=root/'input/restart'
 for name,suffix in [('continuous',''),('split','_split'),('restart','_restart')]:
  d=suite/name;require((d/'case.cfg').read_bytes()==(archived/f'cmbhf_EBA_PICES_P4{suffix}.cfg').read_bytes(),'restart cfg provenance')
  require((d/'refstate_EBA_PICES_P4.txt').read_bytes()==(archived/'refstate_EBA_PICES_P4.txt').read_bytes(),'restart refstate')
  for f in (archived/'pices_p4_forcing').iterdir():
   if not f.name.startswith('._'):require((d/'pices_p4_forcing'/f.name).read_bytes()==f.read_bytes(),'restart forcing')
 report['restart']=restart(suite,True,'P4',64 if stage=='P7' else 4,32 if stage=='P7' else 2);report['projection']=checks
 if not local:
  for f in [stage.lower()+'_complete.txt','mpi_path.txt','binary_ldd.txt']:require((root/f).stat().st_size>0,'missing P6 provenance '+f)
 report['production_decision']='HOLD_PENDING_SCIENTIFIC_REVIEW_AND_PRODUCTION_PILOT'
 return report
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('--local',action='store_true');p.add_argument('--partial',action='store_true');p.add_argument('--summary',type=Path);a=p.parse_args()
 try:r=verify(a.root,a.local,a.partial)
 except (AssertionError,ValueError,OSError,KeyError) as e:r=dict(status='FAIL',error=str(e))
 if a.summary:a.summary.write_text(json.dumps(r,indent=2)+'\n')
 print(json.dumps({k:v for k,v in r.items() if k in ['status','error','cases_completed','production_decision']},indent=2));raise SystemExit(r['status'] not in ('PASS','PARTIAL_PASS'))
