"""Summarize all-rank replay checkpoint diagnostics, including failed launches."""
import argparse,collections,json
from pathlib import Path

def audit(root,ranks):
 differences=collections.Counter();examples={};missing=[];zero_matches=[];unused={};bad=[]
 for rank in range(ranks):
  stem=root/'DATA'/str(rank)/('global.checkpoint_diag.%d'%rank)
  paths=[Path(str(stem)+s) for s in ('.actual.json','.expected.json','.tsv')]
  if not all(p.is_file() for p in paths):missing.append(rank);continue
  actual,expected=[json.loads(p.read_text()) for p in paths[:2]]
  keys=[k for k in actual.keys()|expected.keys() if actual.get(k)!=expected.get(k)]
  if keys:bad.append(rank)
  for k in keys:
   differences[k]+=1
   if k not in examples:examples[k]=dict(rank=rank,saved=actual.get(k),running=expected.get(k))
  for line in paths[2].read_text().splitlines():
   fields=line.split('\t')
   if fields[0]=='candidate_zero_unused_fields_sha256' and fields[1]==actual['physics_mesh_sha256']:
    zero_matches.append(rank)
   if fields[0]=='scalar' and fields[1].split('.')[-1] in ('SMOOTH','TDEPV_AVE','EQUIVDD','equivddopt','rheol_layers'):
    unused.setdefault(fields[1],collections.Counter())[fields[3]]+=1
 return dict(status='INCOMPLETE' if missing else 'MISMATCH' if bad else 'EXACT',ranks_expected=ranks,ranks_captured=ranks-len(missing),missing_ranks=missing,mismatch_ranks=bad,field_counts=dict(differences),examples=examples,zero_unused_candidate_matches_saved_physics=zero_matches,running_unused_values={k:dict(v) for k,v in unused.items()},source_component_values_available=False,scope='Saved and running metadata compared exactly. Named scalars and block hashes describe running state only; legacy checkpoints lack source component snapshots. Zero-unused candidate is diagnostic only, never an acceptance rule.')

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('case',type=Path);p.add_argument('--ranks',type=int,default=384);a=p.parse_args()
 result=audit(a.case,a.ranks);(a.case/'checkpoint_diagnostic_summary.json').write_text(json.dumps(result,indent=2)+'\n')
 print(result['status'],'captured',result['ranks_captured'],'of',a.ranks,'field_counts',result['field_counts'])
 if result['missing_ranks']:raise SystemExit(2)
if __name__=='__main__':main()
