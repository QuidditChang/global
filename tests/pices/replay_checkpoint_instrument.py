"""Instrument a build-local checkpoint translation unit; keep lib/ untouched."""
import json,re
from pathlib import Path

VERSION=1

def instrument(source):
 marker='#define META_SIZE 4096'
 assert source.count(marker)==1
 source=source.replace(marker,marker+'\n#include "replay_checkpoint_diagnostic.h"')
 start=source.index('static void forcing_hash(');end=source.index('static void metadata(',start)
 part=source[start:end]
 part=part.replace('pices_sha_add(s,digest,64);','replay_diag_forcing(path,digest);replay_diag_add(s,digest,64,"forcing_digests");')
 def scalar_array(m):
  name,body=m.groups();expr=[s.strip() for s in body.split(',')]
  if name=='coupled':
   assert len(expr)==29
   assert [expr[i] for i in (18,20,21,22,23)]==['E->viscosity.SMOOTH','E->viscosity.TDEPV_AVE','E->viscosity.EQUIVDD','E->viscosity.equivddopt','E->viscosity.rheol_layers']
  idx='d' if name=='phase' else '-1'
  return m[0]+'\n'+''.join(' replay_diag_scalar(%s,%s,%s[%d]);\n'%(json.dumps(name+'.'+e),idx,name,i) for i,e in enumerate(expr))
 part=re.sub(r'double (params|ta|thermal|phase|composition|coupled)\[\]=\{(.*?)\};',scalar_array,part,flags=re.S)
 def add(m):
  data,size=m[1].split(',',1)
  return 'replay_diag_add(&s,%s,%s,%s);'%(data,size,json.dumps(data))
 part=re.sub(r'pices_sha_add\(&s,([^;]+)\);',add,part)
 source=source[:start]+part+source[end:]
 old='readtext(E,meta,actual);metadata(E,path,state,expected);'
 assert source.count(old)==1
 new='''readtext(E,meta,actual);
 replay_diag_begin(E);metadata(E,path,state,expected);
 { int different=strcmp(actual,expected)!=0, rejected=different, any_rejected;
   replay_diag_finish(E,actual,expected,different);
   /* Preserve the existing explicit compatibility gate; do not accept the
    * diagnostic zero-unused-fields candidate. Save all ranks before abort. */
   if(rejected) { char compatible[META_SIZE];strcpy(compatible,expected);
     if(coverage_restart_compatible(actual,compatible))rejected=0; }
   MPI_Allreduce(&rejected,&any_rejected,1,MPI_INT,MPI_MAX,E->parallel.world);
   if(E->parallel.me==0)fprintf(stderr,"PICES_REPLAY_CHECKPOINT_DIAG any_rejected=%d files=checkpoint_diag.*\\n",any_rejected);
   if(any_rejected)pices_fail(E,"checkpoint metadata/physics/normalization/commit mismatch (all-rank diagnostics saved)");
 }
'''
 return source.replace(old,new)
