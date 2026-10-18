"""Check replay numeric records through the actual unmodified C parser.
Usage: python3 check_replay_parser.py RUNTIME_CFG ORIGINAL_CFG OUTPUT_DIR
No MPI launch, solver changes, source checkpoint writes or production rebuild.
"""
import argparse,os,subprocess
from pathlib import Path
from motion_window import parameters,restore_input_precision,standalone_runtime,standalone_mesh,replace

HARNESS=r'''
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "global_defs.h"
#include "parsing.h"
/* Only the parser's incidental cfg-copy operation is suppressed. */
void unique_copy_file(struct All_variables *E,char *name,char *tag) {}
int main(int argc,char **argv) {
 struct All_variables *E=calloc(1,sizeof(*E));
 FILE *f;char key[64];int n,i,bad=0,records=0;double expected[128],d[128];float v[128];
 if(argc!=3)return 2;
 setup_parser(E,argv[1]);f=fopen(argv[2],"r");if(!f)return 2;
 while(fscanf(f,"%63s %d",key,&n)==2) {
  if(n<1||n>128)return 2;
  for(i=0;i<n;i++){if(fscanf(f,"%lf",&expected[i])!=1)return 2;d[i]=NAN;v[i]=NAN;}
  input_double_vector(key,n,d,0);input_float_vector(key,n,v,0);
  for(i=0;i<n;i++)if(d[i]!=expected[i]||v[i]!=(float)expected[i]) {
   fprintf(stderr,"MISMATCH %s[%d] expected=%.17g double=%.17g float=%.9g\n",key,i,expected[i],d[i],v[i]);bad++;
  }
  records++;
 }
 fclose(f);shutdown_parser(E);free(E);
 printf("records=%d mismatches=%d\n",records,bad);return bad?1:0;
}
'''

def main():
 p=argparse.ArgumentParser(description=__doc__)
 for key in ('runtime','original','output'):p.add_argument(key,type=Path)
 a=p.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
 raw=a.runtime.read_text();original=a.original.read_text()
 # Reproduce r3's precision restoration verbatim, before the lexical fix.
 broken=raw
 for key,value in parameters(original).items():
  shown=parameters(raw).get(key)
  if shown is None:continue
  try:changed=[float(x) for x in value.split(',')]!=[float(x) for x in shown.split(',')]
  except ValueError:continue
  if changed:broken=replace(broken,key,value)
 broken,_=standalone_runtime(broken);broken,_,_=standalone_mesh(broken)
 fixed,restored=restore_input_precision(raw,original)
 fixed,_=standalone_runtime(fixed);fixed,_,_=standalone_mesh(fixed)
 records=[]
 for key,value in parameters(fixed).items():
  try:values=[float(v) for v in value.split(',')]
  except ValueError:continue
  records.append(key+' '+str(len(values))+' '+' '.join(format(v,'.17g') for v in values))
 (out/'expected.txt').write_text('\n'.join(records)+'\n');(out/'parser_check.c').write_text(HARNESS)
 root=Path(__file__).resolve().parents[2]
 cc=os.environ.get('MPICC','mpicc')
 subprocess.run([cc,'-std=gnu99','-w','-Wno-error=implicit-function-declaration','-Wno-error=implicit-int','-I'+str(root/'lib'),str(root/'lib/Parsing.c'),str(out/'parser_check.c'),'-lm','-o',str(out/'parser_check')],check=True)
 for name,text in [('r3',broken),('fixed',fixed)]:
  path=out/(name+'.cfg');path.write_text(text)
  r=subprocess.run([str(out/'parser_check'),str(path),str(out/'expected.txt')],capture_output=True,text=True)
  (out/(name+'.log')).write_text(r.stdout+r.stderr)
  if name=='r3':assert r.returncode==1 and 'MISMATCH phase_delta_s[1]' in r.stderr,r
  else:assert r.returncode==0,r.stderr
  print(name,r.stdout.strip())
 print('PASS: r3 array truncation reproduced; all fixed numeric records match, float and double; restored',len(restored))

if __name__=='__main__':main()
