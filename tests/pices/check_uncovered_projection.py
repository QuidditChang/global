"""Exercise actual consistent-projection code on covered/Dirichlet/free nodes."""
from pathlib import Path
import subprocess,tempfile
root=Path(__file__).resolve().parents[2]
s=(root/'lib/Pices_projection.c').read_text().replace('#include "global_defs.h"','').replace('#include "pices.h"','')
t=(root/'lib/Pices.c').read_text();valid=t[t.index('int pices_temperature_valid('):t.index('static double *array(')]
stub=r'''
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#define NCS 2
#define TBX 1
#define TBY 2
#define TBZ 4
#define MPI_SUM 1
#define MPI_MIN 2
#define MPI_MAX 3
#define MPI_DOUBLE 1
typedef int MPI_Op;
struct All_variables {
 struct {int nno,nel;} lmesh;
 struct {int me,world;} parallel;
 struct {int levmax;} mesh;
 struct {int solution_cycles;} monitor;
 struct {double Ttop,ref_temperature;} data;
 struct {int ntracers[2],ielement[2][3];} trace;
 int node[2][4];double T[2][4],sx[2][4][4];FILE *fp;
 void (*exchange_node_d)(struct All_variables *,double **,int);
};
void MPI_Allreduce(const void *a,void *b,int n,int t,int op,int comm){memcpy(b,a,n*sizeof(double));}
void share_stub(struct All_variables *E,double **x,int lev){}
void pices_fail(struct All_variables *E,const char *s){fprintf(stderr,"%s\n",s);exit(72);}
void tracer_temperature_weights(struct All_variables *E,int cap,int p,int *n,double *w){for(int a=1;a<=8;a++){n[a]=p;w[a]=a==1?1:0;}}
'''
main=r'''
int main(int argc,char **argv){
 struct All_variables E={0};double values[3]={0,.25,.5},out[4];
 E.fp=stdout;E.lmesh.nno=3;E.lmesh.nel=1;E.data.Ttop=300;E.data.ref_temperature=3000;
 E.trace.ntracers[1]=2;E.trace.ielement[1][1]=E.trace.ielement[1][2]=1;E.exchange_node_d=share_stub;
 E.T[1][1]=.25;E.T[1][2]=.5;E.T[1][3]=.7;E.node[1][3]=TBX;
 if(argc>1 && !strcmp(argv[1],"free"))E.node[1][3]=0;
 if(argc>1 && !strcmp(argv[1],"nan"))values[1]=NAN;
 assert(pices_temperature_valid(&E,-.05));assert(pices_temperature_valid(&E,1.5));
 assert(!pices_temperature_valid(&E,-.2));assert(!pices_temperature_valid(&E,INFINITY));
 pices_consistent_project(&E,values,out,1);
 assert(fabs(out[1]-.25)<1e-12 && fabs(out[2]-.5)<1e-12 && out[3]==.7);
 pices_consistent_project(&E,values,out,0);assert(out[3]==0);
 return 0;
}
'''
with tempfile.TemporaryDirectory() as tmp:
 p=Path(tmp);(p/'test.c').write_text(stub+valid+s+main)
 subprocess.run(['cc','-std=gnu99',str(p/'test.c'),'-lm','-o',str(p/'test')],check=True)
 for mode,expected in [('',0),('free',72),('nan',72)]:
  x=subprocess.run([str(p/'test'),mode],capture_output=True,text=True)
  assert x.returncode==expected,x.stderr
  if mode=='free':assert 'zero particle coverage at free node' in x.stderr
 print('PASS: Dirichlet zero-support absolute/signed projection; free-node failure; invalid temperature checks')
