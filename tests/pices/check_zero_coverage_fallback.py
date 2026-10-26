"""Compile the production lumped projection; test uncovered-node semantics."""
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
s = (ROOT / 'lib/Pices.c').read_text()
valid = s[s.index('int pices_temperature_valid('):s.index('static double *array(')]
project = s[s.index('static void project('):s.index('static double energy(', s.index('static void project('))]
stub = r'''
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
struct All_variables {
 struct {int nno,nel;} lmesh;
 struct {int me;} parallel;
 struct {int solution_cycles;} monitor;
 struct {double Ttop,ref_temperature;} data;
 struct {int consistent_projection;} pices;
 struct {int ntracers[2],ielement[2][3];} trace;
 double T[2][4],sx[2][4][4];int fixed[4];FILE *fp;
};
static int remote,exchange_call,negative;
static void pices_fail(struct All_variables *E,const char *s){fprintf(stderr,"%s\n",s);exit(72);}
static double *array(struct All_variables *E,int n){return calloc(n,sizeof(double));}
static int fixed(struct All_variables *E,int n){return E->fixed[n];}
static void exchange(struct All_variables *E,double *v){
 /* Production calls exchange(out), exchange(den). Simulate a shared node
  * with no LOCAL particle contribution but valid REMOTE coverage. */
 if(remote)v[3]+=(exchange_call++%2==0)?.8:2.;
 if(negative)v[3]=-1.;
}
static void tracer_temperature_weights(struct All_variables *E,int cap,int p,int *n,double *w){
 for(int a=1;a<=8;a++){n[a]=(a==2?2:1);w[a]=a==1?.25:(a==2?.75:0);}
}
static void pices_consistent_project(struct All_variables *E,const double *v,double *out,int b){abort();}
'''
main = r'''
int main(int argc,char **argv){
 struct All_variables E={0};double values[3]={0,.6,.2},out[4],stage[4]={0,.1,.2,.9};
 const char *mode=argc>1?argv[1]:"";
 E.fp=stdout;E.lmesh.nno=3;E.lmesh.nel=1;E.data.Ttop=300;E.data.ref_temperature=3400;
 E.trace.ntracers[1]=1;E.trace.ielement[1][1]=1;
 E.T[1][1]=.1;E.T[1][2]=.2;E.T[1][3]=.7;
 if(!strcmp(mode,"bad_fallback"))E.T[1][3]=-.2;
 if(!strcmp(mode,"nan"))values[1]=NAN;
 if(!strcmp(mode,"negative"))negative=1;
 if(!strcmp(mode,"remote"))remote=1;
 if(!strcmp(mode,"fixed"))E.fixed[3]=1;
 if(!strcmp(mode,"empty"))E.trace.ntracers[1]=0;
 project(&E,values,out,1,E.T[1]);
 assert(fabs(out[3]-(remote?.4:.7))<1e-14);
 if(E.trace.ntracers[1])assert(fabs(out[1]-.6)<1e-14 && fabs(out[2]-.6)<1e-14);
 if(remote)return 0;
 /* Signed corrections must be zero, never a temperature, at an empty node. */
 project(&E,values,out,0,NULL);assert(out[3]==0);
 /* End-of-heat diagnostics retain the latest grid, not the old timestep. */
 project(&E,values,out,0,stage);assert(out[3]==.9);
 /* TA projects in place: preservation must precede clearing accumulators. */
 project(&E,values,E.T[1],1,E.T[1]);assert(E.T[1][3]==.7);
 /* A subsequent grid heat update must survive the next projection. */
 E.T[1][3]+=.02;project(&E,values,E.T[1],1,E.T[1]);assert(fabs(E.T[1][3]-.72)<1e-14);
 /* Independent fallback alias, used by a diagnostic without boundary fixing. */
 project(&E,values,stage,0,stage);assert(stage[3]==.9);
 return 0;
}
'''
with tempfile.TemporaryDirectory() as folder:
    p = Path(folder)
    (p / 'test.c').write_text(stub + valid + project + main)
    subprocess.run(['cc', '-std=gnu11', '-Wall', '-Wextra', '-Wno-unused-parameter',
                    str(p / 'test.c'), '-lm', '-o', str(p / 'test')], check=True)
    for mode in ['', 'fixed', 'empty', 'remote', 'bad_fallback', 'nan', 'negative']:
        got = subprocess.run([str(p / 'test'), mode], capture_output=True, text=True)
        assert got.returncode == (72 if mode in ['bad_fallback', 'nan', 'negative'] else 0), (mode, got.stderr)
        if mode == '':
            assert 'zero_free_nodes=1' in got.stdout
        if mode == 'remote':
            assert 'zero_free_nodes=0' in got.stdout
    print('PASS: actual C projection, in-place TA, signed increments, heated-grid fallback, MPI-shared coverage, Dirichlet and invalid-input checks')
