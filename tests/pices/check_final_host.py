"""Replay rank 287 / step 2240 weights through the actual final-host selector.

The geometry stub replays the measured two wedges, rather than recomputing
compiler/platform-dependent rounding from rounded diagnostic coordinates.
"""
from pathlib import Path
import subprocess
import tempfile

src = (Path(__file__).resolve().parents[2] / 'lib/Full_tracer_advection.c').read_text()
weights = src[src.index('static int pices_simplex('):src.index('void full_get_velocity(')]
selector = src[src.index('static int pices_validated_host('):src.index('int full_iget_element(')]
# All successful locator branches converge on this final validation.
body = src[src.index('int full_iget_element('):]
final = body[body.index(' foundit:'):body.index('return ifinal_iel;')]
assert 'pices_validated_host(E,j,ifinal_iel,theta,phi,rad)' in final
harness = r'''
#include <math.h>
#include <assert.h>
struct All_variables {
 struct {int nel;} lmesh;
 struct {int node[9];} ien[2][3];
 double sx[2][4][9];
};
static int calls;
static void spherical_to_uv(struct All_variables *E,int c,double t,double p,double *u,double *v){
 *u=-.036595286727281705;*v=-.11836816391782183;
}
static void get_2dshape(struct All_variables *E,int c,int e,double u,double v,int w,double *s){
 const double recorded[2][2][3]={
  {{.51829220596224701,-.51923967890672551,1.0009474729444729},
   {-1.0008349704548891e-10,.48087542786221427,.5191245722378709}},
  {{.51912457213763297,.48087542776157616,1.0080114520860661e-10},
   {.99913884631918926,.48164346134220715,-.48078230766138791}}};
 calls++;
 for(int i=1;i<=3;i++) s[i]=recorded[e-1][w-1][i-1];
}
''' + weights + selector + r'''
int main(void){
 struct All_variables E={0};int n[9];double w[9],sum;
 const double t=3.0891580913534531,p=4.7183981175314695,r=.96883360885710978;
 E.lmesh.nel=2;
 for(int a=1;a<=8;a++){
  E.ien[1][1].node[a]=E.ien[1][2].node[a]=a;
  E.sx[1][3][a]=a<=4?.96799600124359131:.97028201818466187;
 }
 /* Element 1 replays rejected host 24563; 2 replays valid host 25587. */
 assert(!pices_element_weights(&E,1,1,t,p,r,n,w));
 assert(pices_validated_host(&E,1,1,t,p,r)==2);
 assert(pices_element_weights(&E,1,2,t,p,r,n,w));
 sum=0;for(int a=1;a<=8;a++){assert(w[a]>=0 && w[a]<=1);sum+=w[a];}
 assert(fabs(sum-1)<1e-14);
 calls=0;assert(pices_validated_host(&E,1,2,t,p,r)==2);assert(calls==1);
 assert(pices_validated_host(&E,1,1,t,p,.99)==-1);
 E.lmesh.nel=1;assert(pices_validated_host(&E,1,1,t,p,r)==-1);
 /* Preserve tiny-roundoff clipping at the original tolerance, reject beyond. */
 double s[4]={0,-.5e-10,.4,.6+.5e-10};assert(pices_simplex(s,3));
 s[1]=-1.0008349704548891e-10;s[2]=.4;s[3]=.6-s[1];
 assert(!pices_simplex(s,3));
 return 0;
}
'''
with tempfile.TemporaryDirectory() as tmp:
    p = Path(tmp)
    (p / 'test.c').write_text(harness)
    subprocess.run(['cc', '-std=gnu11', str(p / 'test.c'), '-lm', '-o', str(p / 'test')], check=True)
    subprocess.run([str(p / 'test')], check=True)
print('PASS: recorded RH20 final host rejection/recovery, fast path, radial rejection, unchanged tolerance')
