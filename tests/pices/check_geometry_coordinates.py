"""Recorded RH20 face: locator and interpolation must make identical decisions."""
from pathlib import Path
import subprocess,tempfile
src=(Path(__file__).resolve().parents[2]/'lib/Full_tracer_advection.c').read_text()
helpers=src[src.index('static int pices_simplex('):src.index('void full_get_velocity(')]
start=src.index('static int icheck_element_column(',src.index('/********  ICHECK ELEMENT COLUMN'))
body=src[start:src.index('/* surface coords of element nodes */',start)]
assert 'atan2' not in body
assert 'double theta, double phi, double rad)' in body
# Compile the actual PICES branch of the locator with the actual weight helper.
branch=body[body.index('if(E->pices.enabled) {'):]
branch=branch[:branch.index('\n    }')+6]
harness=r'''
#include <math.h>
#include <assert.h>
struct All_variables {
 struct {int enabled;} pices;
 struct {int nel;} lmesh;
 struct {int node[9];} ien[2][2];
 double sx[2][4][9];
};
static void spherical_to_uv(struct All_variables *E,int c,double t,double p,double *u,double *v){*u=t;*v=p;}
/* Replay measured shape weights. Equivalent longitude representations used
 * by the old locator can round to opposite sides of the tolerance. */
static void get_2dshape(struct All_variables *E,int c,int e,double u,double v,int w,double *s){
 s[1]=.11350080352017855;s[3]=v>=0?-1.000017846308765e-10:-9.99982319171977e-11;
 s[2]=1-s[1]-s[3];
}
'''+helpers+'''
static int locator(struct All_variables *E,int j,int nel,double x,double y,double z,double theta,double phi,double rad){
'''+branch+r'''
return 0;
}
int main(void){
 struct All_variables E={0};int n[9];double w[9];
 double t=.5880185089527552,p=3.5875096924756065,r=.965150215512386;
 double x=-.48302976747254578,y=-.23090236318601698,z=.80304500553464397;
 E.pices.enabled=1;E.lmesh.nel=1;
 for(int a=1;a<=8;a++){E.ien[1][1].node[a]=a;E.sx[1][3][a]=a<=4?.96342402696609497:.96570998430252075;}
 assert(pices_column_shape(&E,1,1,atan2(hypot(x,y),z),atan2(y,x),w));
 assert(!pices_element_weights(&E,1,1,t,p,r,n,w));
 assert(!locator(&E,1,1,x,y,z,t,p,r));
 p=atan2(y,x);
 assert(locator(&E,1,1,x,y,z,t,p,r));
 assert(pices_element_weights(&E,1,1,t,p,r,n,w));
 return 0;
}
'''
with tempfile.TemporaryDirectory() as d:
 p=Path(d);(p/'test.c').write_text(harness)
 subprocess.run(['cc','-std=gnu11',str(p/'test.c'),'-lm','-o',str(p/'test')],check=True)
 subprocess.run([str(p/'test')],check=True)
print('PASS: measured face weights, unified angular inputs, unchanged tolerance')
