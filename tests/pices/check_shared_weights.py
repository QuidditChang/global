"""Compile shared geometry helpers; replay step 74 and reject true outsiders."""
from pathlib import Path
import subprocess,tempfile,unittest
ROOT=Path(__file__).resolve().parents[2]
class SharedWeights(unittest.TestCase):
 def test_production_face_and_simplex_bounds(self):
  source=(ROOT/'lib/Full_tracer_advection.c').read_text()
  helper=source[source.index('static int pices_simplex('):source.index('void full_get_velocity(')]
  stub=r'''
#include <math.h>
#include <assert.h>
struct All_variables {
 struct {int nel;} lmesh;
 struct {int node[9];} ien[2][2];
 double sx[2][4][9];
};
static void spherical_to_uv(struct All_variables *E,int c,double t,double p,double *u,double *v){*u=t;*v=p;}
static void get_2dshape(struct All_variables *E,int c,int e,double u,double v,int w,double *s){
 const double abc[2][3][3]={
 {{-6.8579748751057954,69.127720382916806,85.109010065372047},
 {.98286530473847034,-138.42815682912706,.2126454205561904},
 {6.875109570367326,69.300436446210256,-85.321655485928233}},
 {{-5.8765082001421414,-69.103451426105778,85.321352888802025},
 {7.8593393864938719,-69.319900120580911,-85.108714849845995},
 {-.9828311863517325,138.42335154668669,-.21263803895602443}}};
 for(int a=1;a<=3;a++)s[a]=abc[w-1][a-1][0]+abc[w-1][a-1][1]*u+abc[w-1][a-1][2]*v;
}
'''
  main=r'''
int main(void){
 struct All_variables E={0};int nodes[9];double weights[9],s[4],sum;
 double u=.0094648666971343977,v=.076540782993689369,r=.65870817770374723;
 E.lmesh.nel=1;
 for(int a=1;a<=8;a++){E.ien[1][1].node[a]=a;E.sx[1][3][a]=a<=4?.65402901172637939:.66538697481155396;}
 assert(pices_column_shape(&E,1,1,u,v,s)==2);
 assert(pices_element_weights(&E,1,1,u,v,r,nodes,weights));
 sum=0;for(int a=1;a<=8;a++){assert(weights[a]>=0&&weights[a]<=1);sum+=weights[a];}
 assert(fabs(sum-1)<1e-14);assert(weights[1]==0 && weights[5]==0);
 assert(fabs(weights[3]-.40512057762307241)<2e-11);
 assert(pices_element_weights(&E,1,1,u,v,E.sx[1][3][1],nodes,weights));
 assert(pices_element_weights(&E,1,1,u,v,E.sx[1][3][5],nodes,weights));
 assert(!pices_element_weights(&E,1,1,u,v,.65,nodes,weights));
 assert(!pices_element_weights(&E,1,1,u,v-1e-5,r,nodes,weights));
 assert(!pices_element_weights(&E,1,-99,u,v,r,nodes,weights));
 for(int k=-10;k<=10;k++) {
   double x=k*1e-11,a[4]={0,x,.3,.7-x};
   assert(pices_simplex(a,3));assert(a[1]>=0);assert(fabs(a[1]+a[2]+a[3]-1)<1e-14);
 }
 double outside[4]={0,-1e-8,.3,.70000001};assert(!pices_simplex(outside,3));
 double nanval[3]={0,NAN,.5};assert(!pices_simplex(nanval,2));
 double diagonal[4]={0,.4,0,.6};assert(pices_simplex(diagonal,3));
 return 0;
}
'''
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);(p/'test.c').write_text(stub+helper+main)
   subprocess.run(['cc','-std=c99',str(p/'test.c'),'-lm','-o',str(p/'test')],check=True)
   subprocess.run([str(p/'test')],check=True)
if __name__=='__main__':unittest.main()
