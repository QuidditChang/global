"""Compile the failure-only snapshot with controlled tangential/radial errors."""
from pathlib import Path
import re,subprocess,tempfile,unittest
SOURCE=Path(__file__).resolve().parents[2]/'lib/Full_tracer_advection.c'
class Diagnostic(unittest.TestCase):
 def test_snapshot_values_invalid_host_and_no_state_mutation(self):
  s=SOURCE.read_text();body=s.split('static void pices_weight_geometry(',1)[1].split('/* PICES uses the same',1)[0]
  stub=r'''
#include <stdio.h>
#include <string.h>
struct All_variables {
 struct {int nel;} lmesh; struct {int me;} parallel;
 struct {double basicq[2][6][2],shape_coefs[2][3][10][2],UV[2][3][9];} trace;
 struct {int node[9];} ien[2][2];
 double sx[2][4][9],x[2][4][9];
};
static void spherical_to_uv(struct All_variables *E,int c,double t,double p,double *u,double *v){*u=t;*v=p;}
static void get_2dshape(struct All_variables *E,int c,int e,double u,double v,int w,double *s){
 for(int a=1;a<=3;a++)s[a]=E->trace.shape_coefs[c][w][3*a-2][e]+E->trace.shape_coefs[c][w][3*a-1][e]*u+E->trace.shape_coefs[c][w][3*a][e]*v;
}
'''
  main=r'''
int main(void){
 struct All_variables E={0},saved;
 E.lmesh.nel=1;E.parallel.me=148;E.trace.basicq[1][2][1]=.699999999;
 for(int a=1;a<=8;a++){E.ien[1][1].node[a]=a;E.sx[1][3][a]=a<=4?.7:.8;}
 for(int w=1;w<=2;w++){E.trace.shape_coefs[1][w][1][1]=-1e-10;E.trace.shape_coefs[1][w][4][1]=.25;E.trace.shape_coefs[1][w][7][1]=.7500000001;}
 saved=E;pices_weight_geometry(&E,1,1,1,"initial");pices_weight_geometry(&E,1,1,-99,"relocated");
 return memcmp(&E,&saved,sizeof(E))!=0;
}
'''
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);(p/'test.c').write_text(stub+'static void pices_weight_geometry('+body+main)
   subprocess.run(['cc','-std=c99',str(p/'test.c'),'-o',str(p/'test')],check=True)
   run=subprocess.run([str(p/'test')],capture_output=True,text=True,check=True)
  out=run.stderr
  self.assertIn('candidate=relocated host=-99 valid=0',out)
  self.assertEqual(out.count('PICES_WEIGHT_NODE '),8)
  self.assertEqual(out.count('PICES_WEIGHT_SHAPE '),2)
  self.assertEqual(out.count('PICES_WEIGHT_COEF '),6)
  radial=[float(x) for x in re.search(r'radial_shape=([^ ]+)',out)[1].split(',')]
  self.assertLess(radial[1],-1e-12)
  tangential=[float(x) for x in re.search(r'tangential=([^ ]+)',out)[1].split(',')]
  self.assertEqual(tangential[0],-1e-10)
  self.assertAlmostEqual(sum(tangential),1.)
if __name__=='__main__':unittest.main()
