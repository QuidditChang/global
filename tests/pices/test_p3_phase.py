"""Production phase split against independent finite-difference derivatives."""
import ctypes as C
import math,subprocess,tempfile,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
def function(s,signature):
 a=s.index(signature);b=s.index('{',a);depth=1;i=b+1
 while depth:depth+=(s[i]=='{')-(s[i]=='}');i+=1
 return s[a:i]
class PhaseSplit(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.tmp=tempfile.TemporaryDirectory();d=Path(cls.tmp.name)
  source='''#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "element_definitions.h"
#include "global_defs.h"
void pices_fail(struct All_variables *E,const char *why) {abort();}
'''+function((ROOT/'lib/Phase_change.c').read_text(),'void phase_change_state(')+function((ROOT/'lib/Pices.c').read_text(),'static void phase_coefficients(')+'''
void evaluate(double t,double velocity,double entropy,double *out) {
 struct All_variables state,*E=&state;struct IEN ien[2];struct Shape_function_dx dx;
 double rho[3]={0,1,2},gravity[3]={0,1,1};float v[9];int a;
 memset(E,0,sizeof(*E));memset(&dx,0,sizeof(dx));E->lmesh.noz=2;
 E->ien[1]=ien;E->refstate.rho=rho;E->refstate.gravity=gravity;E->sphere.cap[1].V[3]=isnan(velocity)?NULL:v;
 E->sphere.ro=1;E->control.surface_temp=.125;
 E->control.phase[0].depth=.25;E->control.phase[0].transT=.5;
 E->control.phase[0].clapeyron=.25;E->control.phase[0].inv_width=4;
 E->control.phase[0].entropy_jump=entropy;
 for(a=1;a<=8;a++) {ien[1].node[a]=a;E->N.vpt[GNVINDEX(a,1)]=.125;v[a]=velocity;
 dx.vpt[GNVXINDEX(2,a,1)]=(a%2?-1:1)/4.;}
 phase_coefficients(E,1,1,&dx,2,t,1.5,1.25,out,isnan(velocity)?NULL:out+1,out+2);
}
'''
  (d/'phase.c').write_text(source)
  subprocess.run(['/usr/local/bin/mpicc','-std=gnu99','-shared','-fPIC','-I'+str(ROOT/'lib'),str(d/'phase.c'),'-lm','-o',str(d/'phase.so')],check=True)
  cls.lib=C.CDLL(str(d/'phase.so'));cls.lib.evaluate.argtypes=[C.c_double]*3+[C.POINTER(C.c_double)]
 @classmethod
 def tearDownClass(cls):cls.tmp.cleanup()
 def test_capacity_pressure_and_velocity_sign(self):
  # At r=.5, rho*g=1.5 and its radial derivative is 1.
  def fraction(t,r):return .5*(1+math.tanh(4*((1-r-.25)*(1.5+r-.5)-.25*(t-.5))))
  for t in [.1,.5,.9]:
   h=1e-6;dt=(fraction(t+h,.5)-fraction(t-h,.5))/(2*h);dr=(fraction(t,.5+h)-fraction(t,.5-h))/(2*h)
   for v in [-2,0,3]:
    out=(C.c_double*5)();self.lib.evaluate(t,v,-.125,out)
    self.assertAlmostEqual(out[0],1.5*(1.25+(t+.125)*(-.125)*dt),places=9)
    self.assertAlmostEqual(out[1],-1.5*(t+.125)*(-.125)*dr*v,places=9)
    self.assertAlmostEqual(out[2],fraction(t,.5),places=14)
 def test_capacity_before_velocity_initialization(self):
  out=(C.c_double*5)();self.lib.evaluate(.5,float("nan"),-.125,out)
  self.assertTrue(math.isfinite(out[0]) and out[0]>0)
 def test_zero_entropy_limit(self):
  out=(C.c_double*5)();self.lib.evaluate(.5,3,0,out)
  self.assertEqual(out[0],1.5*1.25);self.assertEqual(out[1],0)
if __name__=='__main__':unittest.main()
