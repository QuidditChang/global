"""Production phase split against independent finite-difference derivatives."""
import ctypes as C
import math,subprocess,tempfile,unittest
import os,shutil
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
  phase_source=(ROOT/'lib/Phase_change.c').read_text()
  phase='\n'.join(function(phase_source,signature) for signature in (
   'int phase_change_count(', 'double phase_change_entropy(',
   'int phase_change_thermal_active(', 'void phase_change_state('))
  pices_source=(ROOT/'lib/Pices.c').read_text()
  source='''#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "phase_change.h"
void pices_fail(struct All_variables *E,const char *why) {abort();}
'''+phase+function(pices_source,'static void phase_coefficients(')+function(pices_source,'static double phase_fraction_change(')+'''
void evaluate_phase(double t,double velocity,double entropy,int phase_index,
 int ppv_enabled,int latent_enabled,double *out) {
 struct All_variables state,*E=&state;struct IEN ien[2];struct Shape_function_dx dx;
 double rho[3]={0,1,2},gravity[3]={0,1,1};float v[9];int a;
 memset(E,0,sizeof(*E));memset(&dx,0,sizeof(dx));E->lmesh.noz=2;
 E->ien[1]=ien;E->refstate.rho=rho;E->refstate.gravity=gravity;E->sphere.cap[1].V[3]=isnan(velocity)?NULL:v;
 E->sphere.ro=1;E->control.surface_temp=.125;
 E->control.pPv_enabled=ppv_enabled;E->control.pPv_latent_enabled=latent_enabled;
 E->control.pPv_k_ratio=1;
 E->control.phase[phase_index].depth=.25;E->control.phase[phase_index].transT=.5;
 E->control.phase[phase_index].clapeyron=.25;E->control.phase[phase_index].inv_width=4;
 E->control.phase[phase_index].entropy_jump=entropy;
 for(a=1;a<=8;a++) {ien[1].node[a]=a;E->N.vpt[GNVINDEX(a,1)]=.125;v[a]=velocity;
 dx.vpt[GNVXINDEX(2,a,1)]=(a%2?-1:1)/4.;}
 phase_coefficients(E,1,1,&dx,2,t,1.5,1.25,out,isnan(velocity)?NULL:out+1,out+2);
}
void evaluate(double t,double velocity,double entropy,double *out) {
 evaluate_phase(t,velocity,entropy,0,0,0,out);
}
double adaptation(int ppv_enabled,int latent_enabled,double ratio,int phase_index,int gp) {
 struct All_variables state,*E=&state;
 double current[16*PHASE_TRANSITIONS]={0},previous[16*PHASE_TRANSITIONS]={0};
 memset(E,0,sizeof(*E));E->control.pPv_enabled=ppv_enabled;
 E->control.pPv_latent_enabled=latent_enabled;E->control.pPv_k_ratio=ratio;
 current[gp*PHASE_TRANSITIONS+phase_index]=.1;
 return phase_fraction_change(E,current,previous,16);
}
'''
  (d/'phase.c').write_text(source)
  compiler=shutil.which(os.environ.get('MPICC','mpicc'))
  if compiler is None:raise RuntimeError('mpicc is required to test the production C block')
  subprocess.run([compiler,'-std=gnu99','-shared','-fPIC','-I'+str(ROOT/'lib'),str(d/'phase.c'),'-lm','-o',str(d/'phase.so')],check=True)
  cls.lib=C.CDLL(str(d/'phase.so'));cls.lib.evaluate.argtypes=[C.c_double]*3+[C.POINTER(C.c_double)]
  cls.lib.evaluate_phase.argtypes=[C.c_double]*3+[C.c_int]*3+[C.POINTER(C.c_double)]
  cls.lib.adaptation.argtypes=[C.c_int,C.c_int,C.c_double,C.c_int,C.c_int]
  cls.lib.adaptation.restype=C.c_double
 @classmethod
 def tearDownClass(cls):cls.tmp.cleanup()
 def test_capacity_pressure_and_velocity_sign(self):
  # At r=.5, rho*g=1.5 and its radial derivative is 1.
  def fraction(t,r):return .5*(1+math.tanh(4*((1-r-.25)*(1.5+r-.5)-.25*(t-.5))))
  for t in [.1,.5,.9]:
   h=1e-6;dt=(fraction(t+h,.5)-fraction(t-h,.5))/(2*h);dr=(fraction(t,.5+h)-fraction(t,.5-h))/(2*h)
   for v in [-2,0,3]:
    out=(C.c_double*6)();self.lib.evaluate(t,v,-.125,out)
    self.assertAlmostEqual(out[0],1.5*(1.25+(t+.125)*(-.125)*dt),places=9)
    self.assertAlmostEqual(out[1],-1.5*(t+.125)*(-.125)*dr*v,places=9)
    self.assertAlmostEqual(out[2],fraction(t,.5),places=14)
 def test_capacity_before_velocity_initialization(self):
  out=(C.c_double*6)();self.lib.evaluate(.5,float("nan"),-.125,out)
  self.assertTrue(math.isfinite(out[0]) and out[0]>0)
 def test_zero_entropy_limit(self):
  out=(C.c_double*6)();self.lib.evaluate(.5,3,0,out)
  self.assertEqual(out[0],1.5*1.25);self.assertEqual(out[1],0)
 def test_pPv_latent_switch(self):
  legacy=(C.c_double*6)();self.lib.evaluate(.5,3,-.125,legacy)
  for enabled,latent in [(0,1),(1,0),(1,1)]:
   out=(C.c_double*6)();self.lib.evaluate_phase(.5,3,-.125,3,enabled,latent,out)
   if enabled and latent:
    self.assertEqual(out[0],legacy[0]);self.assertEqual(out[1],legacy[1])
   else:
    self.assertEqual(out[0],1.5*1.25);self.assertEqual(out[1],0)
   self.assertEqual(out[5],legacy[2] if enabled else 0)
 def test_phase_fraction_adaptation_uses_active_thermal_phases_and_four_slot_stride(self):
  for enabled,latent,ratio in [(0,1,2),(1,0,1),(1,1,1),(1,0,2)]:
   for gp in [8,15]:
    self.assertEqual(self.lib.adaptation(enabled,latent,ratio,0,gp),2)
    expected=2 if enabled and (latent or ratio!=1) else 0
    self.assertEqual(self.lib.adaptation(enabled,latent,ratio,3,gp),expected)
  self.assertEqual(self.lib.adaptation(1,1,2,3,7),0)
if __name__=='__main__':unittest.main()
