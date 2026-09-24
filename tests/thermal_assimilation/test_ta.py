"""Compile and exercise production TA formulas and the real one-column BC path."""
import math
from pathlib import Path
import shutil, subprocess, tempfile, unittest, re
ROOT=Path(__file__).resolve().parents[2]

def function(s,signature):
 start=s.rindex(signature);i=s.index('{',start);j=i+1;n=1
 masked=re.sub(r'/\*.*?\*/|//[^\n]*|"(?:\\.|[^"\\])*"',lambda m:' '*len(m.group()),s,flags=re.S)
 while n:
  n+=(masked[j]=='{')-(masked[j]=='}');j+=1
 return s[start:j]

class ThermalAssimilationTest(unittest.TestCase):
 def test_production_column_and_defaults(self):
  src=(ROOT/'lib/Lith_age.c').read_text()
  funcs=['static float effective_plate_age_nd(', 'static double lith_age_surface_anomaly(', 'static double lith_age_target_temperature(', 'static double lith_age_old_temperature_weight(', 'static double lith_age_asml_thickness(', 'static void validate_lith_age_asml(', 'void lith_age_temperature_bound_adj(', 'void lith_age_conform_tbc(', 'void assimilate_lith_conform_bcs(', 'void lith_age_update_tbc(', 'void lith_age_construct_tic(']
  prefix=r'''
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "element_definitions.h"
#include "global_defs.h"
static int reads=0;
void parallel_process_termination(void) { abort(); }
float find_age_in_MY(struct All_variables *E) { return 100.; }
void temperatures_conform_bcs(struct All_variables *E) { (void)E; }
void reader(struct All_variables *E,int output) { reads++; }
'''
  main=r'''
static struct All_variables state;
int main(void) {
 struct All_variables *E=&state;
 int k; double h,expected,old,age,w;
 E->control.lith_age=1; E->control.lith_age_time=1; E->control.lith_age_asml=1;
 E->control.eba_formulation=1; E->control.lith_age_depth=.04;
 E->control.max_plate_age_Ma=70.; E->control.lith_age_mantle_temp=.395917143;
 E->sphere.ro=1.; E->sphere.ri=.546225; E->sphere.caps_per_proc=1;
 E->data.scalet=1286263.5; E->data.Ttop=300.;E->data.Tbottom=3700.;
 E->data.ref_temperature=3400.; E->data.radius_km=6371.;
 E->control.TBCtopval=0.;E->control.TBCbotval=1.;
 E->mesh.gridmax=1; E->mesh.levmax=1; E->mesh.noz=5;
 E->mesh.nox=E->mesh.noy=1; E->lmesh.nox=E->lmesh.noy=1;
 E->lmesh.noz=E->lmesh.nno=5; E->lmesh.nxs=E->lmesh.nys=E->lmesh.nzs=1;
 E->solver.lith_age_read_files=reader;
 E->age_t=calloc(2,sizeof(float));E->flag_depth2=calloc(2,sizeof(float));
 E->refstate.Tref=calloc(6,sizeof(double));
 E->T[1]=calloc(6,sizeof(float));E->assim_delta_T[1]=calloc(6,sizeof(double));
 E->node[1]=calloc(6,sizeof(unsigned int));
 for(k=1;k<=3;k++) {
  E->sx[1][k]=calloc(6,sizeof(double));
  E->sphere.cap[1].TB[k]=calloc(6,sizeof(float));
 }
 E->refstate.has_temperature=1;E->refstate.temperature_surface=.388205;
 E->refstate.temperature_cmb=.67;
 for(k=1;k<=5;k++) E->refstate.Tref[k]=.388205+.01*(5-k);
 E->flag_depth2[1]=.2; age=70./E->data.scalet;E->age_t[1]=age;
 h=E->control.lith_age_depth;
 E->sx[1][3][1]=E->sphere.ri;E->sx[1][3][2]=1.-h;
 E->sx[1][3][3]=1.-h/2.;E->sx[1][3][4]=1.-h/4.;E->sx[1][3][5]=1.;
 validate_lith_age_asml(E);
 assert(lith_age_old_temperature_weight(E,0,h)==0.);
 assert(lith_age_old_temperature_weight(E,h,h)==1.);
 assert(fabs(lith_age_old_temperature_weight(E,h/2,h)-.5)<1e-12);
 assert(fabs(lith_age_surface_anomaly(E,h/2,age)-
    E->refstate.temperature_surface*erfc(h/4/sqrt((float)age)))<1e-12);
 assert(fabs(lith_age_surface_anomaly(E,h/2,100./E->data.scalet)-
    lith_age_surface_anomaly(E,h/2,age))<1e-10);
 assert(lith_age_surface_anomaly(E,h/2,0.)<1e-40);
 E->control.lith_age_asml_exp=3.;
 w=(exp(-1.5)-exp(-3.))/(1-exp(-3.));
 assert(fabs((1-lith_age_old_temperature_weight(E,h/2,h))-w)<1e-12);
 assert(lith_age_old_temperature_weight(E,0,h)==0.);
 assert(lith_age_old_temperature_weight(E,h,h)==1.);
 assert(lith_age_old_temperature_weight(E,h/4,h)<lith_age_old_temperature_weight(E,h/2,h));
 /* Same trench support, but initialization applies full HSC anomaly. */
 lith_age_construct_tic(E);
 expected=E->refstate.Tref[3]-lith_age_surface_anomaly(E,h/2,age);
 assert(fabs(E->T[1][3]-expected)<1e-7);
 assert(fabs(E->T[1][2]-(E->refstate.Tref[2]-lith_age_surface_anomaly(E,h,age)))<1e-7);
 E->flag_depth2[1]=.05;
 lith_age_construct_tic(E);
 assert(fabs(E->T[1][3]-E->refstate.Tref[3])<1e-7);
 E->flag_depth2[1]=.08;
 assert(fabs(lith_age_asml_thickness(E,1)/h-(.003+.997*.5))<1e-6);
 E->flag_depth2[1]=.2;
 expected=E->refstate.Tref[3]-lith_age_surface_anomaly(E,h/2,age);
 E->monitor.solution_cycles=1;
 lith_age_conform_tbc(E);
 assert(fabs(E->sphere.cap[1].TB[3][3]-expected)<1e-7);
 lith_age_temperature_bound_adj(E,1);
 assert((E->node[1][2]&(TBX|TBY|TBZ))==0); /* zero-weight bottom remains free */
 assert(E->node[1][3]&TBZ);
 E->T[1][3]=.6;old=E->T[1][3];
 assimilate_lith_conform_bcs(E);
 assert(fabs(E->T[1][3]-((1-w)*old+w*expected))<1e-7);
 assert(fabs(E->assim_delta_T[1][3]-(E->T[1][3]-old))<1e-12);
 /* Legacy depth weights and target retained when switch off. */
 E->control.lith_age_asml=0;
 assert(fabs(lith_age_old_temperature_weight(E,h,h)-.5)<1e-12);
 lith_age_conform_tbc(E);
 expected=E->control.lith_age_mantle_temp*erf(h/4/sqrt((float)age));
 assert(fabs(E->sphere.cap[1].TB[3][3]-expected)<1e-7);
 /* Runtime signed trench mask reduces depth using unchanged .003 factor. */
 E->control.lith_age_asml=1;E->flag_depth2[1]=-.2;
 lith_age_temperature_bound_adj(E,1);
 assert(!(E->node[1][3]&TBZ));assert(E->node[1][5]&TBZ);
 /* Setup callback must not access unloaded reference values, including restart. */
 E->refstate.has_temperature=0;free(E->refstate.Tref);E->refstate.Tref=NULL;
 lith_age_conform_tbc(E);
 assert(reads>=3);
 puts("PASS: target, cap, initialization parity, taper, live flags, legacy, deltaT, pre-reference callback");
 return 0;
}
'''
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);c=p/'test.c';exe=p/'test';c.write_text(prefix+'\n'.join(function(src,x) for x in funcs)+main)
   subprocess.run([shutil.which('mpicc'),'-std=gnu99','-Wno-deprecated-non-prototype','-I'+str(ROOT/'lib'),str(c),'-lm','-o',str(exe)],check=True,text=True)
   r=subprocess.run([str(exe)],check=True,capture_output=True,text=True);print(r.stdout)
 def test_pyre_and_c_defaults(self):
  self.assertIn('input_int("lith_age_asml",&(E->control.lith_age_asml),"0",m)',(ROOT/'lib/Lith_age.c').read_text())
  self.assertIn('pyre.inventory.int("lith_age_asml", default=0)',(ROOT/'CitcomS/Components/Param.py').read_text())
  self.assertIn('getIntProperty(properties, "lith_age_asml", E->control.lith_age_asml, fp)',(ROOT/'module/setProperties.c').read_text())

if __name__=='__main__':unittest.main()
