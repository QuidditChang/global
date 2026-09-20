"""Execute production Qvis heat source and refstate reader, not a formula copy."""
import ctypes
import math
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]

def function(source, signature):
    start = source.index(signature)
    i = source.index('{', start) + 1
    depth = 1
    while depth:
        depth += (source[i] == '{') - (source[i] == '}')
        i += 1
    return source[start:i]

class QvisTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.directory = Path(cls.tmp.name)
        energy=(ROOT/'lib/Advection_diffusion.c').read_text()
        material=(ROOT/'lib/Material_properties.c').read_text()
        prefix='''
#include <math.h>
#include <ctype.h>
#include <float.h>
#include <string.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "qvis_limiter.h"
static double fixture_s;
void strain_rate_2_inv(struct All_variables *E,int m,float *s,int root) {s[1]=fixture_s;}
void parallel_process_termination(void) { exit(19); }
'''
        thermal=prefix+function(energy,'static void process_visc_heating(')+'''
void heat(int mode,double eta,double shear,double pressure,double cohesion,double phi,double *out) {
 struct All_variables state={0},*E=&state;
 float visc[9]; double p[3]={0,pressure,pressure}; int i;
 E->lmesh.nel=E->lmesh.elz=1;
 E->control.eba_formulation=1; E->control.qvis_mode=mode;
 E->control.qvis_cohesion_pa=cohesion; E->control.qvis_friction_angle_rad=phi;
 E->control.disptn_number=E->control.Atemp=1;
 E->data.ref_viscosity=1e21; E->data.kappa0=1e-6; E->data.radius_km=6371;
 E->refstate.lithostatic_pressure_pa=p;
 for(i=1;i<=8;i++) visc[i]=eta;
 E->EVi[1]=visc;
 /* Simple shear du_x/dy=gamma: S=gamma^2. */
 fixture_s=shear*shear;
 double used[2],raw[2],capped[2];
 process_visc_heating(E,1,used,raw,capped);
 out[0]=raw[1];out[1]=capped[1];out[2]=used[1];
}
'''
        f=cls.directory/'heat.c';f.write_text(thermal)
        cc=os.environ.get('MPICC','mpicc')
        flags=[cc,'-std=gnu99','-Wno-deprecated-non-prototype','-I'+str(ROOT/'lib')]
        subprocess.run(flags+['-shared','-fPIC',str(f),'-o',str(cls.directory/'heat.so')],check=True)
        cls.lib=ctypes.CDLL(str(cls.directory/'heat.so'))
        cls.lib.heat.argtypes=[ctypes.c_int]+[ctypes.c_double]*5+[ctypes.POINTER(ctypes.c_double)]
        reader=prefix+function(material,'static int read_refstate_data_line(FILE *fp, char *buffer, int length)\n{')+function(material,'static int read_eba_row(')+function(material,'static void read_refstate(struct All_variables *E)\n{')+'''
int main(int argc,char **argv) {
 struct All_variables state={0},*E=&state;
 double fields[10][8]={{0}};
 E->control.eba_formulation=1; E->mesh.noz=5;
 E->lmesh.noz=3; E->lmesh.nzs=atoi(argv[2]);
 E->refstate.rho=fields[0]; E->refstate.gravity=fields[1];
 E->refstate.temperature=fields[2]; E->refstate.thermal_expansivity=fields[3];
 E->refstate.heat_capacity=fields[4]; E->refstate.dis=fields[5];
 E->refstate.gamma_eff=fields[6]; E->refstate.beta_ala=fields[7];
 E->refstate.lithostatic_pressure_pa=fields[8];
 strcpy(E->refstate.filename,argv[1]);read_refstate(E);
 printf("%d %.17g %.17g\\n",E->refstate.has_lithostatic_pressure,fields[8][1],fields[8][3]);
 return 0;
}
'''
        f=cls.directory/'reader.c';f.write_text(reader)
        subprocess.run(flags+[str(f),'-o',str(cls.directory/'reader')],check=True)

    @classmethod
    def tearDownClass(cls): cls.tmp.cleanup()

    def heat(self,mode,eta=100.,shear=10000.,pressure=0.,cohesion=1e7,phi=.085):
        out=(ctypes.c_double*3)()
        self.lib.heat(mode,eta,shear,pressure,cohesion,phi,out)
        return list(out)

    def test_simple_shear_yield_and_factor_of_two(self):
        raw,cap,used=self.heat(2)
        sigma=6e7*math.cos(.085)/(math.sqrt(3)*(3+math.sin(.085)))
        scale=ctypes.c_float(1e21).value*ctypes.c_float(1e-6).value/6371000.**2
        self.assertEqual(raw,100.*10000.**2)
        self.assertAlmostEqual(cap/(sigma/scale*10000.),1.,places=12)
        self.assertEqual(cap,used)

    def test_diagnostic_does_not_change_heat(self):
        raw,cap,used=self.heat(1)
        self.assertLess(cap,raw)
        self.assertEqual(used,raw)
        self.assertEqual(self.heat(0),[raw]*3)

    def test_no_yield_and_zero_strain(self):
        a=self.heat(2,shear=1.)
        self.assertEqual(a,[a[0]]*3)
        self.assertEqual(self.heat(2,shear=0),[0.]*3)

    def test_pressure_strengthens_and_zero_strength(self):
        self.assertGreater(self.heat(2,pressure=1e9)[1],self.heat(2)[1])
        self.assertEqual(self.heat(2,cohesion=0,phi=0)[1:], [0.,0.])

    def read(self,pressures,start=1,extra=''):
        f=self.directory/'refstate.txt'
        f.write_text(''.join('1 1 0.4 1 1'+(' '+str(p) if p is not None else '')+extra+'\n' for p in pressures))
        return subprocess.run([str(self.directory/'reader'),str(f),str(start)],capture_output=True,text=True)

    def test_radial_mpi_slices_use_global_pressure(self):
        self.assertEqual(self.read([40,30,20,10,0],1).stdout.strip(),'1 40 20')
        self.assertEqual(self.read([40,30,20,10,0],3).stdout.strip(),'1 20 0')
        self.assertEqual(self.read([None]*5).stdout.strip(),'0 0 0')

    def test_rejects_bad_pressure_or_schema_on_all_ranks(self):
        for p in ([40,30,20,10,-1],[40,30,35,10,0],[40,30,20,10,1],
                  [40,30,float('nan'),10,0],[40,30,None,10,0]):
            self.assertNotEqual(self.read(p,3).returncode,0,p)
        self.assertNotEqual(self.read([40,30,20,10,0],extra=' junk').returncode,0)
        self.assertNotEqual(self.read([40,30,20,10,0],extra=' 123').returncode,0)

if __name__=='__main__': unittest.main()
