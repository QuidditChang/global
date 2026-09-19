"""Numerically exercise production CBF geometry and shared physical RHS.

The volume fixture deliberately supplies angular GNx and a nonunit radius;
using them as Cartesian derivatives would fail the independent expectation.
This is an operator test, not a solved spherical-shell convergence benchmark.
"""
import ctypes as C
import math
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]

def function(source, signature):
    start = source.index(signature)
    opening = source.index('{', start)
    end, depth = opening + 1, 1
    while depth:
        depth += (source[end] == '{') - (source[end] == '}')
        end += 1
    return source[start:end]

class CBFKernel(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        energy = (ROOT/'lib/Advection_diffusion.c').read_text()
        # Skip forward declarations; extract actual production definitions.
        energy = energy[energy.index('static void element_thermal_transport(', energy.index('static void pg_solver(', energy.index('static void pg_solver(')+1)):]
        parts = [function(energy, 'static void element_thermal_transport('),
                 function(energy, 'static void element_residual('),
                 function(energy, 'void cbf_element_thermal_residual(')]
        phase = function((ROOT/'lib/Phase_change.c').read_text(), 'void phase_change_state(')
        prefix = r'''
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "cbf_geometry.h"
static double test_velocity;
double conductivity_element_prefactor(struct All_variables *E,int m,int el,double k) {return k;}
double conductivity_temperature_factor(struct All_variables *E,double t) {return 1.0;}
void parallel_process_termination(void) {abort();}
void get_global_1d_shape_fn() {abort();}
void velo_from_element(struct All_variables *E,float v[4][9],int m,int el,int sphere) {
 int a; memset(v,0,4*9*sizeof(float)); for(a=1;a<=8;a++) v[1][a]=test_velocity;
}
void get_global_shape_fn(struct All_variables *E,int el,struct Shape_function *GN,
 struct Shape_function_dx *g,struct Shape_function_dA *w,int p,int sphere,
 double rtf[4][9],int lev,int m) {
 int a,i; memset(g,0,sizeof(*g));
 for(i=1;i<=8;i++) {
  w->vpt[i]=0.125; rtf[1][i]=1.1; rtf[3][i]=2.0;
  for(a=1;a<=8;a++) g->vpt[GNVXINDEX(0,a,i)]=(a%2 ? -1.0:1.0)/8.0;
 }
}
'''
        wrapper = r'''
int ray_value(const double *xyz,const double *q,const double *ray,double *out) {return cbf_face_ray_value((const double (*)[3])xyz,q,ray,out);}
void mass(const double *xyz,double *out) {cbf_face_gll_mass((const double (*)[3])xyz,out);}
double jac(const double *xyz,double u,double v) {return cbf_face_jacobian((const double (*)[3])xyz,u,v);}
void residual(double rate,double velocity,double k,double source,double entropy,double *out) {
 struct All_variables state,*E=&state; struct IEN ien[2];
 double T[10],Tdot[10],rho[3]={0,1,1},cp[3]={0,1,1},gravity[3]={0,1,1};
 double phase[2]={0,123.456},adi[2]={0,0},visc[2]={0,0};
 unsigned int flags[10]={0}; int i,a;
 memset(E,0,sizeof(*E)); E->mesh.nsd=3; E->mesh.dof=3;
 E->lmesh.noz=2; E->lmesh.elz=1; E->lmesh.nno=8;
 E->ien[1]=ien; E->T[1]=T; E->Tdot[1]=Tdot; E->node[1]=flags;
 E->refstate.rho=rho; E->refstate.heat_capacity=cp; E->refstate.gravity=gravity;
 E->heating_phase[1]=phase; E->heating_adi[1]=adi; E->heating_visc[1]=visc;
 E->control.reference_conductivity=k; E->control.Q0=source; E->data.scalev=1000;
 E->sphere.ro=1; E->control.surface_temp=0.1;
 E->control.phase[0].entropy_jump=entropy;
 E->control.phase[0].depth=0.5; E->control.phase[0].transT=0.5;
 E->control.phase[0].clapeyron=0.3; E->control.phase[0].inv_width=2;
 for(a=1;a<=8;a++) {
  ien[1].node[a]=a; T[a]=a%2 ? 0:1; Tdot[a]=rate;
  for(i=1;i<=8;i++) E->N.vpt[GNVINDEX(a,i)]=0.125;
 }
 test_velocity=velocity;
 cbf_element_thermal_residual(E,1,1,out);
 out[0]=phase[1];
}
'''
        src=Path(cls.tmp.name)/'fixture.c';src.write_text(prefix+phase+'\n'+'\n'.join(parts)+wrapper)
        lib=Path(cls.tmp.name)/'fixture.so'
        subprocess.run([shutil.which('mpicc'),'-std=gnu99','-w','-shared','-fPIC','-I'+str(ROOT/'lib'),str(src),'-lm','-o',str(lib)],check=True)
        cls.lib=C.CDLL(str(lib)); ptr=C.POINTER(C.c_double)
        cls.lib.ray_value.argtypes=[ptr,ptr,ptr,ptr];
        cls.lib.mass.argtypes=[ptr,ptr]; cls.lib.jac.argtypes=[ptr,C.c_double,C.c_double];cls.lib.jac.restype=C.c_double
        cls.lib.residual.argtypes=[C.c_double]*5+[ptr]

    @classmethod
    def tearDownClass(cls): cls.tmp.cleanup()

    def masses(self, points):
        x=(C.c_double*12)(*sum(points,[])); y=(C.c_double*4)();self.lib.mass(x,y);return x,list(y)

    def test_affine_mass_rotation_and_orientation(self):
        points=[[0,0,0],[2,0,0],[2,3,0],[0,3,0]]
        for pts in [points,list(reversed(points)),[[p[0],0,p[1]] for p in points]]:
            _,m=self.masses(pts)
            for a in m:self.assertAlmostEqual(a,1.5,places=14)

    def test_warped_face_uses_endpoint_metric(self):
        x,m=self.masses([[0,0,0],[2,0,0],[2,3,1],[0,3,0]])
        expected=[1.5,math.sqrt(10)/2,7/4,math.sqrt(5)*3/4]
        for a,b in zip(m,expected):self.assertAlmostEqual(a,b,places=14)
        g=1/math.sqrt(3)
        gauss_area=sum(self.lib.jac(x,u,v) for u in [-g,g] for v in [-g,g])
        self.assertGreater(abs(sum(m)-gauss_area),0.01)

    def test_shared_rhs_signs_metric_and_no_state_mutation(self):
        out=(C.c_double*9)()
        # grad_theta T = .5 and r^-1 = 2: physical gradient is 1.
        for rate,velocity,k,source in [(0,0,2,0),(3,0,0,0),(0,4,0,0),(0,0,0,5),(3,4,2,5)]:
            self.lib.residual(rate,velocity,k,source,0,out)
            self.assertEqual(out[0],123.456)
            for a in range(1,9):
                physical_grad_N=(-1 if a%2 else 1)/4
                expected=(source-rate-velocity)/8-k*physical_grad_N
                self.assertAlmostEqual(out[a],expected,places=13)

    def test_ray_mapping_pole_seam_and_outside(self):
        x,_=self.masses([[-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]])
        q=(C.c_double*4)(0,2,4,2);out=C.c_double()
        for direction,expected in [([0,0,1],2),([1,1,1],4),([-1,-1,1],0),([0.5,0,1],2.5)]:
            norm=math.sqrt(sum(t*t for t in direction));r=(C.c_double*3)(*(t/norm for t in direction))
            self.assertEqual(self.lib.ray_value(x,q,r,C.byref(out)),1)
            self.assertAlmostEqual(out.value,expected,places=12)
        self.assertEqual(self.lib.ray_value(x,q,(C.c_double*3)(0,0,-1),C.byref(out)),0)
        # A face crossing lon=0; +0 and -0 longitudes must agree.
        x,_=self.masses([[1,-1,-1],[1,1,-1],[1,1,1],[1,-1,1]])
        for phi in [0,2*math.pi]:
            r=(C.c_double*3)(math.cos(phi),math.sin(phi),0)
            self.assertEqual(self.lib.ray_value(x,q,r,C.byref(out)),1)
            self.assertAlmostEqual(out.value,2,places=12)

    def test_nonzero_phase_energy_is_included(self):
        out=(C.c_double*9)();self.lib.residual(3,0,0,0,0.2,out)
        # q=0: dX/dT=-0.5*clapeyron*inv_width=-0.3.
        phase=(0.5+C.c_float(0.1).value)*C.c_float(0.2).value*(-0.5*C.c_float(0.3).value*2)*3
        for a in range(1,9):self.assertAlmostEqual(out[a],-(3+phase)/8,places=12)
        self.assertEqual(out[0],123.456)

if __name__=='__main__':unittest.main()
