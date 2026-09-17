"""Audit actual eight-GP spherical element geometry, energy and phase force.

Uses current strict cfg/reference data. External conductivity is constant and
flux BCs are disabled, while the production transport interpolation, complete
element_residual, phase nodal routines and geometry functions run unchanged.
"""
from __future__ import annotations

import configparser
import ctypes
import os
import re
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

import numpy as np
from test_phase_energy_geometry import function

ROOT=Path(__file__).resolve().parents[2]
RUNS=ROOT.parents[1]/"runs"


class PhaseElementTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        compiler=shutil.which(os.environ.get("MPICC","mpicc"))
        if not compiler:
            raise RuntimeError("mpicc required for actual C element audit")
        def extract(file, signature):
            source=(ROOT/"lib"/file).read_text()
            starts=[m.start() for m in re.finditer(r"^"+re.escape(signature),source,re.M)
                    if not source[source.index(")",m.start())+1:].lstrip().startswith(";")]
            return function(source[starts[-1]:],signature)
        prefix='''
#include <math.h>
#include <string.h>
#include "element_definitions.h"
#include "global_defs.h"
double myatan(double y,double x) {return atan2(y,x);}
void parallel_process_termination(void) {abort();}
void get_global_1d_shape_fn() {abort();} /* never used: FLAGS=NULL */
double conductivity_element_prefactor(struct All_variables *E,int c,int e,double k) {return k;}
double conductivity_temperature_factor(struct All_variables *E,double T) {return 1.;}
static void apply_one_phase(struct All_variables *,double **,int);
static void calc_phase_change(struct All_variables *,int);
static void debug_phase_change(struct All_variables *,int);
'''
        selections={
            "Shape_functions.c": ["double lpoly(","double lpolydash(","void construct_shape_functions("],
            "General_matrix_functions.c": ["double determinant(","double cofactor("],
            "Size_does_matter.c": ["void form_rtf_bc(","void get_global_shape_fn("],
            "Phase_change.c": ["void phase_change_state(","float phase_change_reference_fraction(",
                "static void calc_phase_change(","static void debug_phase_change(",
                "static void apply_one_phase(","void phase_change_apply("],
            "Advection_diffusion.c": ["static void pg_shape_fn(","static void element_thermal_transport(","static void element_residual("],
        }
        matrix=(ROOT/"lib/General_matrix_functions.c").read_text()
        epsilon=matrix[matrix.index("int epsilon["):matrix.index(";",matrix.index("int epsilon["))+1]
        code=prefix+epsilon+"\n"+"\n".join(extract(f,s) for f,ss in selections.items() for s in ss)
        code+=(Path(__file__).with_name("phase_element_harness.c")).read_text()
        cls.temp=tempfile.TemporaryDirectory(prefix="phase-element-")
        cls.addClassCleanup(cls.temp.cleanup)
        folder=Path(cls.temp.name); (folder/"audit.c").write_text(code)
        result=subprocess.run([compiler,"-std=gnu89","-shared","-fPIC",
            "-Wno-deprecated-non-prototype","-I",str(ROOT/"lib"),str(folder/"audit.c"),
            "-lm","-o",str(folder/"audit.so")],capture_output=True,text=True)
        if result.returncode:
            raise RuntimeError(result.stderr)
        cls.dll=ctypes.CDLL(str(folder/"audit.so"))
        cls.call=cls.dll.audit_element
        pointer=np.ctypeslib.ndpointer(dtype=np.float64,flags="C_CONTIGUOUS")
        cls.call.argtypes=[ctypes.c_double]*4+[pointer,pointer]+[ctypes.c_double]*4+[ctypes.c_int,pointer]
        cls.call.restype=None
        cfg=configparser.ConfigParser(interpolation=None); cfg.read(RUNS/"cmbhf_ALA_strict.cfg")
        def vector(key): return np.array([float(v) for v in cfg["CitcomS.solver.phase"][key].split(',')])
        rows=np.stack([vector(k) for k in ("phase_depth","phase_clapeyron","phase_width",
            "phase_transT","phase_delta_s","phase_delta_rho")],axis=1)
        const=cfg["CitcomS.solver.const"]
        rayleigh=cfg.getfloat("CitcomS.solver","rayleigh")
        ra=rayleigh*rows[:,5]/(const.getfloat("rho0")*const.getfloat("alpha0")*const.getfloat("deltaT"))
        cls.parameters=np.r_[np.column_stack([rows,ra]).ravel(),cfg.getfloat("CitcomS.solver","surfaceT")]
        cls.radii=np.loadtxt(RUNS/"GLB.coor.global.dat",skiprows=1,usecols=1)
        cls.ref=np.loadtxt(RUNS/"refstate_ALA_strict.txt")

    def evaluate(self, phase=0, shift=.0, rate=.03, velocity=.2, conductivity=0.,mask=7):
        radius=1.-self.parameters[7*phase]
        index=np.searchsorted(self.radii,radius)-1
        rin,rout=self.radii[index:index+2]
        reference=np.ascontiguousarray(self.ref[index:index+2][:,[0,1,4,2]].ravel())
        out=np.zeros(129)
        self.call(rin,rout,.012,.014,reference,self.parameters,shift,rate,velocity,conductivity,mask,out)
        return out,reference

    def expected_gp(self, out, rate, velocity, mask):
        gp=out[:64].reshape(8,8)
        total=np.zeros(8)
        for p in range(3):
            if not (mask & (1<<p)): continue
            z0,gamma,w,t0,entropy,_,_=self.parameters[7*p:7*p+7].astype(np.float32).astype(float)
            # C first reads width into float, then stores its reciprocal as float.
            invw=float(np.float32(1./w))
            r,T,rho,rg,drg,gradT=(gp[:,k] for k in (1,2,3,4,5,6))
            eps=1e-7
            actual_velocity=float(np.float32(velocity))
            dT=rate+actual_velocity*gradT
            def fraction(sign):
                q=(1.-(r+sign*eps*actual_velocity)-z0)*(rg+sign*eps*actual_velocity*drg)-gamma*(T+sign*eps*dT-t0)
                return .5*(1.+np.tanh(invw*q))
            total+=rho*(T+float(np.float32(self.parameters[21])))*entropy*(fraction(1)-fraction(-1))/(2*eps)
        return total

    def test_complete_geometry_and_residual_phase_contribution(self):
        for phase in range(3):
            for velocity in (-.2,.2):
                for conductivity in (0.,1.):
                    with self.subTest(phase=phase,v=velocity,k=conductivity):
                        out,_=self.evaluate(phase,shift=.01,velocity=velocity,conductivity=conductivity)
                        gp=out[:64].reshape(8,8)
                        np.testing.assert_allclose(gp[:,0]*gp[:,1],1.,rtol=3e-15)
                        self.assertTrue(np.all(gp[:,7]>0.))
                        q=self.expected_gp(out,.03,velocity,7)
                        # Partition of unity: sum of nodal phase residuals = -integral(Qphase).
                        self.assertAlmostEqual(out[64:72].sum(),-np.dot(q,gp[:,7]),delta=1e-7*abs(np.dot(q,gp[:,7])))
                        self.assertAlmostEqual(out[80],q.mean(),delta=1e-7*abs(q.mean()))

    def test_reference_subtraction_and_actual_phase_buoyancy(self):
        for phase in range(3):
            for shift in (0.,-.05,.05):
                out,reference=self.evaluate(phase,shift=shift)
                expected=np.zeros(8)
                for p in range(3):
                    x=out[81+16*p:89+16*p]; xref=out[89+16*p:97+16*p]
                    self.assertTrue(np.all((x>=0)&(x<=1)))
                    z0,gamma,width,t0=self.parameters[7*p:7*p+4].astype(np.float32).astype(float)
                    invwidth=float(np.float32(1./width))
                    ri=np.searchsorted(self.radii,1.-self.parameters[7*phase])-1
                    radii=np.tile(self.radii[ri:ri+2],4)
                    ref=np.tile(reference.reshape(2,4),(4,1))
                    for actual,delta in ((x,shift),(xref,0.)):
                        q=(1.-radii-z0)*ref[:,0]*ref[:,1]-gamma*(ref[:,3]+delta-t0)
                        target=(.5*(1.+np.tanh(invwidth*q))).astype(np.float32).astype(float)
                        np.testing.assert_allclose(actual,target,rtol=1e-7,atol=1e-8)
                    if shift:
                        # Positive Gamma: hotter material has less dense phase;
                        # negative Gamma reverses the temperature response.
                        self.assertTrue(np.all((x-xref)*gamma*shift<=0.))
                    if shift==0.: np.testing.assert_array_equal(x,xref)
                    ra=float(np.float32(self.parameters[7*p+6]))
                    # Production B-Xref is a float subtraction.
                    expected-=ra*(x.astype(np.float32)-xref.astype(np.float32)).astype(float)
                np.testing.assert_allclose(out[72:80],expected,rtol=3e-7,atol=1e-8)
                if shift==0.: np.testing.assert_array_equal(out[72:80],np.zeros(8))

    def test_three_latent_terms_are_added_once(self):
        combined,_=self.evaluate(mask=7)
        pieces=[self.evaluate(mask=1<<p)[0] for p in range(3)]
        np.testing.assert_allclose(combined[64:72],sum(x[64:72] for x in pieces),rtol=1e-12,atol=1e-16)
        self.assertAlmostEqual(combined[80],sum(x[80] for x in pieces),places=12)


if __name__=="__main__": unittest.main()
