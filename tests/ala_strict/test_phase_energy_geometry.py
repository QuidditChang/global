"""Execute the production GP phase-energy block, including its geometry call.

Compile the actual form_rtf_bc, phase_change_state and element_residual phase
loop. Compare material phase derivatives to independent centered differences;
checking phase_change_state alone would not catch the former call-site error.
"""
from __future__ import annotations

import ctypes
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]


def function(source: str, signature: str) -> str:
    start = source.index(signature)
    opening = source.index("{", start)
    depth = 1
    end = opening + 1
    while depth:
        depth += (source[end] == "{") - (source[end] == "}")
        end += 1
    return source[start:end]


class PhaseEnergyGeometryTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        compiler = shutil.which(os.environ.get("MPICC", "mpicc"))
        if compiler is None:
            raise RuntimeError("mpicc is required to test the production C block")
        geometry = function((ROOT / "lib/Size_does_matter.c").read_text(),
                            "void form_rtf_bc(")
        phase = function((ROOT / "lib/Phase_change.c").read_text(),
                         "void phase_change_state(")
        energy = (ROOT / "lib/Advection_diffusion.c").read_text()
        start = energy.index("    for(i=1;i<=vpts;i++) {\n      double material_dT")
        block = energy[start:energy.index("    E->heating_phase[m][el] = 0.0;", start)]
        prefix = r'''
#include <math.h>
#include <string.h>
#include "element_definitions.h"
#include "global_defs.h"
double myatan(double y, double x) { return atan2(y,x); }
'''
        wrapper = r'''
double energy_at(double ro, double radius, double phase_depth,
                 double temperature, double transT, double clapeyron,
                 double inv_width, double density, double density_gradient,
                 double temperature_rate, double temperature_gradient,
                 double radial_velocity, double entropy) {
    static struct All_variables state;
    struct All_variables *E=&state;
    struct Shape_function_dx GNx;
    struct IEN connectivity[2];
    double ref_rho[9], ref_g[9];
    double rtf[4][9]={{0}}, bc[4][4]={{0}}, x[4];
    double dT[9]={0}, v1[9]={0}, v2[9]={0}, v3[9]={0};
    double tx1[9]={0}, tx2[9]={0}, tx3[9]={0};
    double tgp[9]={0}, rho_gp[9]={0}, phase_energy[9]={0};
    double sfn, h=0.001;
    int i,j,node,nz,ends=8,vpts=1,m=1,el=1;
    memset(E,0,sizeof(*E)); memset(&GNx,0,sizeof(GNx));
    E->sphere.ro=ro; E->lmesh.noz=8; E->ien[m]=connectivity;
    E->refstate.rho=ref_rho; E->refstate.gravity=ref_g;
    E->control.surface_temp=0.1;
    E->control.phase[0].depth=phase_depth;
    E->control.phase[0].transT=transT;
    E->control.phase[0].clapeyron=clapeyron;
    E->control.phase[0].inv_width=inv_width;
    E->control.phase[0].entropy_jump=entropy;
    for(j=1;j<=ends;j++) {
        double sign=j<=4 ? -1.0 : 1.0;
        connectivity[el].node[j]=j;
        ref_rho[j]=density+sign*h*density_gradient; ref_g[j]=1.0;
        E->N.vpt[GNVINDEX(j,1)]=0.125;
        GNx.vpt[GNVXINDEX(2,j,1)]=sign/(8.0*h);
    }
    x[1]=radius*sin(1.1)*cos(0.4);
    x[2]=radius*sin(1.1)*sin(0.4); x[3]=radius*cos(1.1);
    form_rtf_bc(1,x,rtf,bc);
    tgp[1]=temperature; rho_gp[1]=density;
    dT[1]=temperature_rate; tx3[1]=temperature_gradient;
    v3[1]=radial_velocity;
'''
        source = prefix + geometry + "\n" + phase + "\n" + wrapper + block + "return phase_energy[1];\n}\n"
        cls.directory = tempfile.TemporaryDirectory(prefix="phase-geometry-")
        cls.addClassCleanup(cls.directory.cleanup)
        path = Path(cls.directory.name)
        (path / "test.c").write_text(source)
        subprocess.run([compiler, "-std=gnu89", "-shared", "-fPIC",
                        "-Wno-deprecated-non-prototype", "-I", str(ROOT / "lib"),
                        str(path / "test.c"), "-lm", "-o", str(path / "test.so")],
                       check=True, capture_output=True, text=True)
        cls.library = ctypes.CDLL(str(path / "test.so"))
        cls.energy = cls.library.energy_at
        cls.energy.argtypes = [ctypes.c_double] * 13
        cls.energy.restype = ctypes.c_double

    @staticmethod
    def fraction(ro, radius, z0, temperature, t0, gamma, width, rho):
        import math
        return 0.5 * (1.0 + math.tanh(((ro-radius-z0)*rho-gamma*(temperature-t0))/width))

    def test_phase_centers_and_material_derivative(self):
        # Use both unit and non-unit outer radii, either velocity direction,
        # nonzero d(rho*g)/dr and nonzero temperature advection/time derivative.
        for ro in (1.0, 0.97):
            for z0, gamma, entropy in ((410/6371, .02, -.03),
                                       (520/6371, .03, -.02),
                                       (660/6371, -.01, .04)):
                # The production phase parameters are floats, even though
                # geometry and the phase derivative arithmetic use doubles.
                z0, gamma, entropy = (ctypes.c_float(x).value
                                      for x in (z0, gamma, entropy))
                for offset in (-.002, 0., .002):
                    for velocity in (-.2, .2):
                        with self.subTest(ro=ro, z0=z0, offset=offset, v=velocity):
                            radius=ro-z0+offset
                            rho, drho, T, T0, width = 1.2, -.7, .5, .5, .004
                            dT, gradT = .03, -.15
                            actual = self.energy(ro,radius,z0,T,T0,gamma,1/width,
                                                 rho,drho,dT,gradT,velocity,entropy)
                            eps=1e-7
                            plus=self.fraction(ro,radius+eps*velocity,z0,
                                T+eps*(dT+velocity*gradT),T0,gamma,width,
                                rho+eps*velocity*drho)
                            minus=self.fraction(ro,radius-eps*velocity,z0,
                                T-eps*(dT+velocity*gradT),T0,gamma,width,
                                rho-eps*velocity*drho)
                            expected=rho*(T+ctypes.c_float(.1).value)*entropy*(plus-minus)/(2*eps)
                            self.assertAlmostEqual(actual,expected,delta=2e-7*abs(expected))

    def test_zero_entropy_and_static_material_state(self):
        args=[1.,1.-410/6371,410/6371,.5,.5,.02,250.,1.2,-.7,.03,-.15,.2,0.]
        self.assertEqual(self.energy(*args),0.)
        args[-1]=-.03; args[9]=0.; args[11]=0.
        self.assertEqual(self.energy(*args),0.)

    def test_former_inverse_radius_call_suppresses_center_transition(self):
        z0=410/6371; radius=1.-z0
        correct=self.fraction(1.,radius,z0,.5,.5,.02,.004,1.2)
        former=self.fraction(1.,1./radius,z0,.5,.5,.02,.004,1.2)
        self.assertAlmostEqual(correct,.5)
        self.assertLess(former,1e-12)
        self.assertGreater(abs(self.energy(1.,radius,z0,.5,.5,.02,250.,
                                           1.2,0.,0.,0.,.2,-.03)),.1)


if __name__ == "__main__":
    unittest.main()
