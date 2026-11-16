"""Exercise pPv in the production transport, nodal output and viscosity paths.

The unchanged background conductivity factors are fixed at one. The pPv
helpers, Gauss-point transport, nodal diagnostics and viscosity update are
compiled from production C and compared with independent analytic values.
"""
import ctypes
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]


def function(source, signature):
    start = source.index(signature)
    opening = source.index("{", start)
    depth, end = 1, opening + 1
    while depth:
        depth += (source[end] == "{") - (source[end] == "}")
        end += 1
    return source[start:end]


class PPvTransport(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        compiler = shutil.which(os.environ.get("MPICC", "mpicc"))
        if compiler is None:
            raise RuntimeError("mpicc is required to test production C")
        phase_source = (ROOT / "lib/Phase_change.c").read_text()
        phase = "\n".join(function(phase_source, signature) for signature in (
            "void phase_change_state(", "double pPv_fraction_at_state(",
            "double pPv_fraction_at_node(", "double pPv_fraction_at_gp(",
            "double pPv_conductivity_factor(", "double pPv_viscosity_factor("))
        advection = (ROOT / "lib/Advection_diffusion.c").read_text()
        # The source has a forward declaration with the same signature.
        transport = function(advection[advection.index(
            "void thermal_transport_at_gp(", advection.index(
                "void thermal_transport_at_gp(") + 1):],
            "void thermal_transport_at_gp(")
        material = function((ROOT / "lib/Material_properties.c").read_text(),
                            "double nodal_conductivity_diagnostic(")
        viscosity = function((ROOT / "lib/Viscosity_structures.c").read_text(),
                             "static void visc_from_pPv(")
        source = r'''
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "phase_change.h"
#include "material_properties.h"
void myerror(struct All_variables *E,char *message) {abort();}
void parallel_process_termination(void) {abort();}
double conductivity_depth_factor(struct All_variables *E,double depth) {return 1.;}
double conductivity_temperature_factor(struct All_variables *E,double t) {return 1.;}
double conductivity_composition_factor(struct All_variables *E,double c) {return 1.;}
double conductivity_element_prefactor(struct All_variables *E,int c,int e,double ref) {return ref;}
''' + phase + transport + material + viscosity + r'''
void evaluate(int enabled,double k_ratio,double eta_ratio,double temperature,
              int poison_geometry,double *out) {
    static struct All_variables state;
    struct All_variables *E=&state;
    struct IEN ien[2];
    double rho_ref[3]={0,1.5,1.5},cp_ref[3]={0,1.25,1.25},gravity[3]={0,1,1};
    double coordinates[4][9]={{0}},radius[9]={0},T[9]={0};
    double tgp[9]={0},rho[9]={0},cp[9]={0},k[9]={0},kap[9]={0};
    float eta[9]={0},*evisc[NCS]={0},stale_fraction[9]={0};
    int a,d,gp;
    memset(E,0,sizeof(*E));
    E->mesh.nsd=3;E->lmesh.noz=2;E->lmesh.nel=1;E->sphere.caps_per_proc=1;
    E->sphere.ro=1.;E->control.reference_conductivity=2.;
    E->control.pPv_enabled=enabled;E->control.pPv_k_ratio=k_ratio;
    E->control.pPv_eta_ratio=eta_ratio;
    E->control.phase[pPv_PHASE_INDEX].depth=.5;
    E->control.phase[pPv_PHASE_INDEX].transT=.5;
    E->control.phase[pPv_PHASE_INDEX].clapeyron=.25;
    E->control.phase[pPv_PHASE_INDEX].inv_width=4.;
    E->refstate.rho=rho_ref;E->refstate.gravity=gravity;E->refstate.heat_capacity=cp_ref;
    E->ien[1]=ien;E->T[1]=T;E->sx[1][3]=radius;evisc[1]=eta;
    E->phase_B[pPv_PHASE_INDEX][1]=stale_fraction;
    for(d=1;d<=3;d++)E->x[1][d]=poison_geometry?NULL:coordinates[d];
    for(a=1;a<=8;a++) {
        ien[1].node[a]=a;T[a]=temperature;eta[a]=1000.;stale_fraction[a]=1.;
        coordinates[1][a]=.3;coordinates[2][a]=(a%2?-.1:.1);coordinates[3][a]=.4;
        radius[a]=sqrt(.26);
        for(gp=1;gp<=8;gp++)E->N.vpt[GNVINDEX(a,gp)]=.125;
        tgp[a]=temperature;
    }
    thermal_transport_at_gp(E,1,1,tgp,2.,rho,cp,k,kap);
    out[0]=k[1];out[1]=kap[1];
    out[2]=nodal_conductivity_diagnostic(E,1,1,CONDUCTIVITY_K_TOTAL);
    out[3]=nodal_conductivity_diagnostic(E,1,1,CONDUCTIVITY_KAPPA_EFF);
    visc_from_pPv(E,evisc);out[4]=eta[1];out[5]=rho[1];out[6]=cp[1];
}
'''
        cls.directory = tempfile.TemporaryDirectory(prefix="ppv-transport-")
        cls.addClassCleanup(cls.directory.cleanup)
        path = Path(cls.directory.name)
        (path / "transport.c").write_text(source)
        subprocess.run([compiler, "-std=gnu99", "-shared", "-fPIC",
                        "-I" + str(ROOT / "lib"), str(path / "transport.c"),
                        "-lm", "-o", str(path / "transport.so")],
                       check=True, capture_output=True, text=True)
        cls.library = ctypes.CDLL(str(path / "transport.so"))
        cls.library.evaluate.argtypes = [ctypes.c_int] + [ctypes.c_double] * 3 + [
            ctypes.c_int, ctypes.POINTER(ctypes.c_double)]

    def evaluate(self, enabled, k_ratio, eta_ratio, temperature, poison=False):
        out = (ctypes.c_double * 7)()
        self.library.evaluate(enabled, k_ratio, eta_ratio, temperature, poison, out)
        return list(out)

    def test_current_gauss_and_node_fractions_control_transport(self):
        for temperature in (.1, .5, .9):
            out = self.evaluate(1, 2., .25, temperature)
            fraction_gp = .5 * (1 + math.tanh(.5 - temperature))
            fraction_node = .5 * (1 + math.tanh(4 * (
                (.5 - math.sqrt(.26)) * 1.5 - .25 * (temperature - .5))))
            self.assertAlmostEqual(out[0], 2 * (1 + fraction_gp), places=13)
            self.assertAlmostEqual(out[1], out[0] / 1.875, places=13)
            self.assertAlmostEqual(out[2], 2 * (1 + fraction_node), places=13)
            self.assertAlmostEqual(out[3], out[2] / 1.875, places=13)
            self.assertAlmostEqual(out[4], 1000 * .25**fraction_gp, delta=6e-5)
            self.assertEqual(out[5:], [1.5, 1.25])

    def test_disabled_or_neutral_ratios_do_not_need_phase_geometry(self):
        for enabled, k_ratio, eta_ratio in [(0, 2., .25), (1, 1., 1.)]:
            out = self.evaluate(enabled, k_ratio, eta_ratio, .5, poison=True)
            self.assertEqual(out, [2., 2 / 1.875, 2., 2 / 1.875, 1000., 1.5, 1.25])


if __name__ == "__main__":
    unittest.main()
