"""Run the actual C parser, phase validation and buoyancy kernel for pPv."""
import ctypes
import math
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]


class PhaseKernelTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory(prefix="pPv-phase-")
        cls.addClassCleanup(cls.tmp.cleanup)
        cls.path = Path(cls.tmp.name)
        source = r'''
#include <setjmp.h>
#include <string.h>
#include "Phase_change.c"
#include "Parsing.c"
static jmp_buf failure;
void parallel_process_termination(void) { longjmp(failure, 1); }
void unique_copy_file(struct All_variables *E, char *a, char *b) {}
int read_case(char *path, double *out) {
    static struct All_variables state;
    struct All_variables *E=&state;
    int i;
    memset(E,0,sizeof(*E));
    E->sphere.ro=1; E->sphere.ri=0.546225;
    E->data.rho0=4000; E->data.alpha0=3e-5;
    E->data.ref_temperature=3400; E->data.radius_km=6371;
    E->data.Ttop=300; E->data.Cp0=1280; E->data.g0=9.81;
    E->control.Atemp=1e8; E->parallel.me=1;
    if(setjmp(failure)) return 1;
    /* Production initializes the parser once per process; isolate test calls. */
    memset(&ext_par,0,sizeof(ext_par));
    setup_parser(E,path); phase_change_input(E); shutdown_parser(E);
    validate_phase_parameters(E);
    out[0]=phase_change_count(E); out[1]=E->control.pPv_enabled;
    out[2]=E->control.pPv_latent_enabled;
    out[3]=E->control.pPv_k_ratio; out[4]=E->control.pPv_eta_ratio;
    for(i=0;i<PHASE_TRANSITIONS;i++) {
      out[5+4*i]=E->control.phase[i].depth;
      out[6+4*i]=E->control.phase[i].density_jump;
      out[7+4*i]=E->control.phase[i].Ra;
      out[8+4*i]=phase_change_entropy(E,i);
    }
    return 0;
}
void response(double T, int enabled, int latent, double *out) {
    static struct All_variables state;
    struct All_variables *E=&state;
    double rho[3]={0,1,1},g[3]={0,1,1},tref[3]={0,.6,.6};
    double radius[3]={0,.6,.7},buoyancy[3]={0},*buoy[2]={0,buoyancy};
    double temp[3]={0,T,T}; float b[4][3]={{0}},boundary[4][2]={{0}};
    int i;
    memset(E,0,sizeof(*E));
    E->sphere.ro=1; E->sphere.caps_per_proc=1;
    E->lmesh.nno=E->lmesh.noz=2; E->lmesh.nox=E->lmesh.noy=1;
    E->sx[1][3]=radius; E->T[1]=temp;
    E->refstate.rho=rho; E->refstate.gravity=g; E->refstate.Tref=tref;
    E->control.pPv_enabled=enabled; E->control.pPv_latent_enabled=latent;
    E->control.pPv_k_ratio=2; E->control.pPv_eta_ratio=.01;
    for(i=0;i<4;i++) {
      E->phase_B[i][1]=b[i]; E->phase_boundary[i][1]=boundary[i];
      E->control.phase[i].depth=.4; E->control.phase[i].transT=.6;
      E->control.phase[i].clapeyron=.1; E->control.phase[i].inv_width=100;
      E->control.phase[i].Ra=i==0 ? 2 : (i==3 ? 3 : 0);
      E->control.phase[i].entropy_jump=-.01;
    }
    phase_change_apply(E,buoy);
    out[0]=pPv_fraction_at_node(E,1,1,T);
    out[1]=phase_change_density_fraction(E,0,1,1,b[0][1]);
    out[2]=phase_change_density_fraction(E,3,1,1,b[3][1]);
    out[3]=buoyancy[1]; out[4]=phase_change_entropy(E,3);
    out[5]=pPv_conductivity_factor(E,out[0]);
    out[6]=pPv_viscosity_factor(E,out[0]);
}
'''
        (cls.path / "fixture.c").write_text(source)
        subprocess.run([shutil.which("mpicc"), "-std=gnu89", "-shared", "-fPIC",
                        "-Wno-implicit-function-declaration", "-Wno-deprecated-non-prototype",
                        "-I" + str(ROOT / "lib"), str(cls.path / "fixture.c"),
                        "-lm", "-o", str(cls.path / "fixture.so")],
                       capture_output=True, text=True, check=True)
        cls.lib = ctypes.CDLL(str(cls.path / "fixture.so"))
        cls.lib.read_case.argtypes = [ctypes.c_char_p, ctypes.POINTER(ctypes.c_double)]
        cls.lib.response.argtypes = [ctypes.c_double, ctypes.c_int, ctypes.c_int,
                                    ctypes.POINTER(ctypes.c_double)]

    def read_case(self, extra=""):
        cfg = self.path / "case.cfg"
        cfg.write_text("phase_depth=.06,.08,.10\nphase_delta_rho=180,58,459\n"
                       "phase_delta_s=-.03,-.02,.01\nphase_width=.002,.002,.002\n"
                       "phase_clapeyron=.05,.08,-.02\nphase_transT=.4,.5,.6\n" + extra)
        out = (ctypes.c_double * 21)()
        status = self.lib.read_case(str(cfg).encode(), out)
        return status, list(out)

    PPV = ("pPv_enabled=on\npPv_depth=.43164\npPv_transT=.64706\n"
           "pPv_delta_rho=66.4\npPv_width=.00306\npPv_clapeyron=.1088\n"
           "pPv_delta_s=-.014\n")

    def test_old_three_entry_input_defaults_off(self):
        status, out = self.read_case()
        self.assertEqual(status, 0)
        self.assertEqual(out[:5], [3, 0, 0, 1, 1])
        self.assertEqual(out[6:16:4], [180, 58, 459])
        self.assertEqual(out[17:], [0, 0, 0, 0])

    def test_scalar_input_and_density_derived_rayleigh(self):
        status, out = self.read_case(self.PPV)
        self.assertEqual(status, 0)
        self.assertEqual(out[:5], [4, 1, 0, 1, 1])
        self.assertAlmostEqual(out[18], 66.4, places=4)
        self.assertAlmostEqual(out[19] / (1e8*66.4/(4000*3e-5*3400)), 1, places=6)
        self.assertEqual(out[20], 0)
        status, out = self.read_case(self.PPV + "pPv_latent_enabled=on\n")
        self.assertEqual(status, 0)
        self.assertAlmostEqual(out[20], -.014, places=8)

    def test_invalid_enabled_inputs_rejected(self):
        for entry in ("pPv_depth=.5", "pPv_width=0", "pPv_delta_rho=-1",
                      "pPv_k_ratio=0", "pPv_eta_ratio=-1", "pPv_transT=-1",
                      "pPv_latent_enabled=on\npPv_delta_s=.01"):
            with self.subTest(entry=entry):
                self.assertEqual(self.read_case(self.PPV + entry + "\n")[0], 1)

    def response(self, T, enabled=1, latent=0):
        out = (ctypes.c_double * 7)()
        self.lib.response(T, enabled, latent, out)
        return list(out)

    def test_reference_state_keeps_pPv_direct_density(self):
        out = self.response(.6)
        self.assertAlmostEqual(out[0], .5, places=6)
        self.assertAlmostEqual(out[1], 0, places=6)
        self.assertAlmostEqual(out[2], .5, places=6)
        self.assertAlmostEqual(out[3], -1.5, places=5)
        self.assertEqual(out[4], 0)
        self.assertAlmostEqual(out[5], 1.5, places=6)
        self.assertAlmostEqual(out[6], .1, places=6)
        self.assertAlmostEqual(self.response(.6, latent=1)[4], -.01, places=8)

    def test_clapeyron_sign_and_buoyancy_residual(self):
        for T in (.55, .65):
            out = self.response(T)
            expected = .5*(1+math.tanh(-10*(T-.6)))
            self.assertAlmostEqual(out[0], expected, places=6)
            self.assertAlmostEqual(out[3], -2*(expected-.5)-3*expected, places=5)
        off = self.response(.6, enabled=0)
        self.assertEqual([off[i] for i in (0, 2, 4, 5, 6)], [0, 0, 0, 1, 1])
        self.assertAlmostEqual(off[3], 0, places=6)


if __name__ == "__main__":
    unittest.main()
