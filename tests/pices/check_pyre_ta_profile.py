"""Exercise the production TA profile bridge and checkpoint fingerprint.

The bridge compiles the real Param setter/getters with Python 3 compatibility
shims for the installed Python 2 API. The fingerprint uses the real checkpoint
implementation; neither test replaces a full Pyre/MPI installation check.
"""
import ast
import json
import os
from pathlib import Path
import re
import subprocess
import sysconfig
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
COMPILER = os.environ.get("MPICC", "/usr/local/bin/mpicc")

SHIM = r'''
#define PyCObject_AsVoidPtr(o) PyCapsule_GetPointer(o,NULL)
#define PyInt_Check PyLong_Check
#define PyInt_AsLong PyLong_AsLong
static int PyString_AsStringAndSize(PyObject *p,char **b,Py_ssize_t *n) {
 *b=(char*)PyUnicode_AsUTF8AndSize(p,n);return *b?0:-1;
}
'''

BRIDGE_MAIN = r'''
FILE *get_output_stream(PyObject *out,struct All_variables *E){return stdout;}
int main(void) {
 static struct All_variables E;PyObject *p,*args,*answer,*cap;
 Py_Initialize();
 SETUP
 p=PyDict_GetItemString(PyModule_GetDict(PyImport_AddModule("__main__")),"p");
 cap=PyCapsule_New(&E,NULL,NULL);args=PyTuple_Pack(3,cap,p,Py_None);
 answer=pyCitcom_Param_set_properties(NULL,args);
 if(!answer){PyErr_Print();return 2;}Py_DECREF(answer);
 if(E.control.lith_age_asml_profile!=0 ||
    E.control.lith_age_asml_logistic_center!=(float)(80.0/127.42) ||
    E.control.lith_age_asml_logistic_width!=(float)(4.0/127.42))return 3;
 PyObject_SetAttrString(p,"lith_age_asml_profile",PyLong_FromLong(1));
 PyObject_SetAttrString(p,"lith_age_asml_logistic_center",PyFloat_FromDouble(.7));
 PyObject_SetAttrString(p,"lith_age_asml_logistic_width",PyFloat_FromDouble(.05));
 answer=pyCitcom_Param_set_properties(NULL,args);
 if(!answer){PyErr_Print();return 4;}Py_DECREF(answer);
 if(E.control.lith_age_asml_profile!=1 ||
    E.control.lith_age_asml_logistic_center!=(float).7 ||
    E.control.lith_age_asml_logistic_width!=(float).05)return 5;
 PyObject_SetAttrString(p,"lith_age_asml_profile",PyUnicode_FromString("wrong_type"));
 answer=pyCitcom_Param_set_properties(NULL,args);
 if(answer || !PyErr_Occurred())return 6;
 PyErr_Clear();puts("TA_PROFILE_BRIDGE_PASS");return 0;
}
'''

FINGERPRINT_MAIN = r'''
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "global_defs.h"
static int one_barrier(MPI_Comm world){return 0;}
static int one_allgather(const void *send,int n,MPI_Datatype a,
                        void *recv,int m,MPI_Datatype b,MPI_Comm world) {
 assert(n==65 && m==65);memcpy(recv,send,65);return 0;
}
#define MPI_Barrier one_barrier
#define MPI_Allgather one_allgather
void pices_fail(struct All_variables *E,const char *why) {
 fprintf(stderr,"PICES_ERROR: %s\n",why);exit(72);
}
int pices_motion_enabled(struct All_variables *E){return 0;}
float find_age_in_MY(struct All_variables *E){return 0;}
#include "Pices_checkpoint.c"
int main(void) {
 struct All_variables *E=calloc(1,sizeof(*E));int d;
 char legacy[65],unused[65],curve[65],center[65],width[65],restored[65];
 E->fp=stdout;E->lmesh.noz=E->lmesh.nox=E->lmesh.noy=2;
 E->lmesh.nno=8;E->lmesh.nel=1;E->lmesh.nxs=E->lmesh.nys=1;E->mesh.nox=2;
 E->pices.p4=E->pices.eba=1;E->control.RH=1;E->viscosity.RHEOL=1;
 E->advection.fixed_timestep=1;E->control.lith_age=E->control.lith_age_asml=1;
 E->control.lith_age_asml_logistic_center=(float)(80.0/127.42);
 E->control.lith_age_asml_logistic_width=(float)(4.0/127.42);
 E->refstate.rho=calloc(3,sizeof(double));
 E->refstate.thermal_expansivity=calloc(3,sizeof(double));
 E->refstate.gravity=calloc(3,sizeof(double));
 E->refstate.heat_capacity=calloc(3,sizeof(double));
 E->refstate.Tref=calloc(3,sizeof(double));
 E->age_t=calloc(5,sizeof(float));E->flag_depth2=calloc(5,sizeof(float));
 for(d=1;d<=3;d++) {
  E->x[1][d]=calloc(9,sizeof(double));
  E->sphere.cap[1].TB[d]=calloc(9,sizeof(float));
 }
 fingerprint(E,legacy,0);
 E->control.lith_age_asml_logistic_center=.7f;
 E->control.lith_age_asml_logistic_width=.05f;
 fingerprint(E,unused,0);assert(!strcmp(legacy,unused));
 E->control.lith_age_asml_profile=1;
 fingerprint(E,curve,0);assert(strcmp(legacy,curve));
 E->control.lith_age_asml_logistic_center=.6f;
 fingerprint(E,center,0);assert(strcmp(curve,center));
 E->control.lith_age_asml_logistic_center=.7f;
 E->control.lith_age_asml_logistic_width=.04f;
 fingerprint(E,width,0);assert(strcmp(curve,width));
 E->control.lith_age_asml_logistic_width=.05f;
 fingerprint(E,restored,0);assert(!strcmp(curve,restored));
 puts("TA_PROFILE_FINGERPRINT_PASS");return 0;
}
'''


class TAProfilePlumbing(unittest.TestCase):
    def test_real_param_bridge_defaults_values_and_generated_configuration(self):
        source = (ROOT / "module/setProperties.c").read_text()
        head = source[source.index("#include <Python.h>"):
                      source.index("char pyCitcom_Advection_diffusion_set_properties__doc__")]
        function = source[source.index("PyObject * pyCitcom_Param_set_properties("):
                          source.index("char pyCitcom_Phase_set_properties__doc__")]
        helpers = source[source.index(
            "int _getStringProperty(PyObject* properties, char* attribute,",
            source.index("FILE *get_output_stream(PyObject *out, struct All_variables*E)\n{")):]
        inventory = (ROOT / "CitcomS/Components/Param.py").read_text()
        defaults = {name: ast.literal_eval(value) for name, value in re.findall(
            r'pyre\.inventory\.\w+\("(\w+)",\s*default=([^\n]+)\)', inventory)}
        setup = "from types import SimpleNamespace\np=SimpleNamespace(**" + repr(defaults) + ")\n"
        main = BRIDGE_MAIN.replace(
            "SETUP", "if(PyRun_SimpleString(" + json.dumps(setup) + ")!=0)return 1;")
        code = head.replace("#include <Python.h>", "#include <Python.h>\n" + SHIM)
        with tempfile.TemporaryDirectory(prefix="ta-profile-bridge-") as folder:
            source_path, executable = Path(folder) / "bridge.c", Path(folder) / "bridge"
            source_path.write_text(code + function + helpers + main)
            flags = [COMPILER, "-std=gnu99", "-I" + str(ROOT / "lib"),
                     "-I" + str(ROOT / "module"), "-I" + sysconfig.get_paths()["include"]]
            subprocess.run(flags + [str(source_path), "-L" + sysconfig.get_config_var("LIBDIR"),
                           "-lpython" + sysconfig.get_config_var("VERSION"), "-lm",
                           "-o", str(executable)], check=True)
            run = subprocess.run([str(executable)], capture_output=True, text=True)
            self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
            for value in ("lith_age_asml_profile=0", "lith_age_asml_profile=1",
                          "lith_age_asml_logistic_center=0.7",
                          "lith_age_asml_logistic_width=0.05", "TA_PROFILE_BRIDGE_PASS"):
                self.assertIn(value, run.stdout)

    def test_real_checkpoint_fingerprint_tracks_only_active_curve_parameters(self):
        with tempfile.TemporaryDirectory(prefix="ta-profile-checkpoint-") as folder:
            source_path, executable = Path(folder) / "fingerprint.c", Path(folder) / "fingerprint"
            source_path.write_text(FINGERPRINT_MAIN)
            subprocess.run([COMPILER, "-std=gnu99", "-I" + str(ROOT / "lib"),
                            str(source_path), "-lm", "-o", str(executable)], check=True)
            run = subprocess.run([str(executable)], capture_output=True, text=True)
            self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
            self.assertIn("TA_PROFILE_FINGERPRINT_PASS", run.stdout)


if __name__ == "__main__":
    unittest.main()
