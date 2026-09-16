"""Exercise production Pyre bridge/settings/fatal diagnostics using a C harness.

Python API/property lookup is mocked (no legacy Python 2/Pyre installation needed).
The full viscosity bridge function and rheol7 diagnostic helpers are extracted
unchanged from production sources. MPI termination uses two real processes.
Run: python3 tests/test_rheol7_runtime.py
"""
from pathlib import Path
import subprocess
import tempfile
import os
import shutil

root = Path(__file__).resolve().parents[1]
mpicc = shutil.which(os.environ.get('MPICC', 'mpicc'))
assert mpicc, 'MPI C compiler required'
# Keep launcher and compiler from the same installation (e.g. avoid ParaView MPI).
mpiexec = os.environ.get('MPIEXEC', str(Path(mpicc).resolve().with_name('mpiexec')))
visc = (root / 'lib/Viscosity_structures.c').read_text()
bridge = (root / 'module/setProperties.c').read_text()
helpers = visc[visc.index('void validate_rheol7_settings('):
               visc.index('static double strict_rheology_reference_temperature(')]
bridge = bridge[bridge.index('PyObject * pyCitcom_Visc_set_properties('):
                bridge.index('char pyCitcom_Incompressible_set_properties__doc__')]
prefix = r'''
#include <assert.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "Steinberger_nuref.h"
typedef void PyObject;
static struct All_variables state;
static float requested_scale;
static int float_reads;
#define Py_None ((PyObject *)&state)
#define Py_INCREF(x) ((void)0)
#define PyArg_ParseTuple(args,fmt,obj,properties,out) \
    (*(obj)=Py_None, *(properties)=Py_None, *(out)=Py_None, 1)
#define PyCObject_AsVoidPtr(obj) (&state)
#define get_output_stream(out,E) stdout
#define PUTS(s) fprintf(fp,s)
#define getStringProperty(p,k,v,f) strcpy(v,"system")
#define getIntProperty(p,k,v,f) \
    (v = !strcmp(k,"rheol") ? 7 : (!strcmp(k,"TDEPV") || !strcmp(k,"num_mat")))
#define getFloatVectorProperty(p,k,v,n,f) ((void)0)
#define getFloatProperty(p,k,v,f) do { \
    v = !strcmp(k,"cold_scale") ? requested_scale : 0; \
    if(!strcmp(k,"cold_scale")) { float_reads++; fprintf(f,"cold_scale=%.9g\n",v); } \
} while(0)
void myerror(struct All_variables *E, char *s) { abort(); }
'''
main = r'''
int main(int argc, char **argv) {
    struct All_variables *E = &state;
    int rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    E->parallel.world = MPI_COMM_WORLD;
    E->parallel.me = rank;
    requested_scale = strtod(argv[1],NULL);
    E->viscosity.cold_scale = NAN; /* malloc state must be overwritten. */
    if(!strcmp(argv[2],"invalid") && rank == 0) {
        MPI_Barrier(MPI_COMM_WORLD);
        abort();
    }
    assert(pyCitcom_Visc_set_properties(NULL,NULL) == Py_None);
    assert(float_reads == 1 && E->viscosity.cold_scale == requested_scale);
    if(!strcmp(argv[2],"fatal") && rank == 1) {
        E->monitor.solution_cycles = 4;
        E->data.Ttop = 300.; E->data.ref_temperature = 3400.;
        E->data.ref_viscosity = 1e21;
        E->mesh.nsd = 3;
        E->ien[1] = calloc(2,sizeof(*E->ien[1]));
        E->T[1] = calloc(9,sizeof(*E->T[1]));
        for(int a=1;a<=8;a++) E->ien[1][1].node[a] = a;
        rheol7_failure(E,"test invalid temperature",1,1,1,7.,-0.1,1600.,-1.);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
'''
with tempfile.TemporaryDirectory() as tmp:
    path = Path(tmp)
    source, exe = path / 'runtime.c', path / 'runtime'
    source.write_text(prefix + helpers + bridge + main)
    subprocess.run([mpicc, '-std=gnu99', '-I'+str(root/'lib'), str(source),
                    '-lm', '-o', str(exe)], check=True)
    for scale in ('0', '0.25', '0.5', '1'):
        result = subprocess.run([str(exe), scale, 'normal'], capture_output=True,
                                text=True, timeout=20, check=True)
        assert 'cold_scale='+scale in result.stdout
        assert 'rheol=7: TDEPV=1 cold_scale='+scale in result.stderr
    for scale in ('-1', 'nan', 'inf'):
        result = subprocess.run([mpiexec, '-n', '2', str(exe), scale, 'invalid'],
                                capture_output=True, text=True, timeout=20)
        assert result.returncode != 0
        assert 'cold_scale must be finite and >= 0' in result.stderr, result.stderr
    result = subprocess.run([mpiexec, '-n', '2', str(exe), '0.5', 'fatal'],
                            capture_output=True, text=True, timeout=20)
    assert result.returncode != 0
    for token in ('rank=1 step=4 cap=1 element=1 gp=1', 'T_K=', 'Tref_K=',
                  'cold_scale=0.5', 'ln_eta=', 'node=8', 'rheol7_clip_diag_v1',
                  'Ttop_K=300', 'DeltaT_K=3400', 'clip_lower_nd=0',
                  'bounds_unordered=0', 'clipped_T_nd=0', 'weight_nonfinite=0'):
        assert token in result.stderr, result.stderr
print('PASS: production bridge transfers 0/0.25/0.5/1; startup logging; '
      'invalid settings and runtime failure abort both MPI ranks with diagnostics.')
