"""Exercise unchanged production Pyre/standalone readers and MPI diagnostics.

Python API/property lookup is mocked without requiring legacy Python2/Pyre.
The vector helper, viscosity bridge, standalone reader and runtime diagnostics
are extracted from production. Invalid runtime settings abort two real ranks.
"""
from pathlib import Path
import os
import shutil
import subprocess
import tempfile

root=Path(__file__).resolve().parents[1]
mpicc=shutil.which(os.environ.get('MPICC','mpicc'))
assert mpicc, 'MPI C compiler required'
mpiexec=os.environ.get('MPIEXEC',str(Path(mpicc).resolve().with_name('mpiexec')))
visc=(root/'lib/Viscosity_structures.c').read_text()
bridge=(root/'module/setProperties.c').read_text()
helpers=visc[visc.index('void validate_rheol7_settings('):visc.index('static double strict_rheology_reference_temperature(')]
reader=visc[visc.index('static void read_rheol7_vector('):visc.index('void viscosity_system_input(')]
bridge=bridge[bridge.index('static int rheol7_get_vector('):bridge.index('char pyCitcom_Incompressible_set_properties__doc__')]
prefix=r'''
#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <stdarg.h>
#include <string.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "Steinberger_nuref.h"
typedef struct FakeObject {int size,bad;double value;struct FakeObject *items[4];} PyObject;
static struct All_variables state;
static float requested_scale;
static int requested_mode,requested_tdep,float_reads,vector_reads,py_error;
static const char *scenario,*requested_text;
static PyObject vectors[4],numbers[4][4];
#define Py_None ((PyObject *)&state)
#define PyExc_ValueError ((PyObject *)&state)
#define Py_INCREF(x) ((void)0)
#define Py_DECREF(x) ((void)0)
#define PyArg_ParseTuple(args,fmt,obj,properties,out) (*(obj)=Py_None,*(properties)=Py_None,*(out)=Py_None,1)
#define PyCObject_AsVoidPtr(obj) (&state)
#define get_output_stream(out,E) stdout
#define PUTS(s) fprintf(fp,s)
#define getStringProperty(p,k,v,f) strcpy(v,"system")
#define getIntProperty(p,k,v,f) (v=!strcmp(k,"rheol")?requested_mode:(!strcmp(k,"TDEPV")?requested_tdep:!strcmp(k,"num_mat")))
#define getFloatVectorProperty(p,k,v,n,f) ((void)0)
#define getFloatProperty(p,k,v,f) do {v=!strcmp(k,"cold_scale")?requested_scale:0; \
 if(!strcmp(k,"cold_scale")){float_reads++;fprintf(f,"cold_scale=%.9g\n",v);}}while(0)
static PyObject *PyObject_GetAttrString(PyObject *unused,const char *key) {
 int k=!strcmp(key,"rheol7_layer_factors")?0:(!strcmp(key,"rheol7_poly_a")?1:(!strcmp(key,"rheol7_poly_b")?2:3));
 vector_reads++;
 if(!strcmp(scenario,"py_missing")){py_error=1;return NULL;}
 return &vectors[k];
}
static PyObject *PySequence_Fast(PyObject *object,const char *message) {
 if(!strcmp(scenario,"py_sequence")){py_error=1;return NULL;}return object;
}
#define PySequence_Fast_GET_SIZE(obj) ((obj)->size)
#define PySequence_Fast_GET_ITEM(obj,i) ((obj)->items[i])
static double PyFloat_AsDouble(PyObject *obj) {if(obj->bad)py_error=1;return obj->value;}
#define PyErr_Occurred() (py_error)
static PyObject *PyErr_Format(PyObject *kind,const char *format,...) {
 va_list args;py_error=1;va_start(args,format);vfprintf(stderr,format,args);va_end(args);return NULL;
}
static int input_string(char *name,char *output,char *fallback,int rank) {
 if(!strcmp(requested_text,"missing"))return 0;
 strcpy(output,requested_text);return 1;
}
void myerror(struct All_variables *E,char *s){abort();}
'''
main=r'''
int main(int argc,char **argv) {
 struct All_variables *E=&state;int rank,i,k;
 MPI_Init(&argc,&argv);MPI_Comm_rank(MPI_COMM_WORLD,&rank);
 E->parallel.world=MPI_COMM_WORLD;E->parallel.me=rank;
 requested_mode=atoi(argv[1]);requested_scale=strtod(argv[2],NULL);scenario=argv[3];
 requested_tdep=strcmp(scenario,"off")!=0;
 for(k=0;k<4;k++) {
  vectors[k].size=4;
  if(!strcmp(scenario,"py_short"))vectors[k].size=3;
  if(!strcmp(scenario,"py_long"))vectors[k].size=5;
  for(i=0;i<4;i++) {
   numbers[k][i].value=k<2?1.123456789012345+i:0.;
   numbers[k][i].bad=!strcmp(scenario,"py_nonnumeric");
   if(!strcmp(scenario,"py_nonfinite"))numbers[k][i].value=NAN;
   vectors[k].items[i]=&numbers[k][i];
  }
 }
 if(!strcmp(scenario,"invalid_factor")) {
  if(requested_mode==71)numbers[0][2].value=0.;
  else {numbers[1][0].value=1.;numbers[2][0].value=-8.;numbers[3][0].value=8.;}
 }
 if(!strcmp(scenario,"invalid") || !strcmp(scenario,"invalid_factor") || !strcmp(scenario,"parse_invalid")) {
  if(rank==0){MPI_Barrier(MPI_COMM_WORLD);abort();}
 }
 if(!strncmp(scenario,"parse",5)) {
  E->viscosity.RHEOL=requested_mode;requested_text=argv[4];
  steinberger_defaults(E->viscosity.rheol7_layer_factors,E->viscosity.rheol7_poly_a,E->viscosity.rheol7_poly_b,E->viscosity.rheol7_poly_c);
  read_rheol7_vector(E,"rheol7_layer_factors",E->viscosity.rheol7_layer_factors);
  if(!strcmp(requested_text,"missing"))assert(E->viscosity.rheol7_layer_factors[0]==1.);
  else {assert(E->viscosity.rheol7_layer_factors[0]==1.123456789012345);assert(E->viscosity.rheol7_layer_factors[3]==4.);}
 } else {
  E->viscosity.cold_scale=NAN;
  PyObject *result=pyCitcom_Visc_set_properties(NULL,NULL);
  if(!strncmp(scenario,"py_",3) && requested_mode!=7) {
   assert(result==NULL && py_error);MPI_Finalize();return 0;
  }
  assert(result==Py_None && !py_error);
  assert(float_reads==1 && E->viscosity.cold_scale==requested_scale);
  assert(vector_reads==(requested_mode==7?0:(requested_mode==71?1:3)));
  if(requested_mode==71)assert(E->viscosity.rheol7_layer_factors[0]==1.123456789012345);
  if(requested_mode==72)assert(E->viscosity.rheol7_poly_a[3]==4.123456789012345);
  if(!strcmp(scenario,"fatal") && rank==1) {
   E->monitor.solution_cycles=4;E->data.Ttop=300.;E->data.ref_temperature=3400.;E->data.ref_viscosity=1e21;
   E->mesh.nsd=3;E->ien[1]=calloc(2,sizeof(*E->ien[1]));E->T[1]=calloc(9,sizeof(*E->T[1]));
   for(i=1;i<=8;i++)E->ien[1][1].node[i]=i;
   rheol7_failure(E,"test invalid temperature",1,1,1,7.,-.1,1600.,-1.);
  }
 }
 MPI_Barrier(MPI_COMM_WORLD);MPI_Finalize();return 0;
}
'''
with tempfile.TemporaryDirectory() as temporary:
    path=Path(temporary);exe=path/'runtime'
    (path/'runtime.c').write_text(prefix+helpers+reader+bridge+main)
    subprocess.run([mpicc,'-std=gnu99','-I'+str(root/'lib'),str(path/'runtime.c'),'-lm','-o',str(exe)],check=True)
    def run(mode,scale,scenario,*extra,mpi=False):
        command=([mpiexec,'-n','2'] if mpi else [])+[str(exe),str(mode),str(scale),scenario]+list(extra)
        return subprocess.run(command,capture_output=True,text=True,timeout=30)
    for mode in (7,71,72):
        for scale in ('0','0.25','0.5','1'):
            result=run(mode,scale,'normal');assert result.returncode==0,result.stderr
            assert f'rheol={mode}: TDEPV=1 cold_scale={scale}' in result.stderr
        result=run(mode,'.5','off');assert result.returncode==0,result.stderr
        assert f'rheol={mode}: TDEPV=0' in result.stderr
        result=run(mode,'-1','invalid',mpi=True)
        assert result.returncode!=0 and 'cold_scale must be finite and >= 0' in result.stderr,result.stderr
        result=run(mode,'.5','fatal',mpi=True);assert result.returncode!=0
        for token in (f'rheol={mode} fatal:','rank=1 step=4 cap=1 element=1 gp=1','T_K=','Tref_K=',
                      'cold_scale=0.5','ln_eta=','node=8','rheol7_clip_diag_v1','Ttop_K=300',
                      'DeltaT_K=3400','clip_lower_nd=0','bounds_unordered=0','clipped_T_nd=0'):
            assert token in result.stderr,result.stderr
    for mode in (71,72):
        for scenario in ('py_short','py_long','py_missing','py_sequence','py_nonnumeric','py_nonfinite'):
            result=run(mode,'.5',scenario);assert result.returncode==0,result.stderr
        result=run(mode,'.5','invalid_factor',mpi=True)
        assert result.returncode!=0 and 'f(z)>0 throughout' in result.stderr,result.stderr
    result=run(7,'.5','py_missing');assert result.returncode==0,result.stderr
    for scale in ('nan','inf'):
        result=run(72,scale,'invalid',mpi=True)
        assert result.returncode!=0 and 'cold_scale must be finite' in result.stderr,result.stderr
    for value in ('missing','1.123456789012345, 2, 3, 4'):
        result=run(71,'.5','parse',value);assert result.returncode==0,result.stderr
    for value in ('1,2,3','1,2,3,4,5','1,2,nan,4','1,2,3,4junk','1 2 3 4','1,2,3,'):
        result=run(71,'.5','parse_invalid',value,mpi=True)
        assert result.returncode!=0 and 'exactly four finite' in result.stderr,result.stderr
print('PASS: all-mode Pyre bridge transfers double coefficients and cold_scale; malformed lists rejected; '
      'standalone exact-four parser; TDEPV on/off startup; invalid settings/runtime failures abort both MPI ranks.')
