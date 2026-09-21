"""Execute the production Pyre->C setter with a minimal Python API fixture.

This covers the path bypassed by standalone Instructions.c parsing. It is not
an end-to-end Python 2/Pyre installation test.
"""
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
ROOT=Path(__file__).resolve().parents[2]

class BindingTest(unittest.TestCase):
    def test_inventory_boolean_reaches_solver_control(self):
        source=(ROOT/'module/setProperties.c').read_text()
        start=source.index('PyObject * pyCitcom_Solver_set_properties(')
        pos=source.index('{',start)+1;depth=1;end=pos
        while depth:
            depth+=(source[end]=='{')-(source[end]=='}');end+=1
        setter=source[start:end]
        macros=source[source.index('#define PUTS'):source.index('/*==============================================================*/',source.index('#define PUTS'))]
        # Only the real get-property macro definitions are needed; prototypes
        # and surrounding helper implementations use additional Python APIs.
        macros='\n'.join(line for line in macros.splitlines() if line.startswith('#define get') or line.startswith('#define PUTS'))
        prefix=r'''
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include "global_defs.h"
typedef struct {int ignored;} PyObject;
static struct All_variables state;
static int requested;
static PyObject dummy;
#define Py_None (&dummy)
#define PyExc_ValueError (&dummy)
#define Py_INCREF(x) ((void)0)
static void PyErr_SetString(PyObject *e,const char *s) {abort();}
static int PyArg_ParseTuple(PyObject *args,const char *fmt,...) {
    va_list ap; int i; va_start(ap,fmt);
    for(i=0;i<3;i++) *va_arg(ap,PyObject **)=&dummy;
    va_end(ap);return 1;
}
static void *PyCObject_AsVoidPtr(PyObject *x) {return &state;}
static FILE *get_output_stream(PyObject *x,struct All_variables *E) {return NULL;}
static int _getIntProperty(PyObject *p,char *key,int *v,FILE *f) {
    *v=strcmp(key,"temperature_audit")==0 ? requested : 0;return 0;
}
static int _getFloatProperty(PyObject *p,char *key,float *v,FILE *f) {*v=0;return 0;}
static int _getDoubleProperty(PyObject *p,char *key,double *v,FILE *f) {*v=0;return 0;}
static int _getStringProperty(PyObject *p,char *key,char *v,size_t n,FILE *f) {v[0]=0;return 0;}
'''
        main=r'''
int main(void) {
    requested=1; state.control.temperature_audit=0;
    if(!pyCitcom_Solver_set_properties(&dummy,&dummy) || state.control.temperature_audit!=1) return 1;
    requested=0; state.control.temperature_audit=1;
    if(!pyCitcom_Solver_set_properties(&dummy,&dummy) || state.control.temperature_audit!=0) return 2;
    return 0;
}
'''
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);(p/'test.c').write_text(prefix+macros+'\n'+setter+main)
            subprocess.run([os.environ.get('MPICC','mpicc'),'-std=gnu99','-I'+str(ROOT/'lib'),str(p/'test.c'),'-lm','-o',str(p/'test')],check=True)
            subprocess.run([str(p/'test')],check=True)

if __name__=='__main__':unittest.main()
