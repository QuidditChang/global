"""Compile the production TA property setter/getters against local Python 3.
Uses compatibility shims for Python 2 string/capsule APIs. This exercises the
production transfer function, not the full Python 2/Pyre installation.
"""
from pathlib import Path
import re,json,subprocess,sysconfig,sys
root=Path(__file__).resolve().parents[2];s=(root/'module/setProperties.c').read_text()
head=s[s.index('#include <Python.h>'):s.index('char pyCitcom_Advection_diffusion_set_properties__doc__')]
fun=s[s.index('PyObject * pyCitcom_Advection_diffusion_set_properties('):s.index('char pyCitcom_BC_set_properties__doc__')]
helpers=s[s.index('int _getStringProperty(PyObject* properties, char* attribute,',s.index('FILE *get_output_stream(PyObject *out, struct All_variables*E)\n{')):]
defaults={name:__import__('ast').literal_eval(val) for name,val in re.findall(r'prop\.\w+\("(\w+)", default=([^\n]+)\)',(root/'CitcomS/Components/Advection_diffusion/Advection_diffusion.py').read_text())}
setup='from types import SimpleNamespace\np = SimpleNamespace(**'+repr(defaults)+')\np.energy_solver="pices"\np.pices_p4=True\n'
shim='''#define PyCObject_AsVoidPtr(o) PyCapsule_GetPointer(o,NULL)
#define PyInt_Check PyLong_Check
#define PyInt_AsLong PyLong_AsLong
static int PyString_AsStringAndSize(PyObject *p,char **b,Py_ssize_t *n){*b=(char*)PyUnicode_AsUTF8AndSize(p,n);return *b?0:-1;}
'''
main=r'''
FILE *get_output_stream(PyObject *out,struct All_variables *E){return stdout;}
void pices_fail(struct All_variables *E,const char *s){fprintf(stderr,"%s\n",s);exit(99);}
int main(void){
 static struct All_variables E; PyObject *p,*args,*answer,*cap;
 Py_Initialize();
 SETUP
 p=PyDict_GetItemString(PyModule_GetDict(PyImport_AddModule("__main__")),"p");
 cap=PyCapsule_New(&E,NULL,NULL);args=PyTuple_Pack(3,cap,p,Py_None);
 for(int flag=0;flag<=1;flag++){
 PyObject_SetAttrString(p,"pices_ta_before_heat",flag?Py_True:Py_False);
 answer=pyCitcom_Advection_diffusion_set_properties(NULL,args);
 if(!answer){PyErr_Print();return 2;}
 if(E.pices.ta_before_heat!=flag)return 3;
 Py_DECREF(answer);
 }
 PyObject_SetAttrString(p,"pices_p4",Py_False);
 answer=pyCitcom_Advection_diffusion_set_properties(NULL,args);
 if(answer || !PyErr_ExceptionMatches(PyExc_ValueError))return 4;
 PyErr_Clear();puts("PROPERTY_TRANSFER_PASS");return 0;
}
'''.replace('SETUP','if(PyRun_SimpleString('+json.dumps(setup)+')!=0)return 1;')
d=Path('/tmp/pices-property-test');d.mkdir(exist_ok=True);f=d/'test.c';f.write_text(head.replace('#include <Python.h>','#include <Python.h>\n'+shim)+fun+helpers+main)
cmd=['/usr/local/bin/mpicc','-std=gnu99','-I'+str(root/'lib'),'-I'+str(root/'module'),'-I'+sysconfig.get_paths()['include'],str(f),'-L'+sysconfig.get_config_var('LIBDIR'),'-lpython'+sysconfig.get_config_var('VERSION'),'-lm','-o',str(d/'test')]
subprocess.run(cmd,check=True);x=subprocess.run([str(d/'test')],capture_output=True,text=True);(d/'output.txt').write_text(x.stdout+x.stderr);assert x.returncode==0,(x.returncode,x.stderr);assert 'pices_ta_before_heat=0' in x.stdout and 'pices_ta_before_heat=1' in x.stdout;print('PASS: real production setter + real property getters, flag off/on, generated flag, invalid combination')
