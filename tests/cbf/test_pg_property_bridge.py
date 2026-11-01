"""Exercise the production Python caller, C setter, and typed property readers.
Host Python 3 aliases only the obsolete Python 2 C API names; no setter/getter
logic is replaced. The actual cluster Python 2/Pyre stack is checked by the
short HPC run and the startup handshake, not claimed by this host test.
"""
import ast, importlib.util, re, subprocess, sys, sysconfig, tempfile, textwrap
from pathlib import Path
from types import SimpleNamespace, ModuleType
import unittest
from unittest.mock import patch
from test_cbf_kernel import function
ROOT=Path(__file__).resolve().parents[2]

class PropertyBridge(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp=tempfile.TemporaryDirectory();d=Path(cls.tmp.name)
        source=(ROOT/'module/setProperties.c').read_text()
        macros=source[source.index('FILE *get_output_stream'):source.index('void myerror')]
        getters=source[source.index('#define getTYPEProperty _getIntProperty'):source.index('/* $Id:',source.index('#define getTYPEProperty _getIntProperty'))]
        c=r'''
#include <Python.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "element_definitions.h"
#include "global_defs.h"
#define PyInt_Check PyLong_Check
#define PyInt_AsLong PyLong_AsLong
#define PyInt_FromLong PyLong_FromLong
#define PyString_AsStringAndSize PyBytes_AsStringAndSize
#define PyCObject_AsVoidPtr(p) PyCapsule_GetPointer(p,NULL)
#define PyFile_Check(p) 0
#define PyFile_AsFile(p) NULL
'''+macros+getters+'\n'+function(source[source.index('/* helper functions */'):],'FILE *get_output_stream(')+'\n'+function(source[source.index('/* helper functions */'):],'int _getStringProperty(')+'\n'+function(source,'PyObject * pyCitcom_Solver_set_properties(')+r'''
static PyObject *call(PyObject *self,PyObject *props) {
 struct All_variables *E=calloc(1,sizeof(*E));
 PyObject *cap=PyCapsule_New(E,NULL,NULL),*args=PyTuple_Pack(3,cap,props,Py_None);
 PyObject *result=pyCitcom_Solver_set_properties(NULL,args),*out=NULL;
 if(result)out=Py_BuildValue("Oi",result,E->control.pg_ta_before_heat);
 Py_XDECREF(result);Py_DECREF(args);Py_DECREF(cap);free(E);return out;
}
static PyMethodDef methods[]={{"call",call,METH_O,"production setter"},{NULL}};
static struct PyModuleDef mod={PyModuleDef_HEAD_INIT,"pg_bridge",NULL,-1,methods};
PyMODINIT_FUNC PyInit_pg_bridge(void){return PyModule_Create(&mod);}
'''
        src=d/'bridge.c';src.write_text(c);so=d/('pg_bridge'+sysconfig.get_config_var('EXT_SUFFIX'))
        cmd=['mpicc','-std=gnu99','-shared','-fPIC','-I'+sysconfig.get_path('include'),'-I'+str(ROOT/'lib'),'-I'+str(ROOT/'module'),str(src),'-o',str(so)]
        if sys.platform=='darwin':cmd+=['-undefined','dynamic_lookup']
        subprocess.run(cmd,check=True,capture_output=True,text=True)
        spec=importlib.util.spec_from_file_location('pg_bridge',so);cls.lib=importlib.util.module_from_spec(spec);spec.loader.exec_module(cls.lib)
        s=(ROOT/'CitcomS/Solver/Solver.py').read_text();start=s.index('    def setProperties(self, stream):');end=s.index('\n    def ',start+10)
        ns={};exec(textwrap.dedent(s[start:end]),ns);cls.caller=staticmethod(ns['setProperties'])
        cls.setter=function(source,'PyObject * pyCitcom_Solver_set_properties(')
    @classmethod
    def tearDownClass(cls):cls.tmp.cleanup()
    def inventory(self,on):
        values={}
        for typ,name in re.findall(r'get(String|Float|Double|Int)Property\(properties, "([^"]+)"',self.setter):
            values[name]=b'test' if typ=='String' else (0 if typ=='Int' else 0.0)
        values['pg_ta_before_heat']=on
        for name in ['mesher','tsolver','vsolver','bc','const','ic','output','param','phase','tracer','visc']:
            values[name]=SimpleNamespace(setProperties=lambda stream:None)
        return SimpleNamespace(**values)
    def test_actual_c_read_and_return(self):
        for on in [False,True]:self.assertEqual(self.lib.call(self.inventory(on)),(int(on),int(on)))
    def test_missing_property_is_fatal(self):
        inv=self.inventory(True);del inv.pg_ta_before_heat
        with self.assertRaises(AttributeError):self.lib.call(inv)
    def test_actual_python_to_c_chain(self):
        fake=ModuleType('CitcomSLib');fake.Solver_set_properties=lambda state,inv,stream:self.lib.call(inv)[0]
        obj=SimpleNamespace(all_variables=None,inventory=self.inventory(True),communicator=SimpleNamespace(rank=1))
        with patch.dict(sys.modules,{'CitcomSLib':fake}):self.caller(obj,None)
    def test_old_or_mismatched_binary_rejected(self):
        obj=SimpleNamespace(all_variables=None,inventory=self.inventory(True),communicator=SimpleNamespace(rank=1))
        for returned in [None,0]:
            fake=ModuleType('CitcomSLib');fake.Solver_set_properties=lambda *args:returned
            with patch.dict(sys.modules,{'CitcomSLib':fake}),self.assertRaises(RuntimeError):self.caller(obj,None)

if __name__=='__main__':unittest.main()
