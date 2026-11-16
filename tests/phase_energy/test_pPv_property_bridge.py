"""Exercise the actual Phase Python caller, C setter and typed readers.

This is a host bridge test, not a claim that the cluster Python 2/Pyre runtime
has been exercised. Python 3 aliases only obsolete Python 2 C API names. A
small declaration collector reads defaults from the real Phase.py inventory;
the production Python setProperties method and C setter/getters run unchanged.
"""
import copy
import ctypes
import importlib.util
import os
from pathlib import Path
import shutil
import subprocess
import sys
import sysconfig
import tempfile
from types import ModuleType, SimpleNamespace
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[2]
VECTORS = ("phase_depth", "phase_delta_rho", "phase_delta_s", "phase_Ra",
           "phase_width", "phase_clapeyron", "phase_transT")


def function(source, signature):
    start = source.index(signature)
    opening = source.index("{", start)
    depth, end = 1, opening + 1
    while depth:
        depth += (source[end] == "{") - (source[end] == "}")
        end += 1
    return source[start:end]


def declared_phase_class():
    """Load real declarations and methods without replacing their contents."""
    inventory = ModuleType("pyre.inventory")

    def declaration(kind):
        def collect(name, default):
            return SimpleNamespace(name=name, default=copy.deepcopy(default), kind=kind)
        return collect

    for kind in ("array", "bool", "float"):
        setattr(inventory, kind, declaration(kind))
    pyre = ModuleType("pyre")
    pyre.inventory = inventory
    component = ModuleType("CitcomComponent")
    component.CitcomComponent = type("CitcomComponent", (), {"Inventory": type("Inventory", (), {})})
    namespace = {"__name__": "phase_inventory_under_test"}
    with patch.dict(sys.modules, {"pyre": pyre, "pyre.inventory": inventory,
                                 "CitcomComponent": component}):
        path = ROOT / "CitcomS/Components/Phase.py"
        exec(compile(path.read_text(), str(path), "exec"), namespace)
    return namespace["Phase"]


class PPvPropertyBridge(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        compiler = shutil.which(os.environ.get("MPICC", "mpicc"))
        if compiler is None:
            raise RuntimeError("mpicc is required for the production bridge test")
        cls.tmp = tempfile.TemporaryDirectory(prefix="pPv-property-bridge-")
        cls.addClassCleanup(cls.tmp.cleanup)
        directory = Path(cls.tmp.name)
        source = (ROOT / "module/setProperties.c").read_text()
        macros = source[source.index("FILE *get_output_stream"):source.index("void myerror")]
        getter_start = source.index("#define getTYPEProperty _getIntProperty")
        getters = source[getter_start:source.index("/* $Id:", getter_start)]
        helper_source = source[source.index("/* helper functions */"):]
        phase_source = (ROOT / "lib/Phase_change.c").read_text()
        setter = function(source, "PyObject * pyCitcom_Phase_set_properties(")
        c = r'''
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "phase_change.h"
#define PyInt_Check PyLong_Check
#define PyInt_AsLong PyLong_AsLong
#define PyInt_FromLong PyLong_FromLong
#define PyString_AsStringAndSize PyBytes_AsStringAndSize
#define PyCObject_AsVoidPtr(p) PyCapsule_GetPointer(p,NULL)
#define PyFile_Check(p) 0
#define PyFile_AsFile(p) NULL
''' + macros + getters + "\n" + function(
            helper_source, "FILE *get_output_stream(") + "\n" + function(
                phase_source, "void phase_change_set_pPv(") + "\n" + setter + r'''
static PyObject *call(PyObject *self,PyObject *properties) {
    struct All_variables *E=calloc(1,sizeof(*E));
    PyObject *capsule=PyCapsule_New(E,NULL,NULL);
    PyObject *args=PyTuple_Pack(3,capsule,properties,Py_None);
    PyObject *result=pyCitcom_Phase_set_properties(NULL,args),*out=NULL;
    int i;
    if(result) {
        PyObject *phases=PyList_New(PHASE_TRANSITIONS);
        for(i=0;i<PHASE_TRANSITIONS;i++) {
            const struct Phase_transition *p=&E->control.phase[i];
            PyList_SET_ITEM(phases,i,Py_BuildValue("(ddddddd)",
                (double)p->depth,(double)p->density_jump,(double)p->entropy_jump,
                (double)p->Ra,(double)p->clapeyron,(double)p->transT,(double)p->inv_width));
        }
        out=Py_BuildValue("{s:i,s:i,s:d,s:d,s:O,s:O}",
            "enabled",E->control.pPv_enabled,"latent",E->control.pPv_latent_enabled,
            "k_ratio",(double)E->control.pPv_k_ratio,
            "eta_ratio",(double)E->control.pPv_eta_ratio,"phases",phases,
            "capability",result);
        Py_DECREF(phases);
    }
    Py_XDECREF(result);Py_DECREF(args);Py_DECREF(capsule);free(E);
    return out;
}
static PyMethodDef methods[]={{"call",call,METH_O,"production phase setter"},{NULL}};
static struct PyModuleDef module={PyModuleDef_HEAD_INIT,"pPv_property_bridge",NULL,-1,methods};
PyMODINIT_FUNC PyInit_pPv_property_bridge(void){return PyModule_Create(&module);}
'''
        source_path = directory / "bridge.c"
        source_path.write_text(c)
        extension_path = directory / ("pPv_property_bridge" + sysconfig.get_config_var("EXT_SUFFIX"))
        command = [compiler, "-std=gnu99", "-shared", "-fPIC",
                   "-I" + sysconfig.get_path("include"), "-I" + str(ROOT / "lib"),
                   "-I" + str(ROOT / "module"), str(source_path), "-o", str(extension_path)]
        if sys.platform == "darwin":
            command += ["-undefined", "dynamic_lookup"]
        subprocess.run(command, check=True, capture_output=True, text=True)
        spec = importlib.util.spec_from_file_location("pPv_property_bridge", extension_path)
        cls.library = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.library)
        cls.phase_class = declared_phase_class()
        cls.declarations = {
            name: value for name, value in vars(cls.phase_class.Inventory).items()
            if isinstance(value, SimpleNamespace) and hasattr(value, "default")
        }

    def inventory(self, **overrides):
        values = {name: copy.deepcopy(value.default)
                  for name, value in self.declarations.items()}
        values.update(overrides)
        return SimpleNamespace(**values)

    def through_python_caller(self, inventory):
        observed = []
        bridge = ModuleType("CitcomSLib")

        def capture(state, properties, stream):
            snapshot = self.library.call(properties)
            observed.append(snapshot)
            return snapshot["capability"]

        bridge.Phase_set_properties = capture
        instance = SimpleNamespace(all_variables=None, inventory=inventory)
        with patch.dict(sys.modules, {"CitcomSLib": bridge}):
            result = self.phase_class.setProperties(instance, None)
        self.assertIsNone(result)
        self.assertEqual(len(observed), 1)
        return observed[0]

    @staticmethod
    def expected_phase(inventory, index):
        f = lambda value: ctypes.c_float(value).value
        width = f(inventory.phase_width[index])
        return (f(inventory.phase_depth[index]), f(inventory.phase_delta_rho[index]),
                f(inventory.phase_delta_s[index]), f(inventory.phase_Ra[index]),
                f(inventory.phase_clapeyron[index]), f(inventory.phase_transT[index]),
                f(1 / width) if width else 0.)

    def test_real_inventory_defaults_preserve_three_legacy_entries(self):
        inventory = self.inventory()
        for name in VECTORS:
            self.assertEqual(self.declarations[name].kind, "array")
            self.assertEqual(self.declarations[name].name, name)
            self.assertEqual(len(getattr(inventory, name)), 3)
        self.assertEqual(inventory.phase_depth, [410 / 6371, 520 / 6371, 660 / 6371])
        expected = {"pPv_enabled": False, "pPv_latent_enabled": False,
                    "pPv_k_ratio": 1., "pPv_eta_ratio": 1., "pPv_depth": 0.,
                    "pPv_transT": 0., "pPv_clapeyron": 0., "pPv_width": 0.,
                    "pPv_delta_rho": 0., "pPv_delta_s": 0.}
        for name, value in expected.items():
            self.assertEqual(getattr(inventory, name), value)
            self.assertEqual(self.declarations[name].name, name)
        state = self.through_python_caller(inventory)
        self.assertEqual(state["capability"], 1)
        self.assertEqual([state[k] for k in ("enabled", "latent", "k_ratio", "eta_ratio")],
                         [0, 0, 1., 1.])
        self.assertEqual(state["phases"][:3], [self.expected_phase(inventory, i) for i in range(3)])
        self.assertEqual(state["phases"][3], (0.,) * 7)

    def test_all_scalar_properties_and_independent_anchor_temperature(self):
        inventory = self.inventory(
            phase_depth=[.0625, .125, .1875], phase_delta_rho=[180., 58., 459.],
            phase_delta_s=[-.03125, -.015625, .0078125], phase_Ra=[2., 4., 8.],
            phase_width=[.00390625, .0078125, .015625],
            phase_clapeyron=[.0625, .125, -.03125], phase_transT=[.125, .375, .625],
            pPv_enabled=True, pPv_latent_enabled=True, pPv_k_ratio=1.75,
            pPv_eta_ratio=.03125, pPv_depth=.4375, pPv_transT=.8125,
            pPv_clapeyron=.1875, pPv_width=.00390625,
            pPv_delta_rho=88.25, pPv_delta_s=-.125)
        state = self.through_python_caller(inventory)
        self.assertEqual(state["capability"], 1)
        self.assertEqual([state[k] for k in ("enabled", "latent", "k_ratio", "eta_ratio")],
                         [1, 1, 1.75, .03125])
        self.assertEqual(state["phases"][:3], [self.expected_phase(inventory, i) for i in range(3)])
        # pPv density remains kg/m^3, entropy remains the supplied Cp0-normalized
        # scalar, and its temperature anchor is independent of all legacy T0s.
        self.assertEqual(state["phases"][3], (.4375, 88.25, -.125, 0., .1875, .8125, 256.))

    def test_wrong_legacy_vector_lengths_are_rejected(self):
        for name in VECTORS:
            for length in (2, 4):
                with self.subTest(property=name, length=length):
                    inventory = self.inventory(**{name: [.25] * length})
                    with self.assertRaisesRegex(ValueError, name):
                        self.through_python_caller(inventory)

    def test_missing_scalar_is_not_silently_ignored(self):
        inventory = self.inventory()
        del inventory.pPv_transT
        with self.assertRaisesRegex(AttributeError, "pPv_transT"):
            self.through_python_caller(inventory)

    def test_old_c_setter_rejected_only_when_pPv_is_enabled(self):
        for capability in (None, 0):
            bridge = ModuleType("CitcomSLib")
            bridge.Phase_set_properties = lambda *args: capability
            for enabled in (False, True):
                instance = SimpleNamespace(all_variables=None,
                                           inventory=self.inventory(pPv_enabled=enabled))
                with self.subTest(capability=capability, enabled=enabled):
                    with patch.dict(sys.modules, {"CitcomSLib": bridge}):
                        if enabled:
                            with self.assertRaisesRegex(RuntimeError, "pPv"):
                                self.phase_class.setProperties(instance, None)
                        else:
                            self.assertIsNone(self.phase_class.setProperties(instance, None))

    def test_real_typed_getters_reject_incorrect_types(self):
        for name, value in (("pPv_enabled", 1.0), ("pPv_k_ratio", "1.5"),
                            ("phase_width", [.01, "bad", .01])):
            with self.subTest(property=name):
                with self.assertRaisesRegex(ValueError, name):
                    self.through_python_caller(self.inventory(**{name: value}))


if __name__ == "__main__":
    unittest.main()
