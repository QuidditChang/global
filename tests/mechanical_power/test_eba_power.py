"""Exercise production EBA work/NPZ code on a manufactured element.

Geometry and the stiffness action are small deterministic fixtures; actual
get_elt_f, gradient/divergence assembly, power bookkeeping and NPZ writing run.
This is a diagnostic regression, not a full spherical MPI solver validation.
"""
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[2]

def function(source, signature):
    start = source.index(signature)
    opening = source.index('{', start)
    end, depth = opening + 1, 1
    while depth:
        depth += (source[end] == '{') - (source[end] == '}')
        end += 1
    return source[start:end]

class MechanicalPowerTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.path = Path(cls.tmp.name)
        drive = (ROOT/'lib/Drive_solvers.c').read_text()
        element = (ROOT/'lib/Element_calculations.c').read_text()
        profile = (ROOT/'lib/Profile_output.c').read_text()
        sources = ['#include "fixture.h"']
        for signature in ['static double **allocate_nodal_field(',
                          'static double **allocate_element_field(',
                          'static double **allocate_equation_field(',
                          'static void free_nodal_field(',
                          'static void free_equation_field(',
                          'static void free_element_field(',
                          'static void assemble_grad_p_unstripped(',
                          'static void add_element_force_to(']:
            sources.append(function(drive, signature))
        sources += [function(element,'void assemble_div_u(struct'),
                    function(element,'void get_elt_f(E,el')]
        sources.append(drive[drive.index('const char *const eba_power_names'):])
        sources += [function(profile,'static int add_f64('),
                    function(profile,'static int write_eba_power_npz('),
                    (ROOT/'tests/mechanical_power/fixture_main.c').read_text()]
        (cls.path/'test.c').write_text('\n'.join(sources))
        subprocess.run([os.environ.get('MPICC','mpicc'),'-std=gnu99',
                        '-Wno-deprecated-non-prototype','-I'+str(ROOT/'lib'),
                        '-I'+str(ROOT/'tests/mechanical_power'),str(cls.path/'test.c'),
                        str(ROOT/'lib/Npz_writer.c'),'-lz','-lm','-o',str(cls.path/'test')],
                       check=True,text=True)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def run_case(self, mode):
        p = self.path/f'{mode}.npz'
        run = subprocess.run([str(self.path/'test'),str(p),str(mode)],
                             check=False,capture_output=True,text=True,
                             env={**os.environ,'OMPI_MCA_btl':'self'})
        if run.returncode: raise RuntimeError(run.stderr)
        return np.load(p),run.stdout

    def test_closed_budget_and_serialized_snapshot(self):
        a, log = self.run_case(0)
        self.assertEqual(a['mechanical_valid'],1)
        self.assertEqual(a['mechanical_step'],12)
        self.assertEqual(a['mechanical_elapsed_time'],0.125)
        self.assertAlmostEqual(float(a['mechanical_Roperator']),0,places=11)
        self.assertAlmostEqual(float(a['mechanical_Rbody_split']),0,places=11)
        self.assertGreater(abs(float(a['mechanical_Wpressure'])),0)
        self.assertGreater(abs(float(a['mechanical_Pother'])),0)
        self.assertGreater(abs(float(a['mechanical_Wtraction'])),0)
        self.assertAlmostEqual(float(a['mechanical_Rmechanical']),
                               float(a['mechanical_Rheating_operator']),places=11)
        for key in a.files:
            if key.endswith('_shell_integral'):
                self.assertAlmostEqual(float(a[key].sum()),float(a[key[:-15]]),places=11)
        for line in log.splitlines():
            fields=line.split()
            if len(fields)==2 and 'mechanical_'+fields[0] in a:
                self.assertEqual(float(fields[1]),float(a['mechanical_'+fields[0]]))
        # The producer intentionally advances the state before NPZ output.
        self.assertEqual(a['mechanical_pre_rigid_rotation'],1)

    def test_free_dof_residual_is_visible(self):
        a,_=self.run_case(1)
        self.assertGreater(abs(float(a['mechanical_Roperator'])),1e-6)
        self.assertAlmostEqual(float(a['mechanical_Rbody_split']),0,places=11)

    def test_missing_snapshot_is_not_zero_power(self):
        a,_=self.run_case(2)
        self.assertEqual(a.files,['mechanical_valid'])
        self.assertEqual(a['mechanical_valid'],0)

    def test_zero_dissipation_number(self):
        a,_=self.run_case(3)
        self.assertEqual(a['mechanical_scale_Di_over_Atemp'],0)
        for key in a.files:
            if key.startswith('mechanical_R') or key.endswith('_shell_integral'):
                self.assertTrue(np.all(a[key]==0))

    def test_pressure_gauge_preserves_closed_budget(self):
        a,_=self.run_case(4)
        self.assertAlmostEqual(float(a['mechanical_Roperator']),0,places=10)

    def test_stationary_state_has_no_work(self):
        a,_=self.run_case(5)
        for key in a.files:
            if key.startswith(('mechanical_P','mechanical_W','mechanical_Q',
                               'mechanical_Dvisc','mechanical_R')):
                self.assertTrue(np.all(a[key]==0),key)

if __name__ == '__main__':
    unittest.main()
