"""Exercise production particle-depth TA with its real column/HSC helpers."""
import ctypes
import math
from pathlib import Path
import random
import re
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]


def function(source, signature):
    start = source.rindex(signature)
    opening = source.index('{', start)
    masked = re.sub(r'/\*.*?\*/|//[^\n]*|"(?:\\.|[^"\\])*"',
                    lambda match: ' ' * len(match.group()), source, flags=re.S)
    end, depth = opening + 1, 1
    while depth:
        depth += (masked[end] == '{') - (masked[end] == '}')
        end += 1
    return source[start:end]


class ParticleDepthReaction(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        lith = (ROOT / 'lib/Lith_age.c').read_text()
        pices = (ROOT / 'lib/Pices.c').read_text()
        signatures = [
            'static float effective_plate_age_nd(',
            'static double lith_age_surface_anomaly(',
            'static double lith_age_target_temperature(',
            'static double lith_age_old_temperature_weight(',
            'static double lith_age_asml_thickness(',
            'static double lith_age_relaxation_fraction(',
            'void lith_age_relaxation_at_depth(',
        ]
        code = r'''
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "lith_age_weight.h"
void pices_fail(struct All_variables *E, const char *message) { abort(); }
'''
        code += '\n'.join(function(lith, signature) for signature in signatures)
        code += '\n' + function(pices, 'int pices_temperature_valid(')
        code += '\n' + function(pices, 'static void ta_particle_depth_reaction(')
        code += r'''
double evaluate(double depth, double dt_Ma, double lower_ref, double upper_ref,
                double radial, const double *columns, const double *flags,
                const double *ages_Ma, double old, double *fraction, double *source)
{
    struct All_variables E;
    int column,nodes[9]={0,1,3,5,7,2,4,6,8};
    double weights[9]={0}, reference[3]={0,lower_ref,upper_ref};
    float signed_distance[5]={0}, ages[5]={0};
    memset(&E,0,sizeof(E));
    E.control.lith_age_asml=1;
    E.control.lith_age_asml_profile=1;
    E.control.lith_age_asml_tau_Ma=.25;
    E.control.lith_age_asml_logistic_center=80./127.42;
    E.control.lith_age_asml_logistic_width=4./127.42;
    E.control.lith_age_depth=.02;
    E.control.max_plate_age_Ma=70;
    E.data.scalet=1.e6;
    E.data.Ttop=300;E.data.ref_temperature=3400;
    E.advection.timestep=dt_Ma/E.data.scalet;
    E.mesh.nox=E.mesh.noy=2;E.mesh.noz=2;
    E.lmesh.nox=E.lmesh.noy=2;E.lmesh.noz=2;
    E.lmesh.nxs=E.lmesh.nys=E.lmesh.nzs=1;
    E.refstate.Tref=reference;E.refstate.temperature_surface=.4;
    E.flag_depth2=signed_distance;E.age_t=ages;
    for(column=0;column<4;column++) {
        signed_distance[column+1]=flags[column];
        ages[column+1]=ages_Ma[column]/E.data.scalet;
        weights[column+1]=columns[column]*(1.-radial);
        weights[column+5]=columns[column]*radial;
    }
    ta_particle_depth_reaction(&E,1,nodes,weights,depth,fraction,source);
    return (1.-*fraction)*old+*source;
}
'''
        cls.temp = tempfile.TemporaryDirectory(prefix='pices-logistic-ta-')
        directory = Path(cls.temp.name)
        source = directory / 'reaction.c'
        source.write_text(code)
        compiler = shutil.which('mpicc')
        if not compiler:
            raise unittest.SkipTest('mpicc is required for the production structs')
        subprocess.run([compiler, '-shared', '-fPIC', '-std=gnu99',
                        '-I' + str(ROOT / 'lib'), str(source), '-lm',
                        '-o', str(directory / 'reaction.so')], check=True)
        cls.lib = ctypes.CDLL(str(directory / 'reaction.so'))
        cls.evaluate = cls.lib.evaluate
        cls.evaluate.argtypes = [ctypes.c_double] * 5 + [ctypes.POINTER(ctypes.c_double)] * 3 + [ctypes.c_double] + [ctypes.POINTER(ctypes.c_double)] * 2
        cls.evaluate.restype = ctypes.c_double
        cls.h = ctypes.c_float(.02).value

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def apply(self, depth, dt=.05, lower=.6, upper=.4, radial=.3,
              columns=(.2, .3, .5, 0.), flags=(.2,) * 4,
              ages=(70.,) * 4, old=.2):
        array = ctypes.c_double * 4
        fraction, source = ctypes.c_double(), ctypes.c_double()
        after = self.evaluate(depth, dt, lower, upper, radial,
                              array(*columns), array(*flags), array(*ages), old,
                              ctypes.byref(fraction), ctypes.byref(source))
        return after, fraction.value, source.value

    def test_exact_particle_cutoff_and_no_surface_reaction(self):
        for depth in (-1.e-9, 0., self.h, self.h + 1.e-9):
            self.assertEqual(self.apply(depth), (.2, 0., 0.))
        small = self.apply(self.h * (1. - 1.e-5))[1]
        smaller = self.apply(self.h * (1. - 5.e-6))[1]
        self.assertGreater(small, 0.)
        self.assertAlmostEqual(smaller / small, .25, delta=1.e-4)

    def test_rate_uses_particle_depth_and_target_uses_radial_reference(self):
        depth = self.h / 2
        first = self.apply(depth, radial=.1)
        second = self.apply(depth, radial=.9)
        self.assertAlmostEqual(first[1], second[1], places=15)
        age_nd = ctypes.c_float(70. / 1.e6).value
        anomaly = .4 * math.erfc(depth / (2 * math.sqrt(age_nd)))
        self.assertAlmostEqual(first[2] / first[1], .58 - anomaly, places=14)
        self.assertAlmostEqual(second[2] / second[1], .42 - anomaly, places=14)

    def test_trench_columns_keep_their_own_cutoff(self):
        depth = self.h / 2
        _, full_fraction, full_source = self.apply(depth, columns=(0., 1., 0., 0.))
        _, mixed_fraction, mixed_source = self.apply(
            depth, columns=(.75, .25, 0., 0.), flags=(0., .2, 0., 0.))
        self.assertAlmostEqual(mixed_fraction, .25 * full_fraction, places=15)
        self.assertAlmostEqual(mixed_source, .25 * full_source, places=15)
        self.assertEqual(self.apply(.003 * self.h, flags=(0.,) * 4), (.2, 0., 0.))
        # A sub-grid trench layer still contains valid particles even when
        # there is no eligible interior radial node.
        self.assertGreater(self.apply(.001 * self.h, flags=(0.,) * 4)[1], 0.)

    def test_constant_target_time_split_and_convexity(self):
        depth = self.h * .4
        whole = self.apply(depth, dt=.1, lower=.45, upper=.45, ages=(0.,) * 4)[0]
        split = .2
        for _ in range(10):
            split = self.apply(depth, dt=.01, lower=.45, upper=.45,
                               ages=(0.,) * 4, old=split)[0]
        self.assertAlmostEqual(split, whole, delta=2.e-9)
        randomizer = random.Random(2026)
        for _ in range(200):
            columns = [randomizer.random() for _ in range(3)] + [0.]
            total = sum(columns)
            columns = [weight / total for weight in columns]
            old = randomizer.uniform(0., 1.)
            after, fraction, _ = self.apply(
                randomizer.uniform(0., self.h * 1.1), old=old,
                lower=.45, upper=.45, ages=(0.,) * 4, columns=columns,
                flags=tuple(randomizer.choice((0., .08, .2)) for _ in range(4)))
            self.assertGreaterEqual(fraction, 0.)
            self.assertLessEqual(fraction, 1. + 1.e-14)
            self.assertGreaterEqual(after, min(old, .45) - 1.e-14)
            self.assertLessEqual(after, max(old, .45) + 1.e-14)


if __name__ == '__main__':
    unittest.main()
