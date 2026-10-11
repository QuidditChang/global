"""Exercise the actual C particle TA reaction on positive, sharp fields."""
import ctypes, math, random, subprocess, tempfile, unittest
from pathlib import Path

class Reaction(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        source=(Path(__file__).resolve().parents[2]/'lib/Pices.c').read_text()
        code=source.split('static double ta_particle_relax(',1)[1].split('\nstatic void pices_assimilate',1)[0]
        cls.tmp=tempfile.TemporaryDirectory();p=Path(cls.tmp.name)
        (p/'reaction.c').write_text('double ta_particle_relax('+code)
        subprocess.run(['cc','-shared','-fPIC',str(p/'reaction.c'),'-o',str(p/'reaction.so')],check=True)
        cls.lib=ctypes.CDLL(str(p/'reaction.so'));cls.fn=cls.lib.ta_particle_relax
        cls.fn.argtypes=[ctypes.POINTER(ctypes.c_int)]+[ctypes.POINTER(ctypes.c_double)]*3+[ctypes.c_double]
        cls.fn.restype=ctypes.c_double
    @classmethod
    def tearDownClass(cls):cls.tmp.cleanup()
    def apply(self,old,weights,alpha,target):
        # Test dimensional Kelvin; convex reaction is affine-scale invariant.
        return self.fn((ctypes.c_int*9)(0,*range(1,9)),(ctypes.c_double*9)(0,*weights),(ctypes.c_double*9)(0,*alpha),(ctypes.c_double*9)(0,*target),old)
    def test_cold_particle_must_warm_even_when_grid_cools(self):
        old=9.294574939287713;grid=460.57410993116275;target=400.
        fraction=14.278546137304447/(grid-target)
        self.assertLess(old+fraction*(target-grid),0)
        got=self.apply(old,[.125]*8,[fraction]*8,[target]*8)
        self.assertGreater(got,old);self.assertLess(got,target)
    def test_repeated_exact_relaxation_and_zero_support(self):
        old=10.;tau=.1;dt=.02;a=-math.expm1(-dt/tau)
        actual=old
        for _ in range(100):actual=self.apply(actual,[.125]*8,[a]*8,[300.]*8)
        self.assertAlmostEqual(actual,300.+(old-300.)*math.exp(-100*dt/tau),places=10)
        self.assertEqual(self.apply(old,[.125]*8,[0.]*8,[0.]*8),old)
    def test_variable_support_convex_bounds_and_scaling(self):
        rng=random.Random(63)
        for _ in range(2000):
            weights=[rng.random() for _ in range(8)];total=sum(weights);weights=[w/total for w in weights]
            alpha=[rng.choice([0.,rng.random(),1.]) for _ in range(8)]
            target=[rng.uniform(1,4000) for _ in range(8)];old=rng.uniform(0,4000)
            got=self.apply(old,weights,alpha,target)
            self.assertGreaterEqual(got,min([old]+target)-1e-10)
            self.assertLessEqual(got,max([old]+target)+1e-10)
            nd=self.apply((old-300)/3400,weights,alpha,[(v-300)/3400 for v in target])
            self.assertAlmostEqual(got,300+3400*nd,places=10)

if __name__=='__main__':unittest.main()
