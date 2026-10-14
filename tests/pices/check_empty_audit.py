"""Cross-rank/seam topology must survive rounding-bin boundaries."""
import unittest
from audit_empty_run import topology,tv

def cube(x,epsilon=0):
 return {'xyz':[v for p in ((x,0,0),(x+1,0,0),(x+1,1,0),(x,1,0),(x,0,1),(x+1,0,1),(x+1,1,1),(x,1,1)) for v in (p[0]+epsilon,p[1],p[2])], 'volume_nd':1}
class Audit(unittest.TestCase):
 def test_seam_and_missing_halo(self):
  mesh={(0,1):cube(0),(1,2):cube(1,2e-12),(2,3):cube(4)}
  meta={i:dict(ri=-10,ro=10) for i in range(3)}
  adj,missing=topology(mesh,meta,set(mesh))
  self.assertEqual(adj[0,1],{(1,2)})
  self.assertEqual(adj[1,2],{(0,1)})
  self.assertFalse(adj[2,3]);self.assertEqual(len(missing),16)
 def test_background_and_recovery_contrast(self):
  self.assertAlmostEqual(tv([.8,.2,0],[.1,.2,.7]),.7)
  self.assertEqual(tv([0,1],[0,1]),0)
if __name__=='__main__':unittest.main()
