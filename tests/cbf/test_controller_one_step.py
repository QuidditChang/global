"""Run the production Controller march method with a lightweight fake solver."""
import ast
from pathlib import Path
import unittest
source = Path(__file__).resolve().parents[2] / 'CitcomS/Controller.py'
cls = next(n for n in ast.parse(source.read_text()).body if isinstance(n, ast.ClassDef))
method = next(n for n in cls.body if isinstance(n, ast.FunctionDef) and n.name == 'march')
namespace = {}
exec(compile(ast.Module(body=[method], type_ignores=[]), str(source), 'exec'), namespace)

class Fake:
    march = namespace['march']
    def __init__(self, step):
        self.step = step; self.done = False; self.saved = []; self.advances = 0
    def startTimestep(self): self.step += 1
    def stableTimestep(self): return .1
    def advance(self, dt): self.advances += 1
    def endTimestep(self, totalTime, steps): self.done = self.step >= steps
    def save(self): self.saved.append(self.step)
    def endSimulation(self): self.ended = True

class OneStep(unittest.TestCase):
    def test_restart_one_step(self):
        c = Fake(13600); c.march(steps=13601)
        self.assertEqual(c.advances, 1); self.assertEqual(c.saved, [13601])
        self.assertTrue(c.ended)
    def test_already_done(self):
        c = Fake(13601); c.march(steps=13601)
        self.assertEqual(c.advances, 0)

if __name__ == '__main__': unittest.main()
