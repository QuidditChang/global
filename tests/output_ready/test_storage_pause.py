"""Run the production safe-point method with four concurrent barrier ranks."""
import ast
import json
import os
from pathlib import Path
import tempfile
import threading
import time
from types import SimpleNamespace
import unittest
from unittest.mock import patch

SOURCE = Path(__file__).resolve().parents[2] / 'CitcomS/Controller.py'
tree = ast.parse(SOURCE.read_text())
cls = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == 'Controller')
ns = {}
exec(compile(ast.Module(body=[n for n in cls.body if isinstance(n, ast.FunctionDef)
                             and n.name == 'storagePause'], type_ignores=[]), str(SOURCE), 'exec'), ns)


class PauseTest(unittest.TestCase):
    def test_rank_zero_ack_means_all_writers_at_barrier(self):
        with tempfile.TemporaryDirectory() as td, patch.dict(os.environ, {'CITCOMS_STORAGE_CONTROL': td}):
            root = Path(td); request = root / 'request.json'
            request.write_text(json.dumps({'token': 'migration-1'}))
            barrier = threading.Barrier(4, timeout=5)
            written, resumed, errors = [], [], []
            def rank(i):
                try:
                    time.sleep(.02 * i)
                    written.append(i)
                    c = SimpleNamespace(rank=i, barrier=barrier.wait)
                    ns['storagePause'](SimpleNamespace(solver=SimpleNamespace(communicator=c), step=50))
                    resumed.append(i)
                except Exception as e: errors.append(e)
            threads = [threading.Thread(target=rank, args=(i,), daemon=True) for i in range(4)]
            for t in threads: t.start()
            deadline = time.monotonic()+3
            while not (root/'solver.paused.json').exists() and time.monotonic()<deadline: time.sleep(.01)
            self.assertEqual(sorted(written), list(range(4)))
            self.assertEqual(resumed, [])
            self.assertEqual(json.loads((root/'solver.paused.json').read_text()), {'token':'migration-1','step':50})
            request.unlink()
            for t in threads: t.join(3); self.assertFalse(t.is_alive())
            self.assertFalse(errors)
            self.assertEqual(sorted(resumed), list(range(4)))

    def test_no_request_still_enters_both_collectives(self):
        with tempfile.TemporaryDirectory() as td, patch.dict(os.environ, {'CITCOMS_STORAGE_CONTROL': td}):
            calls=[]
            communicator=SimpleNamespace(rank=0,barrier=lambda:calls.append('barrier'))
            ns['storagePause'](SimpleNamespace(solver=SimpleNamespace(communicator=communicator),step=1))
            self.assertEqual(calls,['barrier','barrier'])


if __name__ == '__main__': unittest.main()
