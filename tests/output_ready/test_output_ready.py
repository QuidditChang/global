"""Exercise production Controller methods without the legacy Python/MPI runtime."""
import ast
import json
import os
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

SOURCE = Path(__file__).resolve().parents[2] / 'CitcomS/Controller.py'
tree = ast.parse(SOURCE.read_text())
controller = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == 'Controller')
ns = {}
exec(compile(ast.Module(body=[n for n in controller.body if isinstance(n, ast.FunctionDef)
                            and n.name in ('save', 'publishOutputReady')], type_ignores=[]),
             str(SOURCE), 'exec'), ns)


class OutputReadyTest(unittest.TestCase):
    def make_controller(self, directory, step=50, rank=0):
        events = []
        def barrier():
            self.assertFalse((directory / ('%d.json' % step)).exists())
            events.append('barrier')
        solver = SimpleNamespace(communicator=SimpleNamespace(rank=rank, barrier=barrier))
        for method in ('save', 'save_profiles', 'checkpoint', 'save_q_CBF'):
            setattr(solver, method, lambda frequency, name=method: events.append(name))
        obj = SimpleNamespace(solver=solver, step=step, inventory=SimpleNamespace(
            monitoringFrequency=50, profileMonitoringFrequency=50,
            checkpointFrequency=50, monitoringFrequency_CBF=1))
        obj.publishOutputReady = lambda freq: ns['publishOutputReady'](obj, freq)
        return obj, events

    def test_publish_after_all_writers_and_collective(self):
        with tempfile.TemporaryDirectory() as temporary, patch.dict(os.environ,
                {'CITCOMS_OUTPUT_READY_DIR': temporary}):
            directory = Path(temporary)
            obj, events = self.make_controller(directory)
            ns['save'](obj)
            self.assertEqual(events, ['save', 'save_profiles', 'checkpoint', 'save_q_CBF', 'barrier'])
            self.assertEqual(json.loads((directory / '50.json').read_text()),
                {'version': 1, 'step': 50, 'streams': ['caps', 'profiles', 'cbf']})
            self.assertFalse((directory / '.50.tmp').exists())

    def test_cbf_only_step_and_nonroot(self):
        with tempfile.TemporaryDirectory() as temporary, patch.dict(os.environ,
                {'CITCOMS_OUTPUT_READY_DIR': temporary}):
            directory = Path(temporary)
            obj, events = self.make_controller(directory, step=51)
            ns['save'](obj)
            self.assertEqual(json.loads((directory / '51.json').read_text())['streams'], ['cbf'])
            obj, events = self.make_controller(directory, step=52, rank=1)
            ns['save'](obj)
            self.assertIn('barrier', events)
            self.assertFalse((directory / '52.json').exists())

    def test_failed_writer_never_publishes(self):
        with tempfile.TemporaryDirectory() as temporary, patch.dict(os.environ,
                {'CITCOMS_OUTPUT_READY_DIR': temporary}):
            directory = Path(temporary)
            obj, events = self.make_controller(directory)
            def fail(frequency):
                raise IOError('write failed')
            obj.solver.save_q_CBF = fail
            with self.assertRaises(IOError):
                ns['save'](obj)
            self.assertNotIn('barrier', events)
            self.assertEqual(list(directory.iterdir()), [])

    def test_opt_in_only(self):
        with tempfile.TemporaryDirectory() as temporary, patch.dict(os.environ, {}, clear=True):
            obj, events = self.make_controller(Path(temporary))
            ns['save'](obj)
            self.assertNotIn('barrier', events)


if __name__ == '__main__':
    unittest.main()
