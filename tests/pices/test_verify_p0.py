"""Negative cases for the output gate, independent synthetic node records."""
import gzip
from pathlib import Path
import tempfile
import unittest

from verify_p0 import read_nodes


class NodeGate(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.path = Path(self.directory.name)/'velo.0.1.gz'
        self.lines = ['1 125 1e-6 300 3700 3400','1 125']
        for i in range(1,126):
            temperature = 3700 if i%5 == 1 else 300 if i%5 == 0 else 2000
            self.lines.append('0 0 0 %d' % temperature)

    def check(self):
        with gzip.open(self.path,'wt') as stream:
            stream.write('\n'.join(self.lines)+'\n')
        return read_nodes(self.path,1)

    def test_valid(self):
        self.assertEqual(len(self.check()),125)

    def test_nan(self):
        self.lines[3]='0 0 0 nan'
        with self.assertRaises(ValueError): self.check()

    def test_nonfinite_velocity(self):
        self.lines[3]='inf 0 0 2000'
        with self.assertRaises(ValueError): self.check()

    def test_truncated(self):
        self.lines.pop()
        with self.assertRaises(ValueError): self.check()

    def test_wrong_time(self):
        self.lines[0]='1 125 0 300 3700 3400'
        with self.assertRaises(ValueError): self.check()

    def test_wrong_bottom(self):
        self.lines[2]='0 0 0 3699'
        with self.assertRaises(ValueError): self.check()

    def test_wrong_top(self):
        self.lines[6]='0 0 0 301'
        with self.assertRaises(ValueError): self.check()


if __name__ == '__main__': unittest.main()
