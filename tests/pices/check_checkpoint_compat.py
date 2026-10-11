"""Compile the real metadata comparator; only the reviewed commit may differ."""
import ctypes
import pathlib
import subprocess
import tempfile
import unittest

SOURCE = pathlib.Path(__file__).resolve().parents[2] / 'lib/Pices_checkpoint.c'
OLD = b'22399749866f452b3c06960e3d00eea5333a37f0'

class Compatibility(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        root = pathlib.Path(cls.tmp.name)
        code = SOURCE.read_text().split('static int compatible_metadata(', 1)[1].split('\nvoid pices_checkpoint_check_state', 1)[0]
        (root/'check.c').write_text('#include <string.h>\n#include <stddef.h>\nint compatible_metadata(' + code)
        subprocess.run(['cc', '-shared', '-fPIC', str(root/'check.c'), '-o', str(root/'check.so')], check=True)
        cls.lib = ctypes.CDLL(str(root/'check.so'))
        cls.compare = cls.lib.compatible_metadata
        cls.compare.argtypes = [ctypes.c_char_p]*3
        cls.compare.restype = ctypes.c_int
        cls.actual = b'{"solver_commit":"'+OLD+b'","physics":"abc","time":1,"checkpoint_sha256":"xyz"}'
        cls.expected = cls.actual.replace(OLD, b'5'*40)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_opt_in_only(self):
        self.assertEqual(self.compare(self.actual, self.actual, None), 1)
        self.assertEqual(self.compare(self.actual, self.expected, None), 0)
        self.assertEqual(self.compare(self.actual, self.expected, OLD), 2)
        self.assertEqual(self.compare(self.actual, self.expected, b'5'*40), 0)
        self.assertEqual(self.compare(self.actual.replace(OLD,b'4'*40), self.expected, OLD), 0)

    def test_every_other_byte_remains_strict(self):
        for field, change in [(b'abc',b'abd'), (b'xyz',b'xy0'), (b'"time":1',b'"time":2')]:
            self.assertEqual(self.compare(self.actual,self.expected.replace(field,change),OLD),0)
        self.assertEqual(self.compare(self.actual+b' ',self.expected,OLD),0)
        self.assertEqual(self.compare(b'{"solver_commit":"short"}',self.expected,OLD),0)

if __name__ == '__main__':
    unittest.main()
