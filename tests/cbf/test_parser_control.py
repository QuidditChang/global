"""Exercise the production default/minimum/maximum parser on constant strings."""
import ctypes as C
from pathlib import Path
import subprocess
import tempfile
import unittest
from test_cbf_kernel import function, ROOT

class ParserControl(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp=tempfile.TemporaryDirectory()
        code=function((ROOT/'lib/Parsing.c').read_text(),'int interpret_control_string(interpret,essential,Default,minvalue,maxvalue)')
        source=Path(cls.tmp.name)/'parser.c'
        source.write_text('#include <stdio.h>\n#include <string.h>\n#define STRANGE_NUM -98765.4321\n'+code)
        lib=Path(cls.tmp.name)/'parser.so'
        subprocess.run(['cc','-std=gnu99','-Wno-deprecated-non-prototype','-shared','-fPIC',str(source),'-o',str(lib)],check=True)
        cls.lib=C.CDLL(str(lib));cls.parse=cls.lib.interpret_control_string
        cls.parse.argtypes=[C.c_char_p,C.POINTER(C.c_int)]+[C.POINTER(C.c_double)]*3
    @classmethod
    def tearDownClass(cls):cls.tmp.cleanup()
    def test_defaults_and_bounds(self):
        missing=-98765.4321
        for text,expected in [('off',(0,0,missing,missing)),('on',(0,1,missing,missing)),
                              ('2,1,nomax',(0,2,1,missing)),('4,nomin,8',(0,4,missing,8)),
                              ('3,-2,5',(0,3,-2,5)),('nodefault,0,10',(0,missing,0,10)),
                              ('essential',(1,missing,missing,missing)),('1,0',(0,1,0,missing))]:
            with self.subTest(text=text):
                essential=C.c_int();values=[C.c_double() for _ in range(3)]
                self.parse(text.encode(),C.byref(essential),*(C.byref(v) for v in values))
                self.assertEqual((essential.value,*(v.value for v in values)),expected)
if __name__=='__main__':unittest.main()
