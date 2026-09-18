"""Exercise build orchestration with command fixtures, not an HPC compiler."""
import importlib.util
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("fix_pyconfig", ROOT/"tools/fix_generated_pyconfig.py")
fix = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fix)


class BuildEntryTests(unittest.TestCase):
    def test_old_generator_preserves_integers_and_expanded_paths(self):
        from lib2to3.refactor import RefactoringTool, get_fixers_from_package
        template=(ROOT/"m4/cit_python.m4").read_text()
        script=template.split('[#!/usr/bin/env python\n',1)[1].split('# end of file]',1)[0]
        script=script.replace(fix.NEW,fix.OLD)
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp)
            p=d/"pyconfig"
            p.write_text(script)
            fix.repair(p)
            first=p.read_text()
            fix.repair(p)
            self.assertEqual(first,p.read_text())
            converted=RefactoringTool(get_fixers_from_package('lib2to3.fixes')).refactor_string(first,'pyconfig')
            p.write_text(str(converted))
            (d/"config.h").write_text('#define HAVE_MPI 1\n')
            (d/"Makefile").write_text('COUNT = 64\nprefix = /tmp/install\nlibdir = $(prefix)/lib\nCOPY = $(COUNT)\n')
            subprocess.run([sys.executable,str(p),'-h',str(d/'config.h'),'-m',str(d/'Makefile'),'-o',str(d/'config.py')],check=True,capture_output=True)
            namespace={}
            exec((d/'config.py').read_text(),namespace)
            self.assertEqual(namespace['makefile']['COUNT'],64)
            self.assertEqual(namespace['makefile']['COPY'],64)
            self.assertEqual(namespace['makefile']['libdir'],'/tmp/install/lib')

    def test_build_success_and_failures_do_not_stamp_partial_install(self):
        for failure in ('none','configure','compile','install','receipt'):
            with self.subTest(failure=failure), tempfile.TemporaryDirectory() as tmp:
                d=Path(tmp)
                (d/'tools').mkdir()
                (d/'commands').mkdir()
                shutil.copy2(ROOT/'config_script',d/'config_script')
                shutil.copy2(ROOT/'tools/fix_generated_pyconfig.py',d/'tools')
                (d/'pyconfig').write_text(fix.OLD)
                (d/'frozen_current_build.json').write_text('stale')
                configure=d/'configure'
                configure.write_text('#!/bin/bash\n[ "$FAILURE" != configure ]\n')
                configure.chmod(0o755)
                make=d/'commands/make'
                make.write_text('''#!/bin/bash
echo "$*" >> calls
case "$1" in
  distclean) exit 0 ;;
  -j4) [ "$FAILURE" != compile ] ;;
  install) [ "$FAILURE" != install ] || exit 9; touch installed ;;
  *) exit 8 ;;
esac
''')
                make.chmod(0o755)
                (d/'tools/strict_ala_frozen_current.py').write_text('''import os, pathlib
assert pathlib.Path('installed').exists()
assert os.environ['FAILURE'] != 'receipt'
pathlib.Path('frozen_current_build.json').write_text('fresh')
''')
                env=dict(os.environ,FAILURE=failure,HDF5ROOT='/fixture/hdf5',
                         STRICT_ALA_PYTHON=sys.executable,PATH=str(d/'commands')+':'+os.environ['PATH'])
                result=subprocess.run(['bash',str(d/'config_script')],env=env,text=True,capture_output=True)
                receipt=d/'frozen_current_build.json'
                if failure=='none':
                    self.assertEqual(result.returncode,0,result.stderr)
                    self.assertEqual(receipt.read_text(),'fresh')
                    self.assertIn('BUILD_COMPLETE',result.stdout)
                    self.assertIn('mkdir_p=/bin/mkdir -p',(d/'calls').read_text())
                else:
                    self.assertNotEqual(result.returncode,0)
                    self.assertFalse(receipt.exists())
                    self.assertIn('phase='+failure,result.stderr)


if __name__ == '__main__':
    unittest.main()
