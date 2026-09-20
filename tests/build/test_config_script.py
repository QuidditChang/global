"""Build-driver regression: replace stale configure and reject unresolved AR."""
from pathlib import Path
import os
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]

class BuildScript(unittest.TestCase):
    def run_driver(self, bad):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);tools=root/'tools';tools.mkdir()
            (root/'config_script').write_text((ROOT/'config_script').read_text())
            (root/'configure').write_text('#!/bin/bash\necho stale >> "$BUILD_TRACE"\nexit 99\n')
            (root/'configure').chmod(0o755)
            (tools/'make').write_text('#!/bin/bash\necho "make $*" >> "$BUILD_TRACE"\n')
            (tools/'autoreconf').write_text('''#!/bin/bash
printf 'autoreconf\\n' >> "$BUILD_TRACE"
cat > configure <<'CONFIGURE'
#!/bin/bash
printf 'configure\\n' >> "$BUILD_TRACE"
mkdir -p lib
if [ "$BAD_ARCHIVER" = 1 ]; then
    printf 'AR = @AR@\\nRANLIB = @RANLIB@\\n' > lib/Makefile
else
    printf 'AR = ar\\nRANLIB = ranlib\\nCC = cc\\nCXX = c++\\n' > lib/Makefile
fi
CONFIGURE
chmod +x configure
''')
            for p in tools.iterdir():p.chmod(0o755)
            env=dict(os.environ,PATH=str(tools)+':/usr/bin:/bin',
                BUILD_TRACE=str(root/'trace'),BAD_ARCHIVER=str(int(bad)),HDF5ROOT=str(root/'hdf5'))
            result=subprocess.run(['/bin/bash','config_script'],cwd=root,env=env,capture_output=True,text=True)
            return result,(root/'trace').read_text()

    def test_regenerates_existing_configure(self):
        result,trace=self.run_driver(False)
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertNotIn('stale',trace)
        self.assertLess(trace.index('autoreconf'),trace.index('\nconfigure'))
        self.assertIn('make -j4',trace)
        self.assertIn('make install SUBDIRS=lib CitcomS etc bin module',trace)
        self.assertNotIn(' -i',trace)

    def test_unresolved_archiver_stops_before_compilation(self):
        result,trace=self.run_driver(True)
        self.assertEqual(result.returncode,2,result.stderr)
        self.assertIn('unresolved tool substitutions',result.stderr)
        self.assertNotIn('make -j4',trace)
        self.assertNotIn('make install',trace)

if __name__=='__main__':unittest.main()
