"""Compile the Python binding sources against the actual solver/Python headers.

Complements the standalone numerical build, which does not compile module/.
This does not link Python or emulate the HPC Intel/HDF5 installation.
"""
import argparse
import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--python-include', type=Path, required=True)
    parser.add_argument('--source', action='append', help='specific module source; default: all')
    parser.add_argument('--pythia-include', type=Path, help='directory containing mpi/pympi.h')
    args = parser.parse_args()
    if not (args.python_include / 'Python.h').is_file():
        parser.error('--python-include must contain Python.h')
    makefile = (ROOT / 'module/Makefile.am').read_text()
    sources = re.search(r'^sources\s*=([\s\S]*?)(?:\n\n|\Z)', makefile, re.M)[1]
    sources = [name for name in sources.replace('\\', '').split() if name.endswith('.c')]
    assert sources, 'empty binding source list'
    if args.source:
        if not set(args.source)<=set(sources):
            parser.error('unknown module source')
        sources=args.source
    command = shlex.split(os.environ.get('MPICC', 'mpicc')) + [
        '-std=gnu11', '-Wno-error=implicit-function-declaration',
        '-Wno-error=int-conversion', '-I' + str(ROOT / 'lib'),
        '-I' + str(args.python_include.resolve())]
    if args.pythia_include:
        command += ['-I' + str(args.pythia_include.resolve())]
    with tempfile.TemporaryDirectory(prefix='pices-python-bindings-') as tmp:
        for name in sources:
            result = subprocess.run(command + ['-c', str(ROOT / 'module' / name),
                '-o', str(Path(tmp) / (Path(name).stem + '.o'))],
                capture_output=True, text=True)
            if result.returncode:
                raise SystemExit(name + '\n' + result.stdout + result.stderr)
    print('PASS: compiled %d Python binding sources' % len(sources))


if __name__ == '__main__':
    main()
