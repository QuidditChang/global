"""Build a real MPI initialization-only driver from unmodified solver sources.

Only a generated copy of bin/Citcom.c is changed: it snapshots and exits after
initial_conditions, before the first Stokes solve. GNU ld --wrap records the
production pre-age call. No solver source is rewritten or stubbed.
"""
import argparse
import concurrent.futures
import hashlib
import json
import os
import platform
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--build-dir', required=True, type=Path)
    p.add_argument('--optimization', choices=['O0', 'O1', 'O2'], default='O1')
    p.add_argument('--jobs', type=int, default=4)
    a = p.parse_args()
    build = a.build_dir.resolve()
    build.mkdir(parents=True, exist_ok=True)
    listing = (ROOT / 'lib/Makefile.am').read_text().split('sources =', 1)[1].split('EXTRA_DIST', 1)[0]
    sources = [ROOT / 'lib' / v for v in listing.replace('\\', '').split() if v.endswith('.c')]
    source = (ROOT / 'bin/Citcom.c').read_text()
    anchor = '      initial_conditions(E);'
    assert source.count(anchor) == 1, 'initialization driver anchor changed'
    source = source.replace('int main(argc,argv)',
        'void preage_test_final(struct All_variables *E);\nint main(argc,argv)', 1)
    source = source.replace(anchor, anchor + '\n      preage_test_final(E);\n'
        '      MPI_Finalize();\n      return 0;', 1)
    driver = build / 'PreageInitialization.c'
    driver.write_text(source)
    sources += [ROOT / 'bin/CitcomSFull.c', driver, HERE / 'fixture.c']
    compiler = os.environ.get('MPICC', 'mpicc')
    revision = subprocess.check_output(['git', '-C', str(ROOT), 'rev-parse', 'HEAD'], text=True).strip()
    flags = [compiler, '-std=gnu99', '-w', '-Wno-error=implicit-function-declaration',
             '-Wno-error=implicit-int', '-Wno-error=int-conversion',
             '-Wno-error=incompatible-pointer-types', '-' + a.optimization, '-g', '-DUSE_GZDIR',
             '-DPICES_SOLVER_COMMIT="' + revision + '"', '-I' + str(ROOT / 'lib')]
    def compile_one(src):
        obj = build / (src.stem + '.o')
        extra = []
        # Apple ld has no GNU --wrap. Rename only the implementation and
        # wrapper symbols at compile time; callers still enter the wrapper.
        if platform.system() == 'Darwin':
            if src.name == 'Thermal_preage.c':
                extra = ['-Dthermal_preage_run=__real_thermal_preage_run']
            elif src == HERE / 'fixture.c':
                extra = ['-D__wrap_thermal_preage_run=thermal_preage_run']
        result = subprocess.run(flags + extra + ['-c', str(src), '-o', str(obj)], capture_output=True, text=True)
        if result.returncode:
            raise RuntimeError(str(src) + '\n' + result.stderr)
        return str(obj)
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.jobs) as pool:
        objects = list(pool.map(compile_one, sources))
    exe = build / 'PreageInitialization'
    wrap = [] if platform.system() == 'Darwin' else ['-Wl,--wrap=thermal_preage_run']
    subprocess.run([compiler, *objects, *wrap, '-lz', '-lm', '-o', str(exe)], check=True)
    (build / 'VALIDATION_BUILD.txt').write_text(
        f'source={ROOT}\nrevision={revision}\noptimization={a.optimization}\n'
        'parser_workaround=False\nproduction_sources_modified=False\n'
        'driver=initial_conditions then snapshot then MPI_Finalize\n')
    tracked = sources + list((ROOT / 'lib').glob('*.h'))
    hashes = {str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path):
              hashlib.sha256(path.read_bytes()).hexdigest() for path in tracked}
    (build / 'SOURCE_SHA256.json').write_text(json.dumps(hashes, indent=2, sort_keys=True) + '\n')
    print(exe)


if __name__ == '__main__':
    main()
