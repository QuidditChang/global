"""Compare disabled-control production output to the committed parent solver.

Use a fresh full build made by tests/cbf/build_validation.py. Reuse unchanged
objects and rebuild the four touched production translation units from a
specified pre-control commit;
omit the new opt-in control object. This avoids creating another worktree or
altering any source file. The reference must precede the opt-in control hooks.
"""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
OLD_UNITS = ('Pices', 'Pices_volume_audit', 'Tracer_setup', 'Full_tracer_advection')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('build', type=Path)
    parser.add_argument('runs', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--baseline-ref', required=True)
    args = parser.parse_args()
    build, output = args.build.resolve(), args.output.resolve()
    parent = output / 'parent_build'
    parent.mkdir(parents=True)
    compiler = os.environ.get('MPICC', '/usr/local/bin/mpicc')
    commit = subprocess.check_output(['git', '-C', str(ROOT), 'rev-parse', args.baseline_ref], text=True).strip()
    changed = set(subprocess.check_output(['git', '-C', str(ROOT), 'diff', '--name-only',
                                          commit, '--', 'lib', 'bin'], text=True).splitlines())
    allowed = {'lib/' + name + '.c' for name in OLD_UNITS}
    allowed.update(('lib/Makefile.am', 'lib/pices_volume_audit.h',
                    'lib/Pices_volume_control.c', 'lib/pices_volume_control.h',
                    'lib/pices_volume_control_math.h'))
    if changed - allowed:
        raise ValueError('Reference differs in other solver units; build it independently: ' +
                         ', '.join(sorted(changed - allowed)))
    flags = [compiler, '-std=gnu99', '-w', '-Wno-error=implicit-function-declaration',
             '-Wno-error=implicit-int', '-Wno-error=int-conversion', '-O1', '-g',
             '-DUSE_GZDIR', '-DPICES_SOLVER_COMMIT="' + commit + '"',
             '-I' + str(ROOT / 'lib'), '-I' + str(ROOT / 'tests/cbf')]
    for name in OLD_UNITS:
        source = parent / (name + '.c')
        source.write_bytes(subprocess.check_output(
            ['git', '-C', str(ROOT), 'show', commit + ':lib/' + name + '.c']))
        subprocess.run(flags + ['-c', str(source), '-o', str(parent / (name + '.o'))], check=True)
    objects = []
    for obj in sorted(build.glob('*.o')):
        if obj.stem in ('Pices_volume_control', 'Citcom_manufactured'):
            continue
        objects.append(str(parent / obj.name if obj.stem in OLD_UNITS else obj))
    subprocess.run([compiler] + objects + ['-lz', '-lm', '-o', str(parent / 'CitcomSFull')], check=True)
    integration = ROOT / 'tests/pices/check_volume_audit_integration.py'
    for name, executable_dir in (('before', parent), ('after', build)):
        subprocess.run([sys.executable, str(integration), str(executable_dir),
                        str(args.runs.resolve()), str(output / name)], check=True)
    spec = importlib.util.spec_from_file_location('audit_integration', integration)
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    compared = 0
    for audit in ('off', 'on'):
        for rank in range(12):
            for step in (0, 2):
                for suffix in ('', '.pices.state'):
                    rel = Path(audit) / 'DATA' / str(rank) / ('PICES_P4.chkpt.%d.%d%s' % (rank, step, suffix))
                    assert helper.physical_payload(output / 'before' / rel) == helper.physical_payload(output / 'after' / rel), rel
                    compared += 1
    result = dict(status='PASS', parent_commit=commit,
                  identical_physical_checkpoint_payloads=compared,
                  checks=['control disabled: production output unchanged from committed solver',
                          'phase-one audit off/on does not change production output'])
    (output / 'results.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
