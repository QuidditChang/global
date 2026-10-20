"""Check production physics parity and the actual LSF runtime snapshot.

Usage: python3 tests/pices/check_production_release.py /path/to/runs
No HPC submission, MPI execution or installed runtime modification.
"""
import configparser
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile


def config(path):
    parser = configparser.ConfigParser(interpolation=None)
    parser.read(path)
    return {(section, key): value for section in parser.sections()
            for key, value in parser.items(section)}


def check(runs):
    name = 'cmbhf_EBA_Q0_30_rheol7_scold1.0_LLSVPsV0.04_B24_0.3_PICES'
    production = config(runs / (name + '.cfg'))
    reference = config(runs / (name + '_motion300.cfg'))
    operational = {('CitcomS', 'steps'), ('CitcomS.job', 'name'),
                   ('CitcomS.job', 'walltime'), ('CitcomS.solver', 'datadir'),
                   ('CitcomS.solver', 'datadir_old')}
    assert production.keys() == reference.keys()
    assert {k: v for k, v in production.items() if k not in operational} == {
        k: v for k, v in reference.items() if k not in operational}
    assert production['CitcomS', 'steps'] == '30000'
    assert production['CitcomS.solver.ic', 'restart'] == 'off'
    assert production['CitcomS.solver.tracer', 'tracers_per_element'] == '27'
    assert production['CitcomS.solver.tsolver', 'pices_projection'] == 'lumped'
    assert production['CitcomS.solver.tsolver', 'pices_ta_before_heat'] == 'on'

    lsf = runs / (name + '.lsf')
    subprocess.run(['bash', '-n', str(lsf)], check=True)
    text = lsf.read_text()
    assert '--code-dir "${RUNTIME_INSTALL}"' in text
    cleanup = text.split('set -eo pipefail\n', 1)[1].split('RUNS_DIR=', 1)[0]
    env = dict(os.environ, PICES_MOTION_AUDIT='1', PICES_EMPTY_AUDIT='1',
               PICES_RESTART_SOURCE_COMMIT='old', PICES_REPLAY_CHECKPOINT_ONLY='1',
               PICES_COVERAGE_PROBE='1', PICES_TEST_MODE='rotation')
    cleaned = subprocess.check_output(['bash', '-c', cleanup +
        '\nfor variable in ${!PICES_@}; do printf "%s=%s\\n" "$variable" "${!variable}"; done'],
        env=env, text=True).strip()
    assert cleaned == 'PICES_COVERAGE_STATS=1', cleaned
    snapshot = text.split('<<\'PY\'\n', 1)[1].split('\nPY\n', 1)[0]
    with tempfile.TemporaryDirectory(prefix='pices-release-') as folder:
        root = Path(folder)
        producer, run = root / 'producer', root / 'run'
        site = producer / 'install/lib/python2.6/site-packages'
        component = site / 'CitcomS/Components/Advection_diffusion/Advection_diffusion.py'
        component.parent.mkdir(parents=True)
        component.write_text('version = "original"\n')
        component.with_suffix('.pyc').write_bytes(b'stale bytecode')
        binary = producer / 'install/bin/mpipycitcoms'
        binary.parent.mkdir(parents=True)
        binary.write_bytes(b'original executable')
        egg = producer / 'python/dependency.egg'
        egg.parent.mkdir()
        egg.write_bytes(b'original dependency')
        # Directory links must be copied as independent files, too.
        shared = producer / 'shared'
        shared.mkdir()
        (shared / 'data').write_text('linked data')
        (site / 'linked').symlink_to(shared, target_is_directory=True)
        run.mkdir()
        (run / 'solver_commit.txt').write_text('a' * 40 + '\n')
        subprocess.run([sys.executable, '-c', snapshot, str(producer), str(run)], check=True)
        runtime = run / 'Runtime'
        frozen_component = runtime / component.relative_to(producer)
        frozen_binary = runtime / binary.relative_to(producer)
        frozen_egg = Path((runtime / 'eggs.path').read_text())
        assert not frozen_component.with_suffix('.pyc').exists()
        assert not (runtime / site.relative_to(producer) / 'linked').is_symlink()
        component.write_text('version = "rebuilt"\n')
        binary.write_bytes(b'rebuilt executable')
        egg.write_bytes(b'rebuilt dependency')
        (shared / 'data').write_text('changed')
        assert frozen_component.read_text() == 'version = "original"\n'
        assert frozen_binary.read_bytes() == b'original executable'
        assert frozen_egg.read_bytes() == b'original dependency'
        assert (runtime / site.relative_to(producer) / 'linked/data').read_text() == 'linked data'
        assert json.loads((runtime / 'manifest.json').read_text())['solver_commit'] == 'a' * 40
        # Existing snapshots cannot be overwritten by another launch.
        retry = subprocess.run([sys.executable, '-c', snapshot, str(producer), str(run)],
                               capture_output=True, text=True)
        assert retry.returncode != 0 and 'FileExistsError' in retry.stderr
    return dict(physical_configuration_matches_motion300=True,
                lsf_syntax='PASS', inherited_diagnostics_cleared=True,
                runtime_snapshot_survives_producer_rebuild=True,
                symlink_targets_copied=True, stale_bytecode_excluded=True,
                existing_snapshot_rejected=True)


if __name__ == '__main__':
    print(json.dumps(check(Path(sys.argv[1]).resolve()), indent=2))
