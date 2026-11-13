"""Run real EBA PICES fixed-volume controls or a coupled short segment.

Requires a fresh standalone full solver build. Inputs are manufactured small
meshes and cannot establish full-resolution RH20 remap feasibility.
"""
import argparse
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import struct

from check_volume_audit_integration import physical_payload


def native_restart_record(root, step, ranks, mode):
    """Authenticate the actual split checkpoint before writing its certificate."""
    hashes, manifests, commits, clocks = [], [], set(), set()
    for rank in range(ranks):
        stem = root / 'DATA' / str(rank) / ('PICES_P4.chkpt.%d.%d' % (rank, step))
        raw = Path(str(stem) + '.pices.json').read_bytes()
        metadata = json.loads(raw)
        for key, value in dict(magic='CITCOMS_EBA_PICES', schema=7, step=step,
                               rank=rank, mpi_size=ranks, phase='accepted', flavors=25).items():
            assert metadata.get(key) == value, (stem, key)
        assert metadata['composition_volume']['mode'] == mode
        assert any(q['name'] == 'volume_weight' for q in metadata['extraq'])
        for suffix, key in (('', 'checkpoint_sha256'), ('.pices.state', 'accepted_velocity_sha256')):
            assert hashlib.sha256(Path(str(stem) + suffix).read_bytes()).hexdigest() == metadata[key], (stem, suffix)
        hashes.append(hashlib.sha256(raw).hexdigest().encode() + b'\0')
        manifests.append(json.loads(Path(str(stem) + '.pices.manifest').read_text()))
        commits.add(metadata['solver_commit'])
        clocks.add((metadata['time'], metadata['dt'], metadata['total_timesteps']))
    digest = hashlib.sha256(b''.join(hashes)).hexdigest()
    expected = dict(magic='CITCOMS_EBA_PICES_COMPLETE', schema=1,
                    mpi_size=ranks, metadata_set_sha256=digest)
    assert len(commits) == len(clocks) == 1 and all(m == expected for m in manifests)
    return dict(step=step, schema=7, mpi_size=ranks, solver_commit=commits.pop(),
                metadata_set_sha256=digest, payload_hash_verification='local_SHA256_and_solver_MPI_preflight',
                operator_transition='none')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('build', type=Path)
    p.add_argument('runs', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--radial-partitions', type=int, choices=(1, 2), default=1)
    p.add_argument('--mode', choices=('shadow', 'active'), default='shadow')
    p.add_argument('--controls', action='store_true')
    p.add_argument('--stop', type=int, default=2)
    p.add_argument('--restart', type=int)
    p.add_argument('--hops', type=int, default=2)
    p.add_argument('--iterations', type=int, default=1000)
    p.add_argument('--control-steps', type=int, default=12)
    p.add_argument('--control-nodes', type=int, choices=(9, 17), default=9)
    p.add_argument('--timeout', type=int, default=3600)
    p.add_argument('--report', type=Path)
    p.add_argument('--compare-checkpoints', type=Path,
                   help='Independent continuous active run; compare physical endpoint payloads on every rank.')
    a = p.parse_args()
    spec = importlib.util.spec_from_file_location('control_fixture', Path(__file__).with_name('check_volume_control_integration.py'))
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    root = a.output.resolve()
    if a.restart is None:
        helper.fixture(a.runs / 'pices_p8a', root, a.radial_partitions)
    text = (root / 'case.cfg').read_text()
    values = dict(pices_volume_mode=a.mode, pices_volume_hops=str(a.hops),
                  pices_volume_iterations=str(a.iterations), pices_volume_tolerance='1e-12',
                  maxstep=str(a.stop), maxtotstep=str(a.stop + 1), checkpointFrequency='1')
    if a.controls:
        # The analytic support begins at r=.56. Keep one fully motionless
        # radial element next to the CMB, as on the refined production mesh.
        # Otherwise this tiny test's FE interpolation spreads velocity into
        # the boundary cushion, invalidating the no-projection control.
        nodes = a.control_nodes
        coarse = (.55, .555, .58, .61, .64, .73, .82, .91, 1.)
        radii = list(coarse) if nodes == 9 else [v for i in range(8)
            for v in (coarse[i], (coarse[i] + coarse[i + 1]) / 2)] + [1.]
        radii = [struct.unpack('f', struct.pack('f', r))[0] for r in radii]
        values.update(nodex=str(nodes), nodey=str(nodes), nodez=str(nodes),
                      mgunitx=str(nodes - 1), mgunity=str(nodes - 1),
                      mgunitz=str((nodes - 1) // a.radial_partitions),
                      coor='1', coor_file='control_radial.coor',
                      tracers_per_element='64')
        (root / 'control_radial.coor').write_text('nodez %d\n' % nodes + ''.join(
            '%d %.17g\n' % (i + 1, radius) for i, radius in enumerate(radii)))
        (root / 'refstate_EBA_PICES_P4.txt').write_text(''.join(
            '%.17g 1 %.17g 1 %.17g\n' % (1.1 - .1 * f, .8 - .4 * f, 1.2 - .2 * f)
            for f in ((r - radii[0]) / (radii[-1] - radii[0]) for r in radii)))
        for path in (root / 'pices_p4_forcing').glob('age.*'):
            path.write_text('50\n' * (nodes * nodes))
    if a.restart is not None:
        record = native_restart_record(root, a.restart, 12 * a.radial_partitions, a.mode)
        old = root / 'Diagnostic'
        assert not (root / 'Diagnostic.part1').exists(), 'This test supports one split restart only'
        shutil.copy2(root / 'case.cfg', root / 'case.cfg.part1')
        old.rename(root / 'Diagnostic.part1')
        (root / 'restart_manifest.json').write_text(json.dumps(record, indent=2) + '\n')
        values.update(restart='on', solution_cycles_init=str(a.restart), zero_elapsed_time='off',
                      datafile_old='PICES_P4', datadir_old='DATA/%RANK')
    for key, value in values.items():
        pattern = r'^' + re.escape(key) + r'\s*=.*$'
        if re.search(pattern, text, re.M):
            text = re.sub(pattern, key + '=' + value, text, flags=re.M)
        else:
            text += '\n' + key + '=' + value + '\n'
    (root / 'case.cfg').write_text(text)
    env = {k: v for k, v in os.environ.items() if not k.startswith('PICES_')}
    env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    if a.mode == 'shadow':
        env.update(PICES_VOLUME_AUDIT='1', PICES_VOLUME_AUDIT_INTERVAL='1')
    if a.controls:
        env.update(PICES_VOLUME_CONTROL='1', PICES_CONTROL_BASE_STEPS=str(a.control_steps),
                   PICES_CONTROL_END_TIME_ND='4e-6', PICES_VOLUME_AUDIT='1',
                   PICES_VOLUME_AUDIT_INTERVAL=str(a.control_steps))
    ranks = 12 * a.radial_partitions
    with (root / 'stdout').open('w') as out, (root / 'stderr').open('w') as err:
        result = subprocess.run([os.environ.get('MPIEXEC', '/usr/local/bin/mpiexec'),
            '--oversubscribe', '-n', str(ranks), str(a.build.resolve() / 'CitcomSFull'),
            'case.cfg'], cwd=root, env=env, stdout=out, stderr=err, timeout=a.timeout)
    assert result.returncode in ((0,) if a.controls else (0, 8)), (
        result.returncode, (root / 'stderr').read_text()[-6000:])
    cases = ('A_static', 'B_exact', 'B_dt', 'C_half', 'C_quarter', 'D_return') if a.controls else ('.',)
    diagnostics = []
    for name in cases:
        directory = root / 'Cases' / name if a.controls else root
        with (directory / 'Diagnostic/VolumeConserve/steps.csv').open() as stream:
            rows = [r for r in csv.DictReader(stream) if r['status'] == 'accepted']
        assert rows
        expected_first = 0 if a.restart is None else a.restart
        assert [int(r['step']) for r in rows] == list(range(expected_first, int(rows[-1]['step']) + 1))
        failures = [(int(r['step']), r['reason']) for r in rows if r['feasible'] != '1']
        diagnostics.append(dict(case=name, accepted_rows=len(rows), failures=failures,
            max_iterations=max(int(r['iterations']) for r in rows)))
        if not a.controls:
            assert int(rows[-1]['step']) == a.stop
            for rank in range(ranks):
                metadata = json.loads((root / 'DATA' / str(rank) / ('PICES_P4.chkpt.%d.%d.pices.json' % (rank, a.stop))).read_text())
                assert metadata['schema'] == 7
                assert any(q['name'] == 'volume_weight' for q in metadata['extraq'])
    compared = 0
    if a.compare_checkpoints:
        assert not a.controls
        for rank in range(ranks):
            for suffix in ('', '.pices.state'):
                relative = Path('DATA') / str(rank) / ('PICES_P4.chkpt.%d.%d%s' % (rank, a.stop, suffix))
                assert physical_payload(root / relative) == physical_payload(a.compare_checkpoints / relative), relative
                compared += 1
    if a.report:
        command = ['python3', str(a.report.resolve()), str(root), '--stage',
                   'controls' if a.controls else 'restart' if a.restart is not None else a.mode,
                   '--expected-stop', str(a.control_steps if a.controls else a.stop), '--no-plots']
        if a.restart is not None:
            assert a.compare_checkpoints, 'Restart report needs an independent continuous reference'
            command += ['--reference', str(a.compare_checkpoints.resolve())]
        subprocess.run(command, check=True)
    summary = dict(status='EXECUTION_COMPLETE', mode=a.mode, ranks=ranks,
                   controls=a.controls, cases=diagnostics, identical_physical_payloads=compared,
                   note='Execution completion does not certify remap feasibility; require report gate.')
    (root / 'integration_results.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
