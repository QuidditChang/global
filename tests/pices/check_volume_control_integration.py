"""Exercise all five transport controls with the actual solver and MPI.

The synthetic nine-node fixture tests execution and invariants; it is not the
RH20 production mesh and cannot predict RH20 drift. A second invocation with
--radial-partitions 2 exercises migration across radial MPI partitions.
"""
import argparse
from array import array
import csv
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import struct


CASES = ('A_static', 'B_exact', 'B_dt', 'C_half', 'C_quarter')


def fixture(inputs, root, radial_partitions):
    text = (inputs / 'continuous.cfg').read_text()
    changes = dict(RH='20', nodex='9', nodey='9', nodez='9', mgunitx='8',
                   mgunity='8', mgunitz=str(8 // radial_partitions),
                   nprocz=str(radial_partitions), tracers_per_element='64')
    radii = [.55 + .45 * i / 8 for i in range(9)]
    if radial_partitions == 2:
        # The uniform mesh partitions at .775, outside the control support
        # [.56,.72]. Put the fifth radial node through the vortex instead.
        # The coor reader uses float, so use its realized radii consistently
        # in the synthetic reference profile and crossing checks below.
        radii = [struct.unpack('f', struct.pack('f', r))[0]
                 for r in (.55, .5725, .595, .6175, .64, .73, .82, .91, 1.)]
        changes.update(coor='1', coor_file='control_radial.coor', verbose='on')
    for key, value in changes.items():
        pattern = r'^' + re.escape(key) + r'\s*=.*$'
        if re.search(pattern, text, re.M):
            text = re.sub(pattern, key + '=' + value, text, flags=re.M)
        else:
            text += '\n' + key + '=' + value + '\n'
    root.mkdir(parents=True)
    (root / 'DATA').mkdir()
    (root / 'case.cfg').write_text(text)
    if radial_partitions == 2:
        (root / 'control_radial.coor').write_text('nodez 9\n' + ''.join(
            '%d %.17g\n' % (i + 1, radius) for i, radius in enumerate(radii)))
    # Interpolate the existing manufactured linear reference profile onto the
    # refined fixture mesh. Production launchers never synthesize these files.
    (root / 'refstate_EBA_PICES_P4.txt').write_text(''.join(
        '%.17g 1 %.17g 1 %.17g\n' % (1.1 - .1 * fraction,
                                            .8 - .4 * fraction,
                                            1.2 - .2 * fraction)
        for fraction in ((r - radii[0]) / (radii[-1] - radii[0]) for r in radii)))
    shutil.copytree(inputs / 'pices_p4_forcing', root / 'pices_p4_forcing',
                    ignore=shutil.ignore_patterns('._*'))
    for path in (root / 'pices_p4_forcing').glob('age.*'):
        path.write_text('50\n' * 81)
    return radii


def radial_crossings(root, boundary):
    """Prove ownership crossings using existing production migration logging.

    full_lost_souls logs departing particles as (theta,phi,r). With nprocx=y=1,
    Full_parallel_related.c assigns rank % nprocz as the radial partition.
    Check that initial snapshots obey that ownership, then require departures
    outside the owner's radial shell in both directions. Counts are departure
    events, not unique particles, and include predictor/corrector stages.
    """
    initial_counts = {}
    for path in sorted((root / 'Initial').glob('baseline.*.bin')):
        with path.open('rb') as stream:
            header = struct.unpack('8i', stream.read(8 * struct.calcsize('i')))
            _, rank, ranks, caps, _, _, basic, _ = header
            assert ranks == 24 and caps == 1 and basic >= 3, header
            count, = struct.unpack('i', stream.read(struct.calcsize('i')))
            stream.seek(2 * count * struct.calcsize('d'), 1)
            radii = array('d')
            radii.fromfile(stream, count)
        assert all((r <= boundary + 1e-12 if rank % 2 == 0 else r >= boundary - 1e-12)
                   for r in radii), ('initial radial ownership mismatch', rank)
        initial_counts[str(rank)] = count
    assert len(initial_counts) == 24, 'missing initial snapshot rank'
    event = re.compile(r'tracer#=\d+ xx=\([^,]+,[^,]+,([^\)]+)\)')
    totals = dict(inner_to_outer=0, outer_to_inner=0)
    by_rank = {}
    # Legacy %g prints six significant digits. Exclude ambiguous rounded
    # interface positions rather than counting a non-crossing as migration.
    tolerance = 1e-6
    for path in sorted((root / 'DATA').glob('*/*.tracer_log.*')):
        rank = int(path.name.rsplit('.', 1)[1])
        count = 0
        with path.open() as stream:
            for line in stream:
                found = event.search(line)
                if not found:
                    continue
                radius = float(found.group(1))
                if rank % 2 == 0 and radius > boundary + tolerance:
                    totals['inner_to_outer'] += 1
                    count += 1
                elif rank % 2 == 1 and radius < boundary - tolerance:
                    totals['outer_to_inner'] += 1
                    count += 1
        if count:
            by_rank[str(rank)] = count
    assert all(totals.values()), ('radial MPI migration not exercised in both directions', totals)
    return dict(interface_radius_nd=boundary, departure_events=totals,
                departures_by_rank=by_rank, initial_particles_by_rank=initial_counts,
                evidence='initial snapshot shell ownership plus full_lost_souls departure radii',
                counting_note='events across all moving cases; not unique particles')


def read_csv(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('build', type=Path)
    parser.add_argument('runs', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--radial-partitions', type=int, choices=(1, 2), default=1)
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    root = args.output.resolve()
    radii = fixture(args.runs / 'pices_p8a', root, args.radial_partitions)
    ranks = 12 * args.radial_partitions
    env = {k: v for k, v in os.environ.items() if not k.startswith('PICES_')}
    env.update(PICES_VOLUME_CONTROL='1', PICES_VOLUME_AUDIT='1',
               PICES_CONTROL_BASE_STEPS='12', PICES_CONTROL_END_TIME_ND='4e-6',
               PICES_VOLUME_AUDIT_INTERVAL='12', OMP_NUM_THREADS='1')
    command = [os.environ.get('MPIEXEC', '/usr/local/bin/mpiexec'),
               '--oversubscribe', '-n', str(ranks),
               str(args.build.resolve() / 'CitcomSFull'), 'case.cfg']
    with (root / 'stdout').open('w') as out, (root / 'stderr').open('w') as err:
        result = subprocess.run(command, cwd=root, env=env, stdout=out,
                                stderr=err, timeout=900)
    assert result.returncode == 0, (result.returncode,
                                    (root / 'stderr').read_text()[-6000:])
    assert not list((root / 'DATA').glob('*/*.chkpt.*')), 'unexpected production continuation'
    case_rows = read_csv(root / 'cases.csv')
    assert len(case_rows) == 5, case_rows
    baselines = []
    counts = []
    drift = {}
    max_residual = 0.0
    steps = []
    for name in CASES:
        raw = read_csv(root / 'Cases' / name / 'Diagnostic/VolumeAudit/globals.csv')
        assert all(r['status'] in ('pending', 'accepted') for r in raw)
        rows = [r for r in raw if r['status'] == 'accepted']
        by_step = {}
        for row in rows:
            by_step.setdefault(int(row['step']), []).append(row)
        assert sorted(by_step) == list(range(max(by_step) + 1))
        assert {int(r['step']) for r in raw} == set(by_step), 'pending-only step'
        assert all(len(batch) == 25 for batch in by_step.values())
        first = sorted(by_step[0], key=lambda r: int(r['flavor']))
        last = by_step[max(by_step)]
        shell = float(first[0]['shell_volume_nd'])
        baselines.append([float(r['volume_after_nd']) for r in first])
        counts.append([int(r['count_after']) for r in first])
        steps.append(max(by_step))
        for batch in by_step.values():
            assert abs(sum(float(r['volume_after_nd']) for r in batch) / shell - 1) < 1e-12
            for row in batch:
                flavor = int(row['flavor'])
                assert int(row['count_after']) == counts[-1][flavor]
                assert row['count_before'] == row['count_moved'] == row['count_after']
                assert float(row['source_nd']) == 0
                value = abs(float(row['residual_nd'])) / shell
                max_residual = max(max_residual, value)
                assert value < 1e-10
                assert math.isfinite(float(row['volume_after_nd']))
                if name == 'A_static':
                    assert abs(float(row['volume_after_nd']) - baselines[-1][flavor]) / shell < 1e-12
        c24 = next(r for r in last if int(r['flavor']) == 24)
        assert counts[-1][24] > 0
        drift[name] = (float(c24['volume_after_nd']) / baselines[-1][24] - 1) * 100
    assert all(b == baselines[0] for b in baselines), 'independent starts differ'
    assert all(c == counts[0] for c in counts), 'independent initial labels differ'
    assert steps[0] == steps[1] == steps[2]
    assert steps[3] == 2 * steps[2] and steps[4] == 4 * steps[2]
    if args.report:
        subprocess.run([os.environ.get('PYTHON', 'python3'), str(args.report.resolve()),
                        str(root), '--no-plots'], check=True)
    summary = dict(status='PASS', ranks=ranks, steps=steps,
                   initial_particles=sum(counts[0]), C24_relative_drift_percent=drift,
                   normalized_budget_residual_max=max_residual,
                   checks=['all five independent cases completed',
                           'identical all-flavor volume and count baselines',
                           'static reconstruction stable',
                           'all-flavor count/source and shell-volume closure',
                           'N/2N/4N production transport steps',
                           'no production checkpoint or solver continuation'])
    if args.radial_partitions == 2:
        summary['radial_MPI_crossings'] = radial_crossings(root, radii[4])
        summary['checks'].append('explicit radial MPI ownership crossing in both directions')
    (root / 'integration_results.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
