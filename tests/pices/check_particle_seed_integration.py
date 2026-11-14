"""Actual MPI initialization/restart regression for cell/volume-uniform tracers.

Uses the manufactured pices_p8a inputs, never RH20 production data. A frozen
pre-change standalone build may be supplied for a legacy bitwise regression.
"""
import argparse
from array import array
from collections import Counter
import json
import math
import os
from pathlib import Path
import re
import shutil
import struct
import subprocess

from check_volume_audit_integration import physical_payload
from check_volume_control_integration import fixture


def settings(text, **values):
    for key, value in values.items():
        pattern = r'^' + re.escape(key) + r'\s*=.*$'
        if re.search(pattern, text, re.M):
            text = re.sub(pattern, key + '=' + str(value), text, flags=re.M)
        else:
            text += '\n' + key + '=' + str(value) + '\n'
    return text


def particles(path):
    data = path.read_bytes()
    nx, ny, nz, *_ = struct.unpack_from('8i', data)
    nn, ne = nx * ny * nz, (nx - 1) * (ny - 1) * (nz - 1)
    off = 44 + 16 + 16 * (nn + 1) + 16 + 8 + 8 * (ne + 1) + 8 * 3 * nn + 16
    basic, extra, flavors, last, count = struct.unpack_from('5i', data, off)
    assert basic == 12 and extra == 2 and flavors == 25, (basic, extra, flavors)
    off += 20
    values = []
    for _ in range(6 + extra):
        row = array('d')
        row.frombytes(data[off + 8:off + 8 * (count + 1)])
        values.append(row)
        off += 8 * (count + 1)
    hosts = array('i')
    hosts.frombytes(data[off + 4:off + 4 * (count + 1)])
    return values, hosts, count, nz - 1


def compare(left, right, ranks, steps):
    checked = 0
    for rank in range(ranks):
        for step in steps:
            for suffix in ('', '.pices.state'):
                relative = Path('DATA') / str(rank) / ('PICES_P4.chkpt.%d.%d%s' % (rank, step, suffix))
                assert physical_payload(left / relative) == physical_payload(right / relative), relative
                checked += 1
    return checked


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('build', type=Path)
    parser.add_argument('runs', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--radial-partitions', type=int, choices=(1, 2), default=1)
    parser.add_argument('--legacy-build', type=Path)
    parser.add_argument('--mode', choices=('volume_uniform', 'cell_uniform'), default='volume_uniform')
    parser.add_argument('--previous-build', type=Path,
                        help='Frozen C-only build: require C step 0/1 payloads and fingerprints unchanged.')
    parser.add_argument('--skip-legacy', action='store_true',
                        help='Skip the already-established A compatibility checks when testing B.')
    parser.add_argument('--nodes', type=int, choices=(5, 9), default=5)
    args = parser.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True)
    ranks = 12 * args.radial_partitions
    env = {k: v for k, v in os.environ.items() if not k.startswith('PICES_')}
    env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')

    def prepare(name, mode=None, stop=2, density=27):
        dest = root / name
        radii = fixture(args.runs / 'pices_p8a', dest, args.radial_partitions)
        if args.nodes == 5:
            radii = (.55, .6625, .775, .8875, 1.) if args.radial_partitions == 1 else (.55, .595, .64, .82, 1.)
            radii = [struct.unpack('f', struct.pack('f', r))[0] for r in radii]
            (dest / 'control_radial.coor').write_text('nodez 5\n' + ''.join(
                '%d %.17g\n' % (i+1, r) for i, r in enumerate(radii)))
            (dest / 'refstate_EBA_PICES_P4.txt').write_text(''.join(
                '%.17g 1 %.17g 1 %.17g\n' % (1.1-.1*f, .8-.4*f, 1.2-.2*f)
                for f in ((r-radii[0])/(radii[-1]-radii[0]) for r in radii)))
            for path in (dest / 'pices_p4_forcing').glob('age.*'):
                path.write_text('50\n' * 25)
        values = dict(maxstep=stop, maxtotstep=stop + 1, checkpointFrequency=1,
                      tracers_per_element=density, pices_projection='lumped',
                      pices_volume_mode='off', verbose='off')
        if args.nodes == 5:
            values.update(nodex=5, nodey=5, nodez=5, mgunitx=4, mgunity=4,
                          mgunitz=4//args.radial_partitions, coor=1, coor_file='control_radial.coor')
        if mode is not None:
            values.update(tracer_seed_distribution=mode, tracer_seed=20261114)
        (dest / 'case.cfg').write_text(settings((dest / 'case.cfg').read_text(), **values))
        return dest, radii

    def execute(dest, binary=None, tag='', failure=None):
        with (dest / ('stdout' + tag)).open('w') as out, (dest / ('stderr' + tag)).open('w') as err:
            result = subprocess.run([os.environ.get('MPIEXEC', '/usr/local/bin/mpiexec'),
                '--oversubscribe', '-n', str(ranks), str((binary or args.build).resolve() / 'CitcomSFull'),
                'case.cfg'], cwd=dest, env=env, stdout=out, stderr=err, timeout=600)
        errors = (dest / ('stderr' + tag)).read_text()
        if failure:
            assert result.returncode not in (0, 8) and failure in errors, (result.returncode, errors[-5000:])
        else:
            assert result.returncode in (0, 8), (result.returncode, errors[-5000:])
        return errors

    comparisons = 0
    if not args.skip_legacy:
        legacy, _ = prepare('legacy_default')
        execute(legacy)
        explicit, _ = prepare('legacy_explicit', 'legacy')
        execute(explicit)
        comparisons += compare(legacy, explicit, ranks, (0, 2))
        if args.legacy_build:
            archived, _ = prepare('legacy_archived')
            execute(archived, args.legacy_build)
            comparisons += compare(legacy, archived, ranks, (0, 2))

    continuous, radii = prepare('uniform_continuous', args.mode)
    errors = execute(continuous)
    match = re.search(r'^TRACER_SEED_GLOBAL (.*)$', errors, re.M)
    assert match, errors[-5000:]
    report = dict(item.split('=', 1) for item in match.group(1).split())
    assert report['mode'] == args.mode
    total = 12 * (args.nodes-1)**3 * 27
    assert int(report['target']) == int(report['generated']) == total
    assert int(report['empty_cells']) == 0 and int(report['min_cell']) > 0
    assert float(report['quota_error_max']) < 1.000001
    if args.mode == 'cell_uniform':
        assert int(report['min_cell']) == int(report['max_cell']) == 27
        assert float(report['lambda_min']) == float(report['lambda_max']) == 27
        assert float(report['quota_error_max']) == 0
        assert report['quota'] == 'constant_per_cell_v1' and report['expectation'] == 'tracers_per_element'
    assert abs(float(report['host_volume']) - 4 * math.pi / 3 * (radii[-1]**3 - radii[0]**3)) < 1e-7
    assert not (continuous / 'Diagnostic/VolumeConserve').exists()
    xi, counts, n24 = [], [], 0
    for rank in range(ranks):
        path = continuous / 'DATA' / str(rank) / ('PICES_P4.chkpt.%d.0' % rank)
        values, hosts, count, local_radial = particles(path)
        counts.append(count)
        if args.mode == 'cell_uniform':
            histogram = Counter(hosts)
            assert set(histogram) == set(range(1, (args.nodes-1)**2 * local_radial + 1))
            assert set(histogram.values()) == {27}, (rank, histogram)
        for r, host, flavor in zip(values[2], hosts, values[6]):
            layer = (rank % args.radial_partitions) * local_radial + (host - 1) % local_radial
            lo, hi = radii[layer:layer + 2]
            x = (r**3 - lo**3) / (hi**3 - lo**3)
            assert -2e-5 < x < 1 + 2e-5, (rank, host, r, lo, hi)
            xi.append(x)
            n24 += flavor == 24
    assert sum(counts) == total
    if args.mode == 'cell_uniform':
        assert set(counts) == {total // ranks}, counts
    elif args.radial_partitions == 2:
        assert max(counts) > 2 * min(counts), counts
    xi.sort()
    ks = max(max(abs(x - i / total), abs(x - (i + 1) / total)) for i, x in enumerate(xi))
    assert ks < 2.5 / math.sqrt(total), ('nonuniform within-cell radial volume sampling', ks)
    expected24 = (0.5813265437084149**3 - radii[0]**3) / (radii[-1]**3 - radii[0]**3)
    if args.mode == 'cell_uniform':
        expected24 = sum(max(0, min(1, (0.5813265437084149**3 - lo**3) / (hi**3 - lo**3)))
                         for lo, hi in zip(radii, radii[1:])) / (len(radii) - 1)
    assert abs(n24 / total - expected24) < 0.004, (n24 / total, expected24)

    split, _ = prepare('uniform_split', args.mode, stop=1)
    execute(split)
    comparisons += compare(continuous, split, ranks, (0, 1))
    original = (split / 'case.cfg').read_text()
    restart = settings(original, restart='on', solution_cycles_init=1, zero_elapsed_time='off',
                       datafile_old='PICES_P4', datadir_old='DATA/%RANK', maxstep=2, maxtotstep=3)
    (split / 'case.cfg').write_text(restart)
    errors = execute(split, tag='.restart')
    assert 'TRACER_SEED_GLOBAL' not in errors, 'restart must not regenerate particles'
    comparisons += compare(continuous, split, ranks, (2,))
    (split / 'case.cfg').write_text(settings(restart, tracer_seed=20261115))
    execute(split, tag='.bad_seed', failure='checkpoint metadata')
    other_mode = 'volume_uniform' if args.mode == 'cell_uniform' else 'cell_uniform'
    (split / 'case.cfg').write_text(settings(restart, tracer_seed_distribution=other_mode))
    execute(split, tag='.bad_mode', failure='checkpoint metadata')
    if args.mode == 'cell_uniform':
        (split / 'case.cfg').write_text(settings(restart, tracers_per_element=28))
        execute(split, tag='.bad_density', failure='checkpoint metadata')
    (split / 'case.cfg').write_text(restart)
    sparse, _ = prepare('reject_empty_quota', 'volume_uniform', density=1)
    execute(sparse, failure='empty initial cells')
    invalid, _ = prepare('reject_bad_mode', 'unknown')
    execute(invalid, failure='tracer_seed_distribution must be')
    invalid_ic, _ = prepare('reject_tracer_file_restart', args.mode)
    (invalid_ic / 'case.cfg').write_text(settings((invalid_ic / 'case.cfg').read_text(), tracer_ic_method=2))
    execute(invalid_ic, failure='require tracer_ic_method=0')
    previous_comparisons = previous_fingerprints = shared_points = 0
    if args.previous_build:
        current_c, _ = prepare('C_current', 'volume_uniform', stop=1)
        previous_c, _ = prepare('C_previous', 'volume_uniform', stop=1)
        execute(current_c)
        execute(previous_c, args.previous_build)
        previous_comparisons = compare(current_c, previous_c, ranks, (0, 1))
        comparisons += previous_comparisons
        for rank in range(ranks):
            for step in (0, 1):
                rel = Path('DATA') / str(rank) / ('PICES_P4.chkpt.%d.%d.pices.json' % (rank, step))
                left, right = (json.loads((case / rel).read_text()) for case in (current_c, previous_c))
                assert left['physics_mesh_sha256'] == right['physics_mesh_sha256'], (rank, step)
                previous_fingerprints += 1
        if args.mode == 'cell_uniform':
            for rank in range(ranks):
                rel = Path('DATA') / str(rank) / ('PICES_P4.chkpt.%d.0' % rank)
                grouped = []
                for case in (continuous, current_c):
                    values, hosts, _, _ = particles(case / rel)
                    cells = {}
                    for index, host in enumerate(hosts):
                        cells.setdefault(host, []).append(tuple(row[index] for row in values[:6]))
                    grouped.append(cells)
                assert grouped[0].keys() == grouped[1].keys()
                for host in grouped[0]:
                    b, c = grouped[0][host], grouped[1][host]
                    count = min(len(b), len(c))
                    assert b[:count] == c[:count], ('B/C sampler prefix differs', rank, host)
                    shared_points += count
    checks = ['fixed global particle count', 'actual host validation',
              'constant count per actual cell' if args.mode == 'cell_uniform' else 'FE floor/ceil quotas',
              'spherical shell coverage and within-cell radial distribution',
              'repeatable initial and step-1 particle state',
              'continuous/split restart physical equality',
              'changed seed and mode rejected by checkpoint identity',
              'unsupported tracer-file IC mode rejected explicitly',
              'empty volume quotas rejected without minimum-particle correction']
    if not args.skip_legacy:
        checks.append('legacy default/explicit physical payload equality')
    if args.previous_build:
        checks.append('C physical payloads and checkpoint fingerprints unchanged from C-only build')
        if args.mode == 'cell_uniform':
            checks.append('B/C share identical within-cell coordinate stream prefixes')
    if args.mode == 'cell_uniform':
        checks.append('changed per-cell density rejected by B checkpoint identity')
    result = dict(status='PASS', ranks=ranks, nodes=args.nodes, mode=args.mode, particles=total, rank_counts=counts,
                  initial_C24_count_fraction=n24 / total, expected_C24_fraction=expected24,
                  within_cell_r3_KS=ks, initialization=report,
                  compared_physical_payloads=comparisons, archived_legacy_checked=bool(args.legacy_build) and not args.skip_legacy,
                  previous_C_payloads=previous_comparisons, previous_C_fingerprints=previous_fingerprints,
                  shared_B_C_coordinate_prefix_points=shared_points,
                  checks=checks)
    (root / 'results.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
