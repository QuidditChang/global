"""Run the real EBA/RH20 solver with the volume audit off and on.

Requires tests/cbf/build_validation.py and runs/pices_p8a inputs. This is a
small-mesh MPI regression, not a replacement for the 384-rank HPC experiment.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import re
import shutil
import struct
import subprocess


def physical_payload(path):
    """Mask unused Fortran-style index-zero slots, as in check_rh_heating.py."""
    data = bytearray(path.read_bytes())
    if path.name.endswith('.state'):
        return data
    nx, ny, nz, *_ = struct.unpack_from('8i', data)
    nn, ne = nx * ny * nz, (nx - 1) * (ny - 1) * (nz - 1)
    off = 44 + 16
    for _ in range(2):
        data[off:off + 8] = bytes(8)
        off += 8 * (nn + 1)
    off += 16 + 8
    data[off:off + 8] = bytes(8)
    off += 8 * (ne + 1) + 8 * 3 * nn + 16
    basic, extra, flavors, last, count = struct.unpack_from('5i', data, off)
    off += 20
    assert basic == 12
    for _ in range(6 + extra):
        data[off:off + 8] = bytes(8)
        off += 8 * (count + 1)
    data[off:off + 4] = bytes(4)
    off += 4 * (count + 1) + 16
    ncomp = struct.unpack_from('i', data, off)[0]
    off += 4 + 16 * ncomp
    for _ in range(ncomp):
        data[off:off + 8] = bytes(8)
        off += 8 * (ne + 1)
    assert off == len(data), (path, off, len(data))
    return data


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('build', type=Path)
    p.add_argument('runs', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--radial-partitions', type=int, choices=(1, 2), default=1)
    args = p.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    inputs = args.runs / 'pices_p8a'
    text = (inputs / 'continuous.cfg').read_text()
    # The existing regression exercises genuine thermal/particle operators and
    # geological reclassification; RH20 makes its C24 heat coupling active.
    changes = {'RH': '20', 'maxstep': '2', 'maxtotstep': '3'}
    ranks = 12 * args.radial_partitions
    # Keep the same five-row reference state and global mesh; only split the
    # four radial cells between one or two MPI processes per cap.
    changes.update(nprocz=str(args.radial_partitions), mgunitz=str(4 // args.radial_partitions))
    for key, value in changes.items():
        pattern = r'^' + re.escape(key) + r'\s*=.*$'
        if re.search(pattern, text, re.M):
            text = re.sub(pattern, key + '=' + value, text, flags=re.M)
        else:
            text += '\n' + key + '=' + value + '\n'
    for name in ('off', 'on'):
        root = output / name
        root.mkdir()
        (root / 'DATA').mkdir()
        shutil.copy2(inputs / 'refstate_EBA_PICES_P4.txt', root)
        shutil.copytree(inputs / 'pices_p4_forcing', root / 'pices_p4_forcing',
                        ignore=shutil.ignore_patterns('._*'))
        (root / 'case.cfg').write_text(text)
        env = {k: v for k, v in os.environ.items() if not k.startswith('PICES_')}
        if name == 'on':
            env.update(PICES_VOLUME_AUDIT='1', PICES_VOLUME_AUDIT_INTERVAL='1')
        with (root / 'stdout').open('w') as out, (root / 'stderr').open('w') as err:
            result = subprocess.run([
                os.environ.get('MPIEXEC', '/usr/local/bin/mpiexec'),
                '--oversubscribe', '-n', str(ranks),
                str(args.build.resolve() / 'CitcomSFull'), 'case.cfg'],
                cwd=root, env=env, stdout=out, stderr=err, timeout=240)
        assert result.returncode in (0, 8), (name, result.returncode,
                                            (root / 'stderr').read_text()[-4000:])
        assert (root / 'DATA/0/PICES_P4.chkpt.0.2').is_file(), (
            name, 'solver did not finish step 2', (root / 'stderr').read_text()[-4000:])
    assert not (output / 'off/Diagnostic/VolumeAudit').exists()
    for rank in range(ranks):
        for step in (0, 2):
            for suffix in ('', '.pices.state'):
                rel = Path('DATA') / str(rank) / ('PICES_P4.chkpt.%d.%d%s' % (rank, step, suffix))
                assert physical_payload(output / 'off' / rel) == physical_payload(output / 'on' / rel), rel
    audit = output / 'on/Diagnostic/VolumeAudit'
    with (audit / 'globals.csv').open() as f:
        rows = list(csv.DictReader(f))
    accepted = [r for r in rows if r['status'] == 'accepted']
    assert len(accepted) == 3 * 25, len(accepted)
    assert {int(r['step']) for r in accepted} == {0, 1, 2}
    max_residual = 0.0
    for step in (0, 1, 2):
        batch = [r for r in accepted if int(r['step']) == step]
        shell = float(batch[0]['shell_volume_nd'])
        assert abs(sum(float(r['volume_after_nd']) for r in batch) / shell - 1) < 1e-12
        assert abs(sum(float(r['source_nd']) for r in batch)) / shell < 1e-12
        assert sum(int(r['count_before']) for r in batch) == sum(int(r['count_moved']) for r in batch)
        assert sum(int(r['count_moved']) for r in batch) == sum(int(r['count_after']) for r in batch)
        for row in batch:
            residual = abs(float(row['residual_nd'])) / shell
            max_residual = max(max_residual, residual)
            assert residual < 1e-10, row
        c24 = next(r for r in batch if int(r['flavor']) == 24)
        assert float(c24['source_nd']) == 0
        assert c24['count_before'] == c24['count_moved'] == c24['count_after']
    with (audit / 'radial.csv').open() as f:
        radial = [r for r in csv.DictReader(f) if r['status'] == 'accepted']
    for step in (0, 1, 2):
        batch = [r for r in radial if int(r['step']) == step]
        assert {int(r['layer']) for r in batch} == {1, 2, 3, 4}
        global_row = next(r for r in accepted if int(r['step']) == step and int(r['flavor']) == 24)
        shell = float(global_row['shell_volume_nd'])
        assert abs(sum(float(r['shell_volume_nd']) for r in batch) - shell) / shell < 1e-12
        assert abs(sum(float(r['C24_after_nd']) for r in batch) - float(global_row['volume_after_nd'])) / shell < 1e-12
    result = dict(status='PASS', ranks=ranks, accepted_steps=[0, 1, 2],
                  compared_checkpoint_payloads=4 * ranks, normalized_budget_residual_max=max_residual,
                  checks=['audit off/on identical physical checkpoint payloads',
                          'all 25 flavors sum to shell volume', 'classification net balance',
                          'global particle count balance', 'immutable C24',
                          'successful Stokes commit markers', 'radial partition integration'])
    (output / 'results.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
