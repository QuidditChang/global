"""Prove that the complete volume-shadow path leaves coupled physics unchanged.

Run the same frozen native binary on the same manufactured 9^3 EBA/RH20
fixture with volume mode off and shadow. Compare every accepted checkpoint,
including thermal, velocity, particle, composition and persistent history.
Only the explicitly named new volume weight and authenticated schema-7 volume
ledger are excluded. Unused index-zero slots use the established normalizer.
This tests software noninterference, not the remapper's scientific accuracy.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import struct
import subprocess
import sys

from check_volume_audit_integration import physical_payload
from check_volume_control_integration import fixture


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def native_config(path):
    result = {}
    for line in path.read_text().splitlines():
        line = line.split('#', 1)[0].strip()
        if '=' in line:
            key, value = line.split('=', 1)
            result[key.strip()] = value.strip()
    return result


def enabled(value):
    return str(value).lower() in ('1', 'on', 'true', 'yes')


class Reader:
    def __init__(self, data):
        self.data, self.offset, self.fields = data, 0, {}

    def take(self, name, size):
        if size < 0 or self.offset + size > len(self.data):
            raise ValueError('Truncated field: ' + name)
        value = self.data[self.offset:self.offset + size]
        self.fields[name] = value
        self.offset += size
        return value

    def finish(self):
        if self.offset != len(self.data):
            raise ValueError('Unparsed checkpoint bytes: %d' % (len(self.data) - self.offset))
        return self.fields


def native_fields(path, metadata):
    """Parse all native fields, dropping only the named extra quantity."""
    r = Reader(physical_payload(path))
    header = struct.unpack('8i', r.take('mesh_and_step', 32))
    nx, ny, nz, px, py, pz, caps, step = header
    assert list((nx, ny, nz)) == metadata['local_mesh'] and caps == 1
    assert [px, py, pz] == metadata['decomposition'] and step == metadata['step']
    nn, ne = nx * ny * nz, (nx - 1) * (ny - 1) * (nz - 1)
    r.take('clock', 12)
    r.take('energy_sentinel', 16)
    r.take('temperature', 8 * (nn + 1))
    r.take('temperature_derivative', 8 * (nn + 1))
    r.take('momentum_sentinel', 16)
    r.take('velocity_counters', 8)
    r.take('pressure', 8 * (ne + 1))
    r.take('velocity_solution', 24 * nn)
    r.take('tracer_sentinel', 16)
    basic, extra, flavors, last, count = struct.unpack('5i', r.take('tracer_header', 20))
    assert basic == metadata['basic_quantities'] == 12
    assert count == metadata['particles'] and flavors == metadata['flavors'] == 25
    entries = metadata['extraq']
    assert len(entries) == extra and sorted(q['slot'] for q in entries) == list(range(extra))
    names = {q['slot']: q['name'] for q in entries}
    assert len(set(names.values())) == extra
    weights = [q for q in entries if q['name'] == 'volume_weight']
    assert len(weights) == (metadata['schema'] == 7)
    r.fields['tracer_header'] = struct.pack('5i', basic, extra - len(weights), flavors, last, count)
    for index in range(6):
        r.take('particle_basic_%d' % index, 8 * (count + 1))
    for slot in range(extra):
        key = 'particle_extra_' + names[slot]
        raw = r.take(key, 8 * (count + 1))
        if names[slot] == 'volume_weight':
            assert all(math.isfinite(w[0]) and w[0] > 0 for w in struct.iter_unpack('d', raw[8:]))
            del r.fields[key]
    r.take('particle_host', 4 * (count + 1))
    r.take('composition_sentinel', 16)
    ncomp, = struct.unpack('i', r.take('composition_count', 4))
    assert ncomp == 24
    r.take('bulk_composition', 8 * ncomp)
    r.take('initial_bulk_composition', 8 * ncomp)
    for component in range(ncomp):
        r.take('cell_composition_%d' % (component + 1), 8 * (ne + 1))
    return r.finish()


def state_fields(path, metadata, cfg):
    """Parse actual conditional state sections; never infer a prefix by length."""
    r = Reader(path.read_bytes())
    schema = metadata['schema']
    assert schema in (6, 7) and float(cfg['RH']) != 1
    assert enabled(cfg['pices_p4']) and metadata['flavors'] == 25
    nx, ny, nz = metadata['local_mesh']
    nn, ne = nx * ny * nz, (nx - 1) * (ny - 1) * (nz - 1)
    for component in range(3):
        r.take('accepted_velocity_%d' % component, 4 * nn)
    r.take('heat_residual_valid', 4)
    r.take('heat_residual', 64 * ne)
    r.take('trench_visit_age', 4)
    count, = struct.unpack('i', r.take('flavor_history_count', 4))
    expected = int(cfg['nodex']) * int(cfg['nodey']) if enabled(cfg['tracer_reclassify_flavors']) else 0
    assert count == expected
    r.take('flavor_history', 4 * count)
    coupled = (int(cfg['rheol']) == 7 or enabled(cfg.get('file_vbcs', 'off')) or
               int(cfg.get('qvis_mode', '0')) != 0 or float(cfg['fixed_timestep']) <= 0)
    if coupled:
        r.take('accepted_element_viscosity', 32 * ne)
        r.take('accepted_node_viscosity', 4 * nn)
    for component in range(3):
        r.take('TA_boundary_%d' % component, 4 * nn)
    if enabled(cfg['lith_age']):
        r.take('frozen_trench_geometry', 4 * int(cfg['nodex']) * int(cfg['nodey']))
    r.take('RH_state', 32)
    if schema == 7:
        assert metadata['composition_volume']['mode'] == 'shadow'
        slot = next(q['slot'] for q in metadata['extraq'] if q['name'] == 'volume_weight')
        header = struct.unpack('4i', r.take('volume_header', 16))
        assert header == (0x564f4c31, 1, 1, slot)
        baseline = struct.unpack('25d', r.take('volume_baseline', 200))
        sources = struct.unpack('25d', r.take('volume_sources', 200))
        assert all(math.isfinite(v) and v >= 0 for v in baseline) and sum(baseline) > 0
        assert all(math.isfinite(v) for v in sources) and sources[24] == 0
        assert abs(sum(sources)) <= 1e-10 * sum(baseline)
        for name in ('volume_header', 'volume_baseline', 'volume_sources'):
            del r.fields[name]
    else:
        assert 'composition_volume' not in metadata
    return r.finish()


def checkpoint_set(root, step, ranks):
    records, hashes, manifests = [], [], []
    for rank in range(ranks):
        stem = root / 'DATA' / str(rank) / ('PICES_P4.chkpt.%d.%d' % (rank, step))
        metadata_path = Path(str(stem) + '.pices.json')
        metadata = json.loads(metadata_path.read_text())
        assert metadata['rank'] == rank and metadata['step'] == step and metadata['mpi_size'] == ranks
        assert metadata['phase'] == 'accepted' and metadata['magic'] == 'CITCOMS_EBA_PICES'
        state = Path(str(stem) + '.pices.state')
        assert digest(stem) == metadata['checkpoint_sha256']
        assert digest(state) == metadata['accepted_velocity_sha256']
        hashes.append(digest(metadata_path).encode() + b'\0')
        manifests.append(json.loads(Path(str(stem) + '.pices.manifest').read_text()))
        records.append((stem, state, metadata))
    expected = dict(magic='CITCOMS_EBA_PICES_COMPLETE', schema=1, mpi_size=ranks,
                    metadata_set_sha256=hashlib.sha256(b''.join(hashes)).hexdigest())
    assert all(m == expected for m in manifests)
    return records


def equal_fields(left, right, context):
    if left.keys() != right.keys():
        raise ValueError('%s: field set differs: %r' % (context, left.keys() ^ right.keys()))
    for name in left:
        if left[name] != right[name]:
            first = next((i for i, (a, b) in enumerate(zip(left[name], right[name])) if a != b),
                         min(len(left[name]), len(right[name])))
            raise ValueError('%s: %s differs at byte %d' % (context, name, first))
    return len(left)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('build', type=Path)
    p.add_argument('runs', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--radial-partitions', choices=(1, 2), type=int, default=1)
    p.add_argument('--hops', type=int, default=4)
    p.add_argument('--iterations', type=int, default=1000)
    args = p.parse_args(argv)
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    binary = root / 'CitcomSFull'
    before = digest(args.build / 'CitcomSFull')
    shutil.copy2(args.build / 'CitcomSFull', binary)
    assert digest(binary) == before
    ranks = 12 * args.radial_partitions
    env = {k: v for k, v in os.environ.items() if not k.startswith('PICES_')}
    env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    configs = {}
    for mode in ('off', 'shadow'):
        case = root / mode
        fixture(args.runs.resolve() / 'pices_p8a', case, args.radial_partitions)
        text = (case / 'case.cfg').read_text()
        changes = dict(maxstep='2', maxtotstep='3', checkpointFrequency='1',
                       pices_volume_mode=mode, pices_volume_hops=str(args.hops),
                       pices_volume_iterations=str(args.iterations), pices_volume_tolerance='1e-12')
        for key, value in changes.items():
            pattern = '^' + re.escape(key) + r'\s*=.*$'
            text = re.sub(pattern, key + '=' + value, text, flags=re.M) if re.search(pattern, text, re.M) else text + '\n' + key + '=' + value + '\n'
        (case / 'case.cfg').write_text(text)
        configs[mode] = native_config(case / 'case.cfg')
    inputs = [Path('refstate_EBA_PICES_P4.txt')]
    inputs += [path.relative_to(root / 'off') for path in (root / 'off/pices_p4_forcing').rglob('*') if path.is_file()]
    if args.radial_partitions == 2:
        inputs.append(Path('control_radial.coor'))
    input_hashes = {}
    for relative in inputs:
        left, right = digest(root / 'off' / relative), digest(root / 'shadow' / relative)
        assert left == right, ('input mismatch', relative)
        input_hashes[str(relative)] = left
    for mode in ('off', 'shadow'):
        case = root / mode
        current_env = dict(env)
        if mode == 'shadow':
            current_env.update(PICES_VOLUME_AUDIT='1', PICES_VOLUME_AUDIT_INTERVAL='1')
        with (case / 'stdout').open('w') as stdout, (case / 'stderr').open('w') as stderr:
            result = subprocess.run([os.environ.get('MPIEXEC', '/usr/local/bin/mpiexec'),
                                     '--oversubscribe', '-n', str(ranks), str(binary), 'case.cfg'],
                                    cwd=case, env=current_env, stdout=stdout, stderr=stderr, timeout=1800)
        assert result.returncode in (0, 8), (mode, result.returncode, str(case / 'stderr'))
    assert {k: v for k, v in configs['off'].items() if k != 'pices_volume_mode'} == {
        k: v for k, v in configs['shadow'].items() if k != 'pices_volume_mode'}
    compared, native_count, state_count = 0, 0, 0
    for step in range(3):
        off = checkpoint_set(root / 'off', step, ranks)
        shadow = checkpoint_set(root / 'shadow', step, ranks)
        for rank, ((op, os_, om), (sp, ss, sm)) in enumerate(zip(off, shadow)):
            assert om['schema'] == 6 and sm['schema'] == 7
            excluded = {'schema', 'composition_volume', 'extraq', 'checkpoint_sha256',
                        'accepted_velocity_sha256', 'physics_mesh_sha256'}
            assert {k: v for k, v in om.items() if k not in excluded} == {
                k: v for k, v in sm.items() if k not in excluded}
            assert [q for q in sm['extraq'] if q['name'] != 'volume_weight'] == om['extraq']
            context = 'step%d rank%d' % (step, rank)
            native_count += equal_fields(native_fields(op, om), native_fields(sp, sm), context + ' native')
            state_count += equal_fields(state_fields(os_, om, configs['off']),
                                        state_fields(ss, sm, configs['shadow']), context + ' persistent state')
            compared += 2
    with (root / 'shadow/Diagnostic/VolumeConserve/steps.csv').open() as stream:
        rows = [r for r in csv.DictReader(stream) if r['status'] == 'accepted']
    assert [int(r['step']) for r in rows] == [0, 1, 2]
    assert not (root / 'off/Diagnostic/VolumeConserve').exists()
    result = dict(status='PASS', purpose='Coupled shadow noninterference; not scientific remap approval',
                  binary_sha256=before, ranks=ranks, steps=[0, 1, 2], schemas={'off': 6, 'shadow': 7},
                  common_input_sha256=input_hashes, rheology=int(configs['off']['rheol']),
                  identical_physical_payloads=compared, identical_native_fields=native_count,
                  identical_persistent_fields=state_count, all_checkpoint_hashes_verified=True,
                  shadow_capacity_failures=[int(r['step']) for r in rows if r['feasible'] != '1'],
                  exclusions=['Unused native index-zero slots', 'Metadata-named volume_weight extraq and its count',
                              'Parsed and validated schema-7 VOL1 header, baseline and source ledger (416 bytes)'],
                  environments='Shadow additionally enables read-only phase-1 VolumeAudit; off enables neither volume diagnostic')
    (root / 'results.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
    return 0


if __name__ == '__main__':
    sys.exit(main())
