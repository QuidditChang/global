"""Audit the P0 standalone PG result; does not certify PICES or energy closure.

python3 tests/pices/verify_p0.py /path/to/downloaded/job_DIRECTORY
Use --local only for a local validation build without the LSF provenance files.
Requires Python >= 3.8 and only the standard library.
"""
import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import struct
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'cbf'))
from verify_native_outputs import verify as verify_native
from verify_thermal_budget import verify as verify_budget


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_nodes(path, step):
    with gzip.open(path, 'rt') as stream:
        lines = stream.read().splitlines()
    header = lines[0].split()
    require(len(header) == 6, str(path) + ': invalid header')
    require(int(header[0]) == step and int(header[1]) == 125, 'wrong step/node count')
    numbers = list(map(float, header[2:]))
    require(all(map(math.isfinite, numbers)), 'nonfinite header')
    expected_time = struct.unpack('f', struct.pack('f', 1e-6))[0] * step
    require(math.isclose(numbers[0], expected_time, rel_tol=1e-5, abs_tol=1e-14), 'wrong field clock')
    require(numbers[1:] == [300., 3700., 3400.], 'wrong temperature normalization')
    require(list(map(int, lines[1].split())) == [1, 125], 'wrong local cap header')
    require(len(lines) == 127, 'missing/extra node rows: ' + str(path))
    rows = [list(map(float, line.split())) for line in lines[2:]]
    for index, row in enumerate(rows, 1):
        require(len(row) == 4 and all(map(math.isfinite, row)), 'invalid node: ' + str(path))
        require(row[3] >= 0, 'negative Kelvin temperature')
        if index % 5 == 1:
            require(abs(row[3] - 3700.) <= 1e-6, 'bottom Dirichlet violated')
        if index % 5 == 0:
            require(abs(row[3] - 300.) <= 1e-6, 'top Dirichlet violated')
    return rows


def verify(root, local=False):
    root = Path(root)
    status = int((root/'mpi_exit_code.txt').read_text())
    require(status in (0, 8), 'MPI failure: %d' % status)
    cfg = {}
    for line in (root/'cmbhf_EBA_PICES_P0.cfg').read_text().splitlines():
        line = line.split('#', 1)[0].strip()
        if not line:
            continue
        key, value = map(str.strip, line.split('=', 1))
        require(key not in cfg, 'duplicate config key: ' + key)
        cfg[key] = value
    expected = dict(solver='full', Solver='cgrad', nodex='5', nodey='5', nodez='5',
                    nproc_surf='12', nprocx='1', nprocy='1', nprocz='1',
                    maxstep='2', maxtotstep='3', fixed_timestep='1e-6',
                    Q0='0.3', dissipation_number='0.1', tracer='0', lith_age='0',
                    compressible_formulation='eba', filter_temp='off',
                    phase_delta_s='0,0,0', checkpointFrequency='2',
                    datadir='DATA/%RANK', output_format='ascii-gz',
                    temperature_audit='on', monitor_max_T='off',
                    CBF_frequency='1', kd_upper_prefactor='4', kd_lower_prefactor='4',
                    kd_upper_linear='0', kd_upper_quadratic='0', kd_lower_linear='0',
                    kd_lower_quadratic='0', kT_exponent='0', kC_ratio='1')
    for key, value in expected.items():
        require(cfg.get(key) == value, 'unexpected config: ' + key)
    require('energy_solver' not in cfg, 'P0 has no energy_solver parameter')
    reference = [list(map(float, s.split())) for s in
                 (root/'refstate_EBA_PICES_P0.txt').read_text().splitlines() if s.strip()]
    require(reference == [[1., 1., t, 1., 1.] for t in (1., .75, .5, .25, 0.)],
            'unexpected P0 reference state')
    require({p.name for p in (root/'DATA').iterdir() if p.is_dir()} ==
            {str(i) for i in range(12)}, 'missing or unexpected rank directories')
    temperature = {step: [] for step in (0, 1, 2)}
    for rank in range(12):
        directory = root/'DATA'/str(rank)
        log = (directory/'log').read_text()
        exits = re.findall(r'^TEMP_AUDIT stage=thermal_exit (.*)$', log, re.M)
        require(len(exits) == 2, 'wrong accepted-step count at rank %d' % rank)
        for step, record in enumerate(exits, 1):
            values = dict(re.findall(r'(\w+)=([^\s]+)', record))
            require(int(values['step']) == step and int(values['rank']) == rank,
                    'wrong audit step/rank')
            require(values['negative'] == '0' and values['nonfinite'] == '0',
                    'invalid temperature audit')
            require(math.isclose(float(values['time_nd']), step*1e-6,
                                 rel_tol=1e-6, abs_tol=1e-14), 'wrong accepted clock')
            require(math.isclose(float(values['dt_nd']), 1e-6,
                                 rel_tol=1e-6, abs_tol=1e-14), 'wrong accepted dt')
        require('TEMP_AUDIT_NODE' not in log, 'invalid intermediate temperature')
        for step in (0, 1, 2):
            rows = read_nodes(directory/str(step)/('velo.%d.%d.gz' % (rank, step)), step)
            temperature[step].extend(row[3] for row in rows)
        for step in (0, 2):
            path = directory/('PICES_P0.chkpt.%d.%d' % (rank, step))
            require(path.stat().st_size > 44, 'truncated checkpoint')
            with path.open('rb') as stream:
                header = struct.unpack('=8i3f', stream.read(44))
            require(header[:8] == (5,5,5,1,1,1,1,step), 'wrong checkpoint header')
            require(math.isclose(header[8], step*1e-6, rel_tol=1e-6, abs_tol=1e-14),
                    'wrong checkpoint time')
    verify_budget(root, [0,1,2], constant_q0=.3)
    for step in (0,1,2):
        verify_native(root, step, 12)
    provenance = {}
    if not local:
        for filename in ('solver_commit.txt','runs_commit.txt','pices-p0-build.commit'):
            value = (root/filename).read_text().strip()
            require(re.fullmatch('[0-9a-f]{40}',value) is not None, 'invalid ' + filename)
            provenance[filename] = value
        require(provenance['solver_commit.txt'] == provenance['pices-p0-build.commit'],
                'build/source commit mismatch')
        require((root/'runtime_source_diff.txt').read_text() == '', 'P0 runtime source changed')
        require((root/'launcher_exit_code.txt').read_text().strip() == '0', 'launcher failed')
        for filename in ('pices-p0-build.log','binary.sha256','mpi_version.txt','platform.txt',
                         'solver_status.txt','runs_status.txt','submitted.lsf'):
            require((root/filename).is_file(), 'missing provenance: ' + filename)
        for filename in ('pices-p0-build.log','mpi_version.txt','platform.txt','submitted.lsf'):
            require((root/filename).stat().st_size > 0, 'empty provenance: ' + filename)
        require(re.fullmatch(r'[0-9a-f]{64}\s+\S[^\n]*\n?',
                             (root/'binary.sha256').read_text()) is not None,
                'invalid binary SHA256 record')
        entries = (root/'input.sha256').read_text().splitlines()
        require(len(entries) == 2, 'incomplete input hashes')
        for entry, filename in zip(entries, ('cmbhf_EBA_PICES_P0.cfg','refstate_EBA_PICES_P0.txt')):
            digest, recorded_name = entry.split(maxsplit=1)
            require(recorded_name.strip() == filename, 'unexpected input hash name')
            require(hashlib.sha256((root/filename).read_bytes()).hexdigest() == digest,
                    'input hash mismatch')
    return dict(status='PASS', scope='P0 standalone PG baseline only; no PICES/energy-closure certification',
                provenance_checked=not local, mpi_exit_code=status, ranks=12, steps=[0,1,2],
                temperature_K={str(s):dict(min=min(v),max=max(v)) for s,v in temperature.items()},
                provenance=provenance)


if __name__ == '__main__':
    if not __debug__:
        sys.exit('Do not use python -O: upstream CBF validators require assertions.')
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('root', type=Path)
    p.add_argument('--local', action='store_true')
    p.add_argument('--summary', type=Path)
    args = p.parse_args()
    try:
        result = verify(args.root, args.local)
    except (ValueError, OSError, EOFError, AssertionError, KeyError, IndexError, struct.error) as error:
        result = dict(status='FAIL', error=str(error))
    rendered = json.dumps(result, indent=2, ensure_ascii=False) + '\n'
    if args.summary:
        args.summary.write_text(rendered)
    print(rendered)
    sys.exit(0 if result['status'] == 'PASS' else 1)
