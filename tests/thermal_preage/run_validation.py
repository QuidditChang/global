"""Exercise initialization-only MPI driver and independently check conduction.

Requires NumPy/SciPy, an MPI implementation, and build_validation.py output.
Synthetic reference/age/trench inputs are created locally; production forcing
and external reconstruction files are never read. Exit status is nonzero on a
missing snapshot, mutation, wrong BC, rank dependence or numerical regression.
"""
import argparse
import json
import os
from pathlib import Path
import re
import shlex
import subprocess

import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import spsolve

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def parameters():
    result = {}
    for line in (ROOT / 'tests/cbf/smoke.cfg').read_text().splitlines():
        line = line.split('#', 1)[0].strip()
        if line:
            key, value = line.split('=', 1)
            result[key] = value
    result.update({
        'nodex': 5, 'nodey': 5, 'nodez': 17, 'mgunitx': 4, 'mgunity': 4, 'mgunitz': 16,
        'datadir': 'DATA', 'datafile': 'preage_test', 'minstep': 0, 'maxstep': 0,
        'maxtotstep': 0, 'CBF_frequency': 0, 'output_q_surf_CBF': 'off',
        'output_q_botm_CBF': 'off', 'write_q_files': 0, 'energy_solver': 'pices',
        'pices_eba': 'on', 'pices_p4': 'on', 'pices_checkpoint': 'off',
        'tracer': 'on', 'tracer_flavors': 2, 'tracer_reclassify_flavors': 'off',
        'tracers_per_element': 8, 'chemical_buoyancy': 'on', 'buoy_type': 1,
        'buoyancy_ratio': 0, 'z_interface': .75, 'ic_method_for_flavors': 0,
        'lith_age': 1, 'lith_age_time': 1, 'lith_age_asml': 1,
        'lith_age_depth': .04, 'lith_age_asml_tau_Ma': 1, 'max_plate_age_Ma': 70,
        'lith_age_file': 'forcing/age.', 'flag_depth_file': 'forcing/trench.',
        'flag_depth_new_file': 'forcing/trench.', 'tf_file': 'forcing/trench.',
        'temperature_bound_adj': 'off', 'bottom_tbl_thickness': 0,
        'start_age': 0, 'reset_startage': 'off', 'zero_elapsed_time': 'off',
        'kT_exponent': 0, 'kC_ratio': 1, 'kC_primordial_flavor': 1,
        'kd_upper_prefactor': 4, 'kd_lower_prefactor': 4,
        'kd_upper_linear': 0, 'kd_upper_quadratic': 0,
        'kd_lower_linear': 0, 'kd_lower_quadratic': 0,
        'Q0': 0, 'dissipation_number': .1, 'perturbmag': 0,
        'output_optional': '', 'thermal_preage_Ma': 1000,
        'thermal_preage_max_dt_Ma': 250,
    })
    return result


def prepare(directory, changes):
    directory.mkdir(parents=True)
    (directory / 'DATA').mkdir()
    cfg = parameters()
    cfg.update(changes)
    cfg['mgunitz'] = (int(cfg['nodez']) - 1) // int(cfg['nprocz'])
    (directory / 'case.cfg').write_text('\n'.join(f'{k}={v}' for k, v in cfg.items()) + '\n')
    nz = int(cfg['nodez'])
    # Thermodynamic background deliberately differs from actual radial BCs.
    (directory / 'refstate.txt').write_text('\n'.join(
        f'1 1 {0.68 - .28*i/(nz-1):.17g} 1 1' for i in range(nz)) + '\n')
    forcing = directory / 'forcing'
    forcing.mkdir()
    for age in (0, 1):
        (forcing / f'trench.{age}.xyz').write_text('')
        for cap in range(12):
            (forcing / f'age.{age}.{cap}').write_text('70\n' * (int(cfg['nodex']) * int(cfg['nodey'])))
    return cfg


def data(directory, stage, ranks):
    return [np.loadtxt(directory / f'nodes.{stage}.{r}.txt', ndmin=2) for r in range(ranks)]


def merged(rows):
    unique = {}
    maps = []
    values = []
    for block in rows:
        ids = []
        for row in block:
            key = tuple(np.round(row[2:5], 10))
            if key not in unique:
                unique[key] = len(values)
                values.append(row)
            else:
                previous = values[unique[key]]
                assert abs(previous[7] - row[7]) < 2e-9, ('shared temperature mismatch', key)
            ids.append(unique[key])
        maps.append(np.array(ids))
    return np.asarray(values), maps, unique


def field(directory, stage, ranks):
    rows, _, unique = merged(data(directory, stage, ranks))
    return {key: rows[index, 7] for key, index in unique.items()}


def distance(a, b):
    assert a.keys() == b.keys(), 'different global meshes'
    return float(np.max(np.abs([a[key] - b[key] for key in a])))


def inspect(directory, ranks, enabled):
    meta = dict(line.split() for line in (directory / 'metadata.0.txt').read_text().splitlines())
    assert int(meta['enabled']) == enabled
    assert int(meta['step']) == 0 and float(meta['clock']) == 0
    assert int(meta['initialized']) == 1
    final = data(directory, 'final', ranks)
    max_anomaly_error = max(float(np.max(abs(x[:, 8] - (x[:, 7] - x[:, 9])))) for x in final)
    assert max_anomaly_error < 2e-13, ('DataT initialized before final T', max_anomaly_error)
    for rank, block in enumerate(final):
        assert np.isfinite(block[:, 7]).all()
        bottom = block[:, 5] == 1
        top = block[:, 5] == int(meta['noz'])
        assert not bottom.any() or np.max(abs(block[bottom, 7] - 1)) < 1e-14
        assert not top.any() or np.max(abs(block[top, 7])) < 1e-14
        particles = np.loadtxt(directory / f'particles.{rank}.txt', ndmin=2)
        assert len(particles) > 0
        assert np.max(abs(particles[:, 2] - particles[:, 3])) < 2e-12, 'Tp not sampled from final temperature'
    result = {'rank_count': ranks, 'dataT_max_error': max_anomaly_error,
              'particles_match_final_T': True, 'production_clock': float(meta['clock'])}
    if not enabled:
        assert not list(directory.glob('isolation.*.txt')), 'disabled initializer was called'
        return result
    preaged = data(directory, 'preaged', ranks)
    before = data(directory, 'before', ranks)
    for rank, (old, new, last) in enumerate(zip(before, preaged, final)):
        bottom = new[:, 5] == 1
        top = new[:, 5] == int(meta['noz'])
        assert not bottom.any() or np.max(abs(new[bottom, 7] - 1)) < 1e-14
        assert not top.any() or np.max(abs(new[top, 7] - float(meta['surface']))) < 1e-14
        # The top temperature is temporary: the production TB values/flags stay unchanged.
        assert np.array_equal(old[:, 11:15], new[:, 11:15]), 'production BC metadata mutated'
        for line in (directory / f'isolation.{rank}.txt').read_text().splitlines():
            label, start, end = line.split()
            assert start == end, ('pre-age mutation', rank, label)
        deep = new[:, 6] < .95
        assert np.max(abs(new[deep, 7] - last[deep, 7])) < 1e-14, 'shallow overlay changed deep solution'
    logs = '\n'.join(path.read_text() for path in (directory / 'DATA').rglob('log'))
    ends = [dict(re.findall(r'(\w+)=([^\s]+)', line)) for line in logs.splitlines()
            if line.startswith('THERMAL_PREAGE_END ')]
    assert len(ends) >= 1, 'pre-age completion log missing'
    assert all(abs(float(x['age_Ma']) - float(meta['duration_Ma'])) < 1e-7 for x in ends)
    assert all(abs(float(x['balance'])) < 2e-8 for x in ends), 'frozen-capacity heat ledger failed'
    progress = [dict(re.findall(r'(\w+)=([^\s]+)', line)) for line in logs.splitlines()
                if line.startswith('THERMAL_PREAGE_PROGRESS ')]
    assert all(float(x['property_fraction']) <= 1 + 1e-10 for x in progress), 'accepted property bound exceeded'
    assert all(float(x['residual']) <= 1.01e-10 for x in ends), 'nonlinear residual too large'
    result.update(steps=int(ends[0]['steps']), rejected=int(ends[0]['rejected']),
                  max_residual=max(float(x['residual']) for x in ends),
                  balance=float(ends[0]['balance']), isolated_state_groups=8)
    return result


def independent_be(directory, ranks, duration, exponent=0., radial_only=False, phase_capacity=False):
    """Independent global SciPy solve from exported geometry, no production K.

    This deliberately checks full 3-D coupling, constrained solve, shared-node
    assembly, and k(newT) iteration, while relying on the real FE geometry.
    Phase-zero constant rho/cp fixture makes capacity independently verifiable.
    """
    initial, maps, _ = merged(data(directory, 'before', ranks))
    size = len(initial)
    old = initial[:, 7].copy()
    nz = int(initial[:, 5].max())
    fixed = (initial[:, 5] == 1) | (initial[:, 5] == nz)
    old[initial[:, 5] == 1] = 1.
    old[initial[:, 5] == nz] = .4
    free = np.flatnonzero(~fixed)
    bc = np.flatnonzero(fixed)
    ids, shape, gradients, weights, conductivity, rho_cp, radii = [], [], [], [], [], [], []
    for rank in range(ranks):
        quadrature = np.loadtxt(directory / f'quadrature.{rank}.txt', ndmin=2)
        block = quadrature[:, 5:-1].reshape(-1, 8, 7)
        ids.append(maps[rank][block[:, :, 0].astype(int) - 1])
        shape.append(block[:, :, 1])
        grad = block[:, :, 2:5].copy()
        grad[:, :, 0] *= quadrature[:, 3, None]
        grad[:, :, 1] *= quadrature[:, 3, None] / np.sin(quadrature[:, 4, None])
        if radial_only:
            grad[:, :, :2] = 0
        gradients.append(grad)
        weights.extend(quadrature[:, 2])
        conductivity.extend(quadrature[:, -1])
        radii.extend(1 / quadrature[:, 3])
        rho_cp.extend((block[:, :, 1] * block[:, :, 5]).sum(axis=1) *
                      (block[:, :, 1] * block[:, :, 6]).sum(axis=1))
    ids, shape, gradients = map(np.concatenate, (ids, shape, gradients))
    weights, conductivity, rho_cp = map(np.asarray, (weights, conductivity, rho_cp))
    meta = dict(line.split() for line in (directory / 'metadata.0.txt').read_text().splitlines())
    if phase_capacity:
        phases = np.loadtxt(directory / 'phase.0.txt', ndmin=2)
        tg_old = (old[ids] * shape).sum(axis=1)
        for depth, entropy, clapeyron, transT, inv_width in phases:
            q = (float(meta['outer_radius']) - np.asarray(radii) - depth) - clapeyron * (tg_old - transT)
            fraction = .5 * (1 + np.tanh(inv_width * q))
            dX_dT = -2 * inv_width * fraction * (1 - fraction) * clapeyron
            rho_cp += (tg_old + float(meta['surface_temp'])) * entropy * dX_dT
    mass = np.bincount(ids.ravel(), weights=(weights[:, None] * rho_cp[:, None] * shape).ravel(), minlength=size)
    dt = duration / float(meta['scalet'])
    gram = np.einsum('gik,gjk->gij', gradients, gradients)
    rr = np.broadcast_to(ids[:, :, None], gram.shape).ravel()
    cc = np.broadcast_to(ids[:, None, :], gram.shape).ravel()
    trial = old.copy()
    for _ in range(100):
        tg = (trial[ids] * shape).sum(axis=1)
        kt = (300 / np.maximum(300., 300. + 3400. * tg)) ** exponent
        matrix = coo_matrix((((weights * conductivity * kt)[:, None, None] * gram).ravel(), (rr, cc)), shape=(size, size)).tocsr()
        system = diags(mass) + dt * matrix
        updated = old.copy()
        updated[free] = spsolve(system[free][:, free], mass[free] * old[free] - system[free][:, bc] @ old[bc])
        change = np.max(abs(updated - trial))
        trial = updated
        if change < 2e-13:
            break
    else:
        raise AssertionError('independent nonlinear reference did not converge')
    observed, _, unique = merged(data(directory, 'preaged', ranks))
    assert np.array_equal(initial[:, 2:5], observed[:, 2:5])
    return trial, observed[:, 7], {'nodes': size, 'reference_max_error': float(np.max(abs(trial - observed[:, 7])))}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--build-dir', required=True, type=Path)
    p.add_argument('--output', required=True, type=Path)
    p.add_argument('--timeout', type=float, default=240)
    p.add_argument('--prepare-only', action='store_true')
    a = p.parse_args()
    root = a.output.resolve()
    root.mkdir(parents=True, exist_ok=True)
    exe = a.build_dir.resolve() / 'PreageInitialization'
    results = {}
    def run(name, changes=None, ranks=12, env=None, error=False):
        directory = root / name
        cfg = prepare(directory, changes or {})
        if a.prepare_only:
            return directory
        command = shlex.split(os.environ.get('MPIEXEC', 'mpiexec'))
        command += shlex.split(os.environ.get('MPIEXEC_FLAGS', ''))
        command += ['-n', str(ranks), str(exe), 'case.cfg']
        environment = dict(os.environ, OMPI_ALLOW_RUN_AS_ROOT='1', OMPI_ALLOW_RUN_AS_ROOT_CONFIRM='1', OPENBLAS_NUM_THREADS='1')
        environment.update(env or {})
        with (directory / 'stdout').open('w') as out, (directory / 'stderr').open('w') as err:
            completed = subprocess.run(command, cwd=directory, env=environment, stdout=out, stderr=err, timeout=a.timeout)
        stderr = (directory / 'stderr').read_text()
        if error:
            assert completed.returncode != 0 and 'thermal_preage' in stderr, (name, completed.returncode, stderr[-3000:])
            assert 'PREAGE_TEST_DONE' not in stderr
            results[name] = {'status': 'PASS', 'guard_rejected': True, 'exit': completed.returncode}
        else:
            assert completed.returncode == 0 and 'PREAGE_TEST_DONE' in stderr, (name, completed.returncode, stderr[-4000:])
            results[name] = dict(status='PASS', **inspect(directory, ranks, cfg.get('thermal_preage', 'off') == 'on'))
        (root / 'summary.json').write_text(json.dumps(results, indent=2) + '\n')
        print(name, results.get(name), flush=True)
        return directory

    if a.prepare_only:
        run('enabled_12', {'thermal_preage': 'on'})
        run('enabled_24', {'thermal_preage': 'on', 'nprocz': 2})
        return

    default = run('default_off')
    explicit = run('explicit_off', {'thermal_preage': 'off'})
    assert distance(field(default, 'final', 12), field(explicit, 'final', 12)) == 0
    zero = run('zero_duration', {'thermal_preage': 'on', 'thermal_preage_Ma': 0})
    assert distance(field(default, 'final', 12), field(zero, 'final', 12)) < 1e-14

    zero_legacy = run('zero_bypasses_legacy', {'thermal_preage': 'on', 'thermal_preage_Ma': 0, 'bottom_tbl_thickness': .07848})
    assert distance(field(zero, 'final', 12), field(zero_legacy, 'final', 12)) == 0

    phase = {'phase_depth': '.4,.08,.1', 'phase_delta_s': '.1,0,0',
             'phase_clapeyron': '-.1,0,0', 'phase_width': '.06,.01,.01',
             'phase_transT': '.68,.5,.5'}
    # Single-step independent solve, including a lateral composition gradient.
    for kind, changes, environment in (
        ('constant', {}, {}),
        ('lateral', {'kC_ratio': 3}, {'PREAGE_TEST_COMPOSITION': 'lateral'}),
        ('nonlinear', {'kT_exponent': .5}, {}),
        ('phase', phase, {}),
    ):
        directory = run('reference_' + kind, dict(thermal_preage='on', thermal_preage_Ma=100,
            thermal_preage_max_dt_Ma=100, **changes), env=dict(PREAGE_TEST_QUADRATURE='1', **environment))
        assert results['reference_' + kind]['steps'] == 1, 'reference requires one accepted step'
        exact, observed, report = independent_be(directory, 12, 100, exponent=changes.get('kT_exponent', 0), phase_capacity=kind == 'phase')
        assert report['reference_max_error'] < 3e-9, report
        results['reference_' + kind].update(report)
        if kind == 'lateral':
            radial, _, _ = independent_be(directory, 12, 100, radial_only=True)
            separation = float(np.max(abs(exact - radial)))
            assert separation > 1e-8, ('fixture fails to detect independent-column solver', separation)
            results['reference_lateral']['full3d_vs_radial_only'] = separation

    flavor24 = run('primordial_flavor_24', {'thermal_preage': 'on', 'thermal_preage_Ma': 100,
        'thermal_preage_max_dt_Ma': 100, 'tracer_flavors': 25, 'kC_primordial_flavor': 24,
        'kC_ratio': 3, 'buoyancy_ratio': ','.join(['0']*24), 'z_interface': ','.join(['.75']*24)},
        env={'PREAGE_TEST_COMPOSITION': 'lateral'})
    assert distance(field(root / 'reference_lateral', 'preaged', 12), field(flavor24, 'preaged', 12)) < 1e-13

    for kind, changes, environment in (
        ('constant', {}, {}),
        ('nonlinear', {'kT_exponent': .7}, {}),
        ('lateral', {'kT_exponent': .4, 'kC_ratio': 3}, {'PREAGE_TEST_COMPOSITION': 'lateral'}),
        ('phase', phase, {}),
    ):
        runs = [run(f'{kind}_dt_{dt}', dict(thermal_preage='on', thermal_preage_max_dt_Ma=dt, **changes), env=environment)
                for dt in (250, 125, 62.5)]
        fields = [field(d, 'preaged', 12) for d in runs]
        coarse, fine = distance(fields[0], fields[1]), distance(fields[1], fields[2])
        assert 1e-9 < fine < coarse, (kind, coarse, fine)
        assert 1.4 < coarse / fine < 2.6, ('not first-order timestep convergence', kind, coarse / fine)
        results[kind + '_convergence'] = {'status': 'PASS', 'coarse_difference': coarse,
                                        'fine_difference': fine, 'ratio': coarse / fine}
        decomposed = run(kind + '_radial_24', dict(thermal_preage='on', thermal_preage_max_dt_Ma=125,
            nprocz=2, **changes), ranks=24, env=environment)
        disagreement = distance(fields[1], field(decomposed, 'preaged', 24))
        assert disagreement < 2e-9, ('radial MPI mismatch', kind, disagreement)
        results[kind + '_radial_24']['max_12_rank_difference'] = disagreement

    adaptive = run('phase_adaptation', dict(phase, thermal_preage='on', thermal_preage_Ma=1000,
        thermal_preage_max_dt_Ma=1000, phase_depth='.42,.08,.1', phase_width='.005,.01,.01'))
    assert results['phase_adaptation']['rejected'] > 0, 'large transition-crossing step was not rejected'

    fraction_only = run('phase_fraction_guard', dict(phase, thermal_preage='on', thermal_preage_Ma=1000,
        thermal_preage_max_dt_Ma=1000, phase_depth='.42,.08,.1', phase_width='.005,.01,.01',
        phase_delta_s='0,0,0'))
    assert results['phase_fraction_guard']['rejected'] > 0, 'phase-fraction-only safeguard was not exercised'

    heated = run('source_exclusion', {'thermal_preage': 'on', 'thermal_preage_max_dt_Ma': 125, 'Q0': 30})
    source_error = distance(field(root / 'constant_dt_125', 'preaged', 12), field(heated, 'preaged', 12))
    assert source_error == 0, ('radiogenic heating leaked into pre-age', source_error)
    results['source_exclusion']['max_Q0_difference'] = source_error

    for name, changes in (
        ('negative_duration', {'thermal_preage_Ma': -1}),
        ('zero_dt', {'thermal_preage_max_dt_Ma': 0}),
        ('nan_duration', {'thermal_preage_Ma': 'nan'}),
        ('infinite_dt', {'thermal_preage_max_dt_Ma': 'inf'}),
        ('non_pices', {'energy_solver': 'pg', 'pices_p4': 'off'}),
        ('p5_benchmark', {'p5_case': 'cold'}),
    ):
        run('reject_' + name, dict(thermal_preage='on', **changes), error=True)
    for mode in ('nan', 'absolute_zero'):
        run('reject_initial_' + mode, {'thermal_preage': 'on'}, env={'PREAGE_TEST_BAD_T': mode}, error=True)
    (root / 'summary.json').write_text(json.dumps(results, indent=2) + '\n')
    print(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
