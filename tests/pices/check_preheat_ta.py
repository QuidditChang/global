"""Check manufactured TA runs: ordering, closure, radial BCs and restart.

Usage: python3 check_preheat_ta.py BEFORE [--restart RESTART]
Run PICESManufactured with rotation, lumped, pices_ta_before_heat=on,
24 ranks (nprocz=2), using the pices_p8c fixture for four steps.
"""
import argparse
import json
import math
from pathlib import Path


def fields(line):
    return dict(word.split('=', 1) for word in line.split() if '=' in word)


def check(root):
    text = (root / 'DATA/0/log').read_text()
    lines = text.splitlines()
    rows = [fields(x) for x in lines if x.startswith('PICES_DIAGNOSTICS ')]
    assert rows and int(rows[-1]['step']) == 4
    assert 'PICES_ERROR' not in (root / 'solver.stderr').read_text()
    closure = max(abs(float(x['closure_residual_nd'])) for x in rows)
    assert closure < 1e-12, closure
    assert all(float(x['particle_Tmin_K']) > 0 for x in rows)
    assert len({x['particle_count_global'] for x in rows}) == 1
    for row in rows:
        step = row['step']
        events = [x for x in lines if fields(x).get('step') == step]
        ta = [i for i, x in enumerate(events) if x.startswith('PICES_TA ')]
        assert len(ta) == 1 and 'stage=before_heat' in events[ta[0]]
        stages = [(i, fields(x)) for i, x in enumerate(events)
                  if x.startswith('PICES_STAGE ') and 'field=grid' in x]
        assert next(i for i, x in stages if x['stage'] == 'after_remap') < ta[0]
        assert ta[0] < next(i for i, x in stages if x['stage'] == 'heat_before_particle')
    # Physical radial Dirichlet values must survive an in-place TA projection.
    count = 0
    for path in root.glob('nodes.*.4.txt'):
        for line in path.read_text().splitlines():
            v = list(map(float, line.split()))
            if not v[5]:
                continue
            radius = math.sqrt(sum(x*x for x in v[1:4]))
            expected = 0.0 if radius > .9 else 1.0
            assert abs(v[4] - expected) < 1e-12, (path, v)
            count += 1
    assert count > 0
    return {'closure_max': closure, 'fixed_boundary_samples': count,
            'last_step': 4, 'particle_count': int(float(rows[-1]['particle_count_global']))}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('--restart', type=Path)
    args = parser.parse_args()
    result = {'continuous': check(args.root)}
    if args.restart:
        result['restart'] = check(args.restart)
        files = sorted(args.root.glob('DATA/*/*.chkpt.*.4.pices.state'))
        assert len(files) == 24
        for path in files:
            assert path.read_bytes() == (args.restart / path.relative_to(args.root)).read_bytes(), path
        result['restart_particle_states_byte_identical'] = len(files)
    print(json.dumps(result, indent=2))
