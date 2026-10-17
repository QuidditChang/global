"""Read-only global low-population census; probabilities use cell-step exposures."""
import argparse
import json
from pathlib import Path
import re


def audit(root, start, end):
    rows = {}
    empty = {}
    for line in (root / 'DATA/0/global.log').read_text().splitlines():
        f = dict(re.findall(r'(\w+)=([^\s]+)', line))
        if line.startswith('PICES_COVERAGE_STATS '):
            step = int(f['step'])
            assert step not in rows, ('duplicate', step)
            rows[step] = f
        if line.startswith('PICES_EMPTY_MONITOR '):
            empty[int(f['step'])] = int(f['empty'])
    assert sorted(rows) == list(range(start, end + 1)), 'incomplete census'
    previous = None
    exposures = [0]*10
    became_empty = [0]*10
    risks = {n: [] for n in (2, 4, 8)}
    sequence = []
    for step, f in sorted(rows.items()):
        h = list(map(int, f['histogram'].split(',')))
        t = list(map(int, f['transitions'].split(',')))
        assert len(h) == 10 and len(t) == 100 and min(h+t) >= 0
        assert h[0] == empty[step], ('empty mismatch', step)
        assert int(f['transition_valid']) == int(step > start)
        if previous is not None:
            assert sum(h) == sum(previous), 'cell count changed'
            assert [sum(t[i*10:i*10+10]) for i in range(10)] == previous
            assert [sum(t[i::10]) for i in range(10)] == h
            for i in range(10):
                exposures[i] += previous[i]
                became_empty[i] += t[10*i]
        else:
            assert not any(t), 'restart baseline must not invent transitions'
        for n in risks:
            low = sum(h[:n]); add = sum((n-i)*h[i] for i in range(n))
            assert low == int(f['low%d' % n]) and add == int(f['add%d' % n])
            risks[n].append(dict(step=step, deficient_cells=low, additions=add))
        sequence.append(dict(step=step, histogram=h))
        previous = h
    return dict(status='PASS', start=start, end=end, cells=sum(previous),
                bins=list(range(9))+['>=9'], series=sequence,
                next_step_empty_probability=[dict(previous_count=i if i<9 else '>=9',
                    cell_step_exposures=exposures[i], empty_next_step=became_empty[i],
                    probability=became_empty[i]/exposures[i] if exposures[i] else None)
                    for i in range(10)],
                thresholds={n:dict(peak_additions=max(x['additions'] for x in v),
                    mean_snapshot_additions=sum(x['additions'] for x in v)/len(v),
                    final=v[-1], series=v) for n,v in risks.items()},
                caveat='Hypothetical snapshot deficits, not cumulative insertions or a prediction of a controlled run. Bin zero measures persistence, not new empty events.')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('root', type=Path)
    p.add_argument('--start-step', type=int, default=300)
    p.add_argument('--steps', type=int, default=500)
    a = p.parse_args()
    result = audit(a.root, a.start_step, a.steps)
    (a.root / 'coverage_stats.json').write_text(json.dumps(result, indent=2)+'\n')
    print('PASS: coverage counts, transitions and threshold estimates verified')
