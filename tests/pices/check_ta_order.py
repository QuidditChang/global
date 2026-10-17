"""Verify the production runtime configuration and each accepted TA stage."""
import argparse
from pathlib import Path
import re


def check(root, steps, start_step=0):
    configs = list(root.glob('pid*.cfg'))
    assert configs, 'missing generated production configuration'
    for cfg in configs:
        assert re.search(r'^pices_ta_before_heat\s*=\s*1\s*$', cfg.read_text(), re.M), cfg
    lines = (root / 'DATA/0/global.log').read_text().splitlines()
    for step in range(start_step + 1, steps + 1):
        rows = [x for x in lines if re.search(r'\bstep=%d\b' % step, x)]
        ta = [i for i, x in enumerate(rows) if x.startswith('PICES_TA ')]
        assert len(ta) == 1 and 'stage=before_heat' in rows[ta[0]], ('TA stage', step)
        heat = [i for i, x in enumerate(rows) if x.startswith('PICES_STAGE ') and 'stage=heat_before_particle ' in x]
        assert heat and ta[0] < min(heat), ('heat order', step)
    print('PASS: generated flag and before-heat TA verified for %d steps' % (steps-start_step))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('root', type=Path)
    p.add_argument('--steps', type=int, default=300)
    p.add_argument('--start-step', type=int, default=0)
    a = p.parse_args()
    check(a.root, a.steps, a.start_step)
