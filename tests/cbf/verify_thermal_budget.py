"""Verify MPI thermal-budget rows against native boundary files and source sums.

python tests/cbf/verify_thermal_budget.py RUN_DIR --steps 0 1 2 --constant-q0 .3
For constant-q0 validation the fixture must have rho_ref=1 everywhere.
"""
import argparse
import math
from pathlib import Path
import re
import struct


def close(a, b):
    assert math.isclose(a, b, rel_tol=2e-12, abs_tol=1e-12), (a, b)


def verify(root, steps, constant_q0=None):
    log = root/'DATA/0/log'
    if not log.exists():
        log, = (root/'DATA/0').glob('*.log')
    if constant_q0 is not None:
        constant_q0=struct.unpack('f',struct.pack('f',constant_q0))[0]
    records = {}
    current = None
    for line in log.read_text().splitlines():
        match = re.match(r'THERMAL_BUDGET  step=(\d+)', line)
        if match:
            step = int(match[1])
            assert step not in records, ('duplicate budget', step)
            current = records[step] = {}
        elif line.startswith('NET_SOURCE_PLUS_BOUNDARY_W='):
            current = None
        elif current is None:
            continue
        elif line.startswith('BUDGET_SCALES '):
            current['scales'] = dict(re.findall(r'(\w+)=([^\s]+)', line))
        elif re.match(r'^(Qtotal|Qvisc(?:_raw|_capped|_removed|_potential_removed|-Qadi_base)?|Qadi(?:_base)?|Qphase|Qinternal|Qassim|q_surf|q_botm)\s+[+-]', line):
            label, *values = line.split()
            current[label] = list(map(float, values))
    assert set(records) == set(steps), records.keys()
    for step, rows in records.items():
        for label, values in rows.items():
            if label != 'scales':
                assert all(map(math.isfinite, values)), (step,label,values)
                assert values[1] >= values[2], (step,label,values)
        close(rows['Qtotal'][0], sum(sign*rows[label][0] for label,sign in
              [('Qinternal',1),('Qvisc',1),('Qadi',-1),('Qphase',-1),('Qassim',1)]))
        if constant_q0 is not None:
            scales = rows['scales']
            close(rows['Qinternal'][0],constant_q0*float(scales['volume_nd'])*float(scales['power_W']))
            for value in rows['Qinternal'][1:]:
                close(value,constant_q0*float(scales['density_W_m3']))
        for side in ('surf','botm'):
            files = sorted(root.glob('DATA/*/q.%s.*.%d' % (side,step)))
            if not files:
                print('step=%d q_%s logged without native files PASS' % (step,side))
                continue
            values = []; heat = []; integrated=0.
            for path in files:
                nodes={};faces=[];meta={}
                for line in path.read_text().splitlines():
                    if line.startswith('#'):
                        meta.update(re.findall(r'(\w+)=([^\s]+)',line))
                    elif line.startswith('N '):
                        fields=line.split();q=float(fields[-3]);values.append(q)
                        nodes[(int(fields[1]),int(fields[2]))]=q
                    elif line.startswith('F '):
                        faces.append(line.split())
                heat.append(float(meta['global_heat_W']))
                length=float(meta['length_scale_m'])
                for f in faces:
                    integrated += sum(nodes[(int(f[1]),int(n))]*float(w)*length**2
                                      for n,w in zip(f[3:7],f[7:11]))
            row=rows['q_'+side]
            for total in heat:
                close(row[0],total)
            close(row[0],integrated);close(row[1],max(values));close(row[2],min(values))
            print('step=%d q_%s native integral/extrema match PASS' % (step,side))
        if 'Qvisc_removed' in rows:
            close(rows['Qvisc_removed'][0],rows['Qvisc_raw'][0]-rows['Qvisc'][0])
            close(rows['Qvisc_potential_removed'][0],rows['Qvisc_raw'][0]-rows['Qvisc_capped'][0])
    print('All source totals, units and budget rows PASS')


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    p.add_argument('--steps',type=int,nargs='+',required=True)
    p.add_argument('--constant-q0',type=float)
    args=p.parse_args();verify(args.root,args.steps,args.constant_q0)
