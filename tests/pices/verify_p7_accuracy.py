"""Audit P7 accuracy follow-up without treating PG or a finer grid as truth."""
import argparse
import gzip
import json
import math
from pathlib import Path
from verify_p5 import verify as integrity, temps, difference
from verify_p2 import require, fields, config


def name(method, n, particles=64, steps=32):
    return f'assim_{method}_n{n}_p{particles if method == "pices" else 0}_s{steps}'


def common_nodes(values, n):
    stride = (n-1)//4
    return [values[rank*n**3 + y*n*n + x*n + z]
            for rank in range(12) for y in range(0,n,stride)
            for x in range(0,n,stride) for z in range(0,n,stride)]


def verify(root, local=False):
    root = Path(root)
    rows = json.loads((root/'input/matrix.json').read_text())['cases']
    expected = {name(m,n) for n in (5,9,17) for m in ('pg','pices')}
    expected |= {name('pices',9,p) for p in (32,128)}
    expected |= {name(m,17,steps=64) for m in ('pg','pices')}
    require({r['name'] for r in rows} == expected, 'accuracy matrix identity')
    for row in rows:
        require(row['nodes'] in (5,9,17) and row['steps'] in (32,64), 'grid/time design')
        require(row['dt'] == 4e-7/row['steps'], 'matched final time')
        require(row['name'] == name(row['method'],row['nodes'],row['particles'],row['steps']), 'case identity')
        c = config(root/'input/cases'/f'{row["name"]}.cfg')
        require(c['p5_length_scale']=='1', 'length scale')
        if row['method']=='pices':
            require(c.get('pices_projection')=='bounded_consistent', 'PICES projection')
    report = integrity(root, local=local, stage='P7_accuracy')
    final = {r['name']:temps(root/r['name'],r['steps'],r['nodes']) for r in rows}
    initial = {r['name']:temps(root/r['name'],0,r['nodes']) for r in rows}
    pairs = {}
    for n in (5,9,17):
        pg, pic = name('pg',n), name('pices',n)
        require(initial[pg] == initial[pic], 'paired initial fields')
        delta = [a-b for a,b in zip(final[pic],final[pg])]
        pairs[str(n)] = dict(nodal_rms_K=difference(final[pic],final[pg]),
            radial_rms_K=[math.sqrt(sum(v*v for v in delta[z::n])/len(delta[z::n])) for z in range(n)],
            vrms_ratio_minus_one_percent=100*(report['cases'][pic]['final']['vrms']/report['cases'][pg]['final']['vrms']-1))
    # Common logical nodes must also coincide in physical spherical coordinates.
    coordinates = []
    for n in (5,9,17):
        xyz = []
        for rank in range(12):
            with gzip.open(root/name('pg',n)/f'DATA/{rank}/coord.{rank}.gz','rt') as f:
                next(f)
                xyz.extend(tuple(map(float,line.split())) for line in f)
        require(len(xyz)==12*n**3 and all(len(x)==3 for x in xyz), 'coordinate shape')
        coordinates.append(common_nodes(xyz,n))
    coordinate_gap = [max(abs(a-b) for p,q in zip(xyz,coordinates[2]) for a,b in zip(p,q))
                      for xyz in coordinates[:2]]
    aligned = max(coordinate_gap) < 2e-6
    spatial = dict(status='ALIGNED' if aligned else 'NON_NESTED_NOT_COMPARABLE',
                   maximum_spherical_coordinate_difference=coordinate_gap,
                   note='Equal logical indices are not necessarily equal physical locations. '
                        'Do not interpret index-wise differences as spatial convergence.', methods={})
    if aligned:
        for method in ('pg','pices'):
            values = [common_nodes(final[name(method,n)],n) for n in (5,9,17)]
            starts = [common_nodes(initial[name(method,n)],n) for n in (5,9,17)]
            spatial['methods'][method] = dict(common_nodes_per_cap=125,
                initial_rms_K=[difference(starts[i],starts[2]) for i in (0,1)],
                successive_final_rms_K=[difference(values[i],values[i+1]) for i in (0,1)],
                note='Sampled nodes, not a volume-weighted norm; boundary copies retained.')
    particle = {f'{a}_to_{b}':difference(final[name('pices',9,a)],final[name('pices',9,b)]) for a,b in ((32,64),(64,128))}
    temporal = {m:difference(final[name(m,17)],final[name(m,17,steps=64)]) for m in ('pg','pices')}
    pairs['17_half_dt'] = dict(nodal_rms_K=difference(final[name('pices',17,steps=64)],final[name('pg',17,steps=64)]),
        vrms_ratio_minus_one_percent=100*(report['cases'][name('pices',17,steps=64)]['final']['vrms']/report['cases'][name('pg',17,steps=64)]['final']['vrms']-1))
    for row in rows:
        if row['method'] != 'pices': continue
        d = root/row['name']
        for rank in range(12):
            log=(d/f'DATA/{rank}/log').read_text().splitlines()
            coverage=[fields(l) for l in log if l.startswith('PICES_COVERAGE ')]
            require([int(x['step']) for x in coverage]==list(range(1,row['steps']+1)), 'coverage completeness')
            require(all(int(x['empty_elements'])==0 for x in coverage), 'empty element')
        projections=[fields(l) for l in (d/'DATA/0/log').read_text().splitlines() if l.startswith('PICES_PROJECTION ')]
        for s in range(1,row['steps']+1):
            part=[x for x in projections if int(x['step'])==s]
            require(sum(x['kind']=='absolute' for x in part)==1 and sum(x['kind']=='signed' for x in part)>=3, 'projection lifecycle')
        require(all(x['method']=='bounded_consistent_v1' and 0<=float(x['residual'])<=float(x['tolerance']) for x in projections), 'projection residual')
    report.update(paired_grids=pairs, spatial=spatial, particle_successive_rms_K=particle,
                  finest_grid_time_rms_K=temporal,
                  production_decision='REVIEW_REQUIRED: inspect spatial and particle convergence against measured time sensitivity; no automatic production release')
    return report


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path);p.add_argument('--local',action='store_true');p.add_argument('--summary',type=Path)
    a=p.parse_args()
    try: result=verify(a.root,a.local)
    except (AssertionError,ValueError,OSError,KeyError) as e: result=dict(status='FAIL',error=str(e))
    if a.summary:a.summary.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k in ('status','error','paired_grids','spatial','particle_successive_rms_K','finest_grid_time_rms_K','production_decision')},indent=2))
    raise SystemExit(result['status']=='FAIL')
