"""Audit a fresh empty1000 run. Requires NumPy for profile checks; no solver-state edits.

Usage: python3 audit_empty_run.py RUN_ROOT [--steps 1000] [--ranks 384]
Output: RUN_ROOT/empty_audit.json. PASS permits scientific review, not automatic
production approval. Neighbor contrast and recovery jumps are proxies, not an
error estimate: neighboring cells may straddle a real compositional interface.
"""
import argparse
from collections import defaultdict, Counter
import itertools
import json
import math
from pathlib import Path
import re

FACES=((0,1,2,3),(4,5,6,7),(0,1,5,4),(1,2,6,5),(2,3,7,6),(3,0,4,7))
TOL=1e-10

def tv(a,b):
    return .5*sum(abs(x-y) for x,y in zip(a,b))

def topology(mesh,meta,targets):
    """Match actual Cartesian face vertices, including independent cap nodes.
    Adjacent spatial bins avoid decimal-rounding discontinuities at seams.
    Unmatched interior faces are missing evidence, never assumed isolated.
    """
    bins=defaultdict(list)
    faces={}
    for key,m in mesh.items():
        points=list(zip(*[iter(m['xyz'])]*3))
        for i,indices in enumerate(FACES):
            verts=[points[j] for j in indices]
            center=tuple(sum(p[d] for p in verts)/4 for d in range(3))
            bucket=tuple(math.floor(x/TOL) for x in center)
            face=(key,i);faces[face]=(verts,bucket)
            bins[bucket].append(face)
    adjacency=defaultdict(set);missing=[]
    for key in targets:
        for i in range(6):
            verts,bucket=faces[key,i]
            matches=[]
            for shift in itertools.product((-1,0,1),repeat=3):
                for other in bins.get(tuple(x+y for x,y in zip(bucket,shift)),()):
                    if other[0]==key:continue
                    ov=faces[other][0]
                    if all(any(math.dist(p,q)<=TOL for q in ov) for p in verts):matches.append(other[0])
            if len(matches)==1:adjacency[key].add(matches[0])
            elif len(matches)>1:missing.append([*key,i,'ambiguous'])
            else:
                radii=[math.sqrt(sum(x*x for x in p)) for p in verts]
                boundary=any(all(abs(r-meta[key[0]][b])<1e-7 for r in radii) for b in ('ri','ro'))
                if not boundary:missing.append([*key,i,'missing'])
    return adjacency,missing

def fields(line):
    return dict(re.findall(r'(\w+)=([^\s]+)',line))

def audit(root,steps,ranks,profile_every=50):
    meta={};mesh={};frames={};cells=defaultdict(dict);issues=[]
    for rank in range(ranks):
        path=root/'PICES_Audit'/('rank_%06d.jsonl'%rank)
        if not path.exists():issues.append('missing '+str(path));continue
        observed=set();times={}
        with path.open() as f:
            for line in f:
                row=json.loads(line);kind=row['type']
                if kind=='meta':
                    if row['rank']!=rank or row['schema']!=1:raise ValueError('audit identity/schema')
                    meta[rank]=row
                elif kind=='frame':
                    s=row['step']
                    if s in observed:raise ValueError('duplicate frame')
                    observed.add(s);times[s]=row['time_Ma']
                elif kind=='mesh':mesh[rank,row['element']]=row
                elif kind=='cell':
                    c=row['c']
                    if not c or any(not math.isfinite(x) or x < -1e-12 or x>1+1e-12 for x in c) or abs(sum(c)-1)>1e-10:
                        issues.append('invalid composition rank=%d step=%d'%(rank,row['step']))
                    cells[row['step']][rank,row['element']]=row
        if observed!=set(range(steps+1)):issues.append('incomplete frames rank=%d'%rank)
        if rank==0:frames=times
        elif times!=frames:issues.append('inconsistent clocks rank=%d'%rank)
    if len(meta)!=ranks:issues.append('incomplete metadata')
    targets={key for rows in cells.values() for key,c in rows.items() if c['count']==0}
    if not targets.issubset(mesh):raise ValueError('missing cell geometry')
    adjacency,missing=topology(mesh,meta,targets)
    if missing:issues.append('unmatched/ambiguous interior faces')
    episodes=[];active={};series=[];max_neighbor=0.;max_cluster=0
    for step in range(steps+1):
        rows=cells[step];empty={k for k,v in rows.items() if v['count']==0}
        for k in set(active)-empty:
            ep=active.pop(k);ep['end']=step;ep['censored']=False
            if k not in rows:issues.append('missing recovery sample %s step=%d'%(k,step))
            else:ep.update(recovery_count=rows[k]['count'],recovery_tv=tv(ep['held'],rows[k]['c']))
            ep['duration_Ma']=frames.get(step,0)-ep['start_time_Ma'];episodes.append(ep)
        contrast=[]
        for k in empty:
            if k not in active:
                xyz=mesh[k]['xyz'];center=[sum(xyz[d::3])/8 for d in range(3)];r=math.sqrt(sum(x*x for x in center))
                active[k]=dict(rank=k[0],element=k[1],start=step,start_time_Ma=frames.get(step,0),held=rows[k]['c'],max_neighbor_tv=0.,
                    lon_deg=math.degrees(math.atan2(center[1],center[0])),lat_deg=math.degrees(math.atan2(center[2],math.hypot(*center[:2]))),
                    depth_km=(meta[k[0]]['ro']-r)*meta[k[0]]['radius_km'])
            ep=active[k]
            if tv(ep['held'],rows[k]['c'])>1e-12:issues.append('empty composition changed %s step=%d'%(k,step))
            for n in adjacency[k]:
                if n not in rows:issues.append('missing halo sample %s step=%d'%(n,step));continue
                if rows[n]['count']:
                    v=tv(rows[k]['c'],rows[n]['c']);contrast.append(v);ep['max_neighbor_tv']=max(ep['max_neighbor_tv'],v)
        todo=set(empty);clusters=[]
        while todo:
            seed=todo.pop();group={seed};queue=[seed]
            while queue:
                for n in adjacency[queue.pop()] & todo:
                    todo.remove(n);group.add(n);queue.append(n)
            clusters.append(len(group))
        max_cluster=max(max_cluster,max(clusters,default=0));max_neighbor=max(max_neighbor,max(contrast,default=0))
        series.append(dict(step=step,empty=len(empty),largest_face_cluster=max(clusters,default=0),max_neighbor_tv=max(contrast,default=0),
                           volume_nd=sum(mesh[k]['volume_nd'] for k in empty)))
    for ep in active.values():
        ep.update(end=steps,censored=True,duration_Ma=frames.get(steps,0)-ep['start_time_Ma']);episodes.append(ep)
    counts={s['step']:s['empty'] for s in series}
    monitor={};diag={};balance={};closure=[];temperatures=[]
    log=root/'DATA/0/global.log'
    if log.exists():
        for line in log.open(errors='replace'):
            if line.startswith('PICES_EMPTY_MONITOR '):
                f=fields(line);monitor[int(f['step'])]=f
            elif line.startswith('PICES_PARTICLE_BALANCE '):
                f=fields(line);balance[int(f['step'])]=f
            elif line.startswith('PICES_DIAGNOSTICS '):
                f=fields(line);diag[int(f['step'])]=f
    if set(monitor)!=set(range(steps+1)):issues.append('incomplete empty monitor')
    if set(balance)!=set(range(steps+1)):issues.append('incomplete particle balances')
    if set(diag)!=set(range(1,steps+1)):issues.append('incomplete accepted-step diagnostics')
    baseline=balance.get(0,{}).get('after')
    for row in series:
        f=monitor.get(row['step'],{})
        row['volume_fraction']=float(f['volume_fraction']) if f else None
        row['longest_current_steps']=int(f['longest_current']) if f else None
    for s,f in monitor.items():
        if counts.get(s)!=int(f['empty']):issues.append('empty count mismatch step=%d'%s)
    for s,f in balance.items():
        if f['after']!=baseline or any(int(f[x]) for x in ('delta','invalid_hosts','count_mismatch_elements')):issues.append('particle balance failure step=%d'%s)
    for s,f in diag.items():
        numeric={k:float(v) for k,v in f.items() if k not in ('stage','scope','energy','source')}
        if any(not math.isfinite(v) for v in numeric.values()):issues.append('nonfinite diagnostics step=%d'%s)
        if str(int(numeric['particle_count_global']))!=baseline:issues.append('accepted particle count mismatch step=%d'%s)
        value=float(f['closure_residual_nd']);closure.append(abs(value))
        # Roundoff budget relative to the frozen-capacity energy scale.
        scale=max(1.,abs(float(f['energy_before_nd'])),abs(float(f['energy_after_nd'])))
        if not math.isfinite(value) or abs(value)>1e-10*scale:issues.append('energy closure failure step=%d'%s)
        for name in ('particle_Tmin_K','grid_Tmin_K'):
            t=float(f[name]);temperatures.append(t)
            if not math.isfinite(t) or t < -1e-9:issues.append('nonphysical temperature step=%d'%s)
    import numpy as np
    profile_checks=0
    for step in range(profile_every,steps+1,profile_every):
        path=root/'DATA/0'/('profiles_hist_%d.npz'%step)
        if not path.exists():issues.append('missing profile step=%d'%step);continue
        with np.load(path,allow_pickle=False) as z:
            if int(z['pices_diagnostics_valid'])!=1 or int(z['pices_diagnostics_step'])!=step:
                issues.append('invalid profile stage step=%d'%step);continue
            for key in z.files:
                name=key.removeprefix('pices_') if hasattr(str,'removeprefix') else key[6:]
                if key.startswith('pices_') and name in diag.get(step,{}):
                    if float(z[key])!=float(diag[step][name]):issues.append('NPZ/log mismatch step=%d field=%s'%(step,name))
                    profile_checks+=1
            if profile_checks < 25*(step//profile_every):issues.append('missing profile diagnostic fields step=%d'%step)
    for name in ('solver.finished.json','result.json'):
        p=root/'PostProc/Online'/name
        if not p.exists():issues.append('missing '+name);continue
        r=json.loads(p.read_text())
        if name=='solver.finished.json' and r.get('returncode')!=0:issues.append('solver did not finish successfully')
        if name=='result.json' and (r.get('solver_returncode')!=0 or r.get('postproc_returncode')!=0 or r.get('drain_timeout')):issues.append('postprocessing incomplete/failed')
    return dict(status='REVIEW_REQUIRED' if not issues else 'INCOMPLETE_OR_FAILED',issues=sorted(set(issues)),
        scope='Sparse global face-neighbor audit; compositions at post-advection stage. TV is contrast, not numerical error.',
        decision='Do not auto-release from this report. Review growth, persistent clusters, and composition contrast/recovery jumps; hard failures require correction, not automatically particle addition.',
        profile_scalar_comparisons=profile_checks,max_face_cluster=max_cluster,max_neighbor_tv=max_neighbor,episode_count=len(episodes),ongoing=len(active),
        max_abs_energy_closure_nd=max(closure,default=None),minimum_temperature_K=min(temperatures,default=None),
        unmatched_faces=missing,series=series,episodes=episodes)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path);p.add_argument('--steps',type=int,default=1000);p.add_argument('--ranks',type=int,default=384);p.add_argument('--profile-every',type=int,default=50);a=p.parse_args()
    result=audit(a.root,a.steps,a.ranks,a.profile_every);out=a.root/'empty_audit.json';out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('series','episodes','unmatched_faces')},indent=2));print(out)
    return 0 if not result['issues'] else 1
if __name__=='__main__':raise SystemExit(main())
