"""P1 integration tests; requires an existing build and local 12-rank MPI.

Uses real production C operators with a test-only prescribed-field driver.
The independent NumPy reference assembles a global FE system and applies the
published P/Q/subgrid formulas, including all heat substeps.
"""
import argparse,gzip,json,math,os,re,subprocess,tempfile
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('build',type=Path);p.add_argument('cfg',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
a.output.mkdir(parents=True,exist_ok=True)
base=a.cfg.read_text();ref=a.cfg.parent/'refstate_EBA_PICES_P1.txt'
results={}
def run(name,mode,changes=None,ok=True,token=None):
 d=a.output/name;d.mkdir();(d/'DATA').mkdir();cfg=base
 changes=dict(changes or {});changes.setdefault('tracers_per_element','32')
 for k,v in changes.items():
  pattern=r'^'+re.escape(k)+r'=.*$'
  cfg=re.sub(pattern,k+'='+str(v),cfg,flags=re.M) if re.search(pattern,cfg,re.M) else cfg+'\n'+k+'='+str(v)+'\n'
 (d/'case.cfg').write_text(cfg);(d/ref.name).write_bytes(ref.read_bytes())
 with (d/'stdout').open('w') as out,(d/'stderr').open('w') as err:
  r=subprocess.run([os.environ.get('MPIEXEC','/usr/local/bin/mpiexec'),'--oversubscribe','-n','12',str(a.build.resolve()/'PICESManufactured'),'case.cfg'],cwd=d,env=dict(os.environ,PICES_TEST_MODE=mode),stdout=out,stderr=err,timeout=90)
 if ok: assert r.returncode in (0,8),(name,r.returncode,(d/'stderr').read_text()[-1800:])
 else:
  assert r.returncode not in (0,8),(name,r.returncode)
  assert 'PICES_ERROR' in (d/'stderr').read_text() and token in (d/'stderr').read_text(),(name,'missing failure reason')
  assert not any('PICES_STEP' in f.read_text() for f in (d/'DATA').glob('*/log'))
 results[name]={'exit':r.returncode,'status':'PASS'}
 return d

def records(d,rank):
 return [dict(re.findall(r'(\w+)=([^\s]+)',s)) for s in (d/'DATA'/str(rank)/'log').read_text().splitlines() if s.startswith('PICES_STEP ')]

# Constant preservation with genuine multi-substep diffusion and same MPI maps.
d=run('constant','constant',{'fixed_timestep':'.01'})
for rank in range(12):
 for step in (0,1,2):
  rows=np.loadtxt(d/f'nodes.{rank}.{step}.txt');assert np.max(abs(rows[:,4]-.5))<1e-7
  parts=np.loadtxt(d/f'particles.{rank}.{step}.txt');assert np.max(abs(parts[:,0]-.5))<2e-13
 assert all(int(r['substeps'])>1 for r in records(d,rank))

# Advection alone must preserve the multiset of persistent Tp through migration.
d=run('rotation','rotation',{'fixed_timestep':'.005','pices_test_no_diffusion':'on'})
start=np.sort(np.concatenate([np.loadtxt(d/f'particles.{r}.0.txt')[:,0] for r in range(12)]))
end=np.sort(np.concatenate([np.loadtxt(d/f'particles.{r}.2.txt')[:,0] for r in range(12)]))
assert start.shape==end.shape and np.array_equal(start,end)
assert all(float(x['subgrid_max'])==0 for r in range(12) for x in records(d,r))

# Diffusion with zero velocity: independent global matrix and particle update.
d=run('diffusion','diffusion',{'fixed_timestep':'.01','maxstep':'1','maxtotstep':'2'})
lookup={};maps=[];nodedata=[]
for rank in range(12):
 rows=np.loadtxt(d/f'nodes.{rank}.0.txt');ids=[]
 for row in rows:
  key=tuple(np.round(row[1:4],10))
  if key not in lookup:lookup[key]=len(lookup);nodedata.append(row)
  ids.append(lookup[key])
 maps.append(np.array(ids))
nodedata=np.array(nodedata);N=len(lookup);K=np.zeros((N,N));parts=[];weights=[];pids=[];lengths=[]
for rank in range(12):
 elem=np.loadtxt(d/f'elements.{rank}.txt');km=np.loadtxt(d/f'matrix.{rank}.txt').reshape(-1,8,8)
 for row,ke in zip(elem,km):
  assert np.max(abs(ke-ke.T))<1e-13
  assert abs(ke.sum(axis=1)).max()<1e-12
  assert np.linalg.eigvalsh(ke).min()>-1e-12
  ids=maps[rank][row[:8].astype(int)-1];K[np.ix_(ids,ids)]+=ke
 rows=np.loadtxt(d/f'particles.{rank}.0.txt');ids=maps[rank][rows[:,1:17:2].astype(int)-1];ww=rows[:,2:17:2]
 assert (ww>=0).all() and np.max(abs(ww.sum(axis=1)-1))<1e-12
 # Recover host by node tuple (all eight node ids are emitted even at zero weights).
 le={tuple(row[:8].astype(int)):row[8] for row in elem}
 lengths.extend(le[tuple(row[1:17:2].astype(int))] for row in rows)
 parts.extend(rows[:,0]);weights.extend(ww);pids.extend(ids)
parts=np.array(parts);weights=np.array(weights);pids=np.array(pids);lengths=np.array(lengths)
den=np.bincount(pids.ravel(),weights=weights.ravel(),minlength=N)
assert (den>0).all();fixed=nodedata[:,5].astype(bool);mass=nodedata[:,6];assert (mass>0).all()
free=~fixed; sym=K[np.ix_(free,free)]/np.sqrt(mass[free,None]*mass[None,free]);eigen=np.linalg.eigvalsh(sym);assert eigen.min()>-1e-11
heat_dt=float(re.search(r'dt_heat=([^ ]+)',(d/'DATA/0/log').read_text()).group(1));assert heat_dt*eigen.max()<=.8+1e-12
def P(v):return np.bincount(pids.ravel(),weights=(weights*v[:,None]).ravel(),minlength=N)/den
def Q(v):return (weights*v[pids]).sum(axis=1)
g=P(parts);g[fixed]=.5
row=records(d,0)[0];dt=float(row['dt']);ns=int(row['substeps']);assert ns>1
raw_survival=0
for i in range(ns):
 ds=dt/ns;old=g.copy();rhs=-(K@g);g[~fixed]+=ds*rhs[~fixed]/mass[~fixed]
 sub=(Q(g)-parts)*(-np.expm1(-ds/lengths**2))
 raw_survival=max(raw_survival,float(np.max(abs(sub-Q(P(sub))))))
 parts+=sub+Q(g-old-P(sub))
actual=np.concatenate([np.loadtxt(d/f'particles.{r}.1.txt')[:,0] for r in range(12)])
assert np.max(abs(parts-actual))<3e-12,np.max(abs(parts-actual))
for rank in range(12):
 actual_g=np.loadtxt(d/f'nodes.{rank}.1.txt')[:,4]
 assert np.max(abs(actual_g-g[maps[rank]]))<6e-8
assert raw_survival>1e-7,raw_survival
results['diffusion'].update(particle_reference_max_error=float(np.max(abs(parts-actual))),raw_subgrid_survival=raw_survival,heat_substeps=ns)

d=run('flavor_slot','constant',{'tracer_flavors':'1'})
for rank in range(12):
 assert 'Tp_slot=1' in (d/'DATA'/str(rank)/'log').read_text()
 assert np.max(abs(np.loadtxt(d/f'particles.{rank}.2.txt')[:,0]-.5))<2e-13
run('reject_empty','empty',{},False,'zero particle coverage')
run('reject_source','constant',{'Q0':'.1'},False,'requires Q0=0')
run('reject_cfl','rotation',{'fixed_timestep':'1'},False,'particle CFL')
run('reject_substeps','diffusion',{'fixed_timestep':'.01','pices_max_substeps':'1'},False,'substep limit')
run('reject_restart','constant',{'restart':'on'},False,'restart')
(a.output/'summary.json').write_text(json.dumps(results,indent=2)+'\n')
print(json.dumps(results,indent=2))
