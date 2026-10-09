"""Independent P8b plate, clock, heat, composition and output audit."""
import argparse,gzip,json,math,re
from collections import Counter
from pathlib import Path
from verify_p2 import require,fields,sha
from verify_p4 import boundary

def verify(root,local=False):
 rows=json.loads((root/'input/matrix.json').read_text())['cases']
 require([(r['name'],r['last']) for r in rows]==[('coupled',4),('uncapped',4),('cap_probe',4),('static_plate',4),('present_day',1),('fixed_step',4),('multigrid',1)],'matrix')
 scalet=6371000.**2/(4/(3300*1250))/(1e6*365.25*86400)
 scalev=6371000./(4/(3300*1250))/(100*365.25*86400)
 report={};first_raw={};temperatures={}
 for row in rows:
  name=row['name'];d=root/name;last=row['last'];start=.03 if name=='present_day' else 2.13
  require((d/'case.cfg').read_bytes()==(root/'input'/(name+'.cfg')).read_bytes(),'cfg identity')
  for f in [root/'input/refstate.txt',*(root/'input/forcing').iterdir()]:
   if f.is_file() and not f.name.startswith('._'):require((d/f.relative_to(root/'input')).read_bytes()==f.read_bytes(),'forcing identity')
  require((d/'mpi_exit_code.txt').read_text().strip() in ('0','8'),'MPI status')
  hist=[Counter() for _ in range(last+1)];max_cfl=0;max_balance=0;visc_range=[math.inf,-math.inf];temps=[]
  for rank in range(12):
   folder=d/'DATA'/str(rank);log=(folder/'log').read_text().splitlines()
   require(not any('PICES_ERROR' in line for line in log),'runtime failure')
   records=lambda tag:[fields(line) for line in log if line.startswith(tag+' ')]
   steps=records('PICES_STEP');require([int(x['step']) for x in steps]==list(range(1,last+1)),'accepted steps')
   clock=records('PICES_TIMESTEP');require(len(clock)==(0 if name=='fixed_step' else last),'automatic clock count')
   for i,(s,c) in enumerate(zip(steps,clock if clock else [dict(candidate=1,particle_limit=1)]*last),1):
    dt=float(s['dt']);time=float(s['time']);max_cfl=max(max_cfl,float(s['cfl']))
    require(dt>0 and dt<=float(c['candidate'])*(1+1e-6) and dt<=float(c['particle_limit'])*(1+1e-6),'step bounds')
    expected=.05 if name=='fixed_step' else ([.05,.05,.03,.05][i-1] if name!='present_day' else .03)
    require(abs(dt*scalet-expected)<1e-6,'forcing knot/max step')
    require(float(s['cfl'])<=.24*(1+1e-6),'particle CFL')
   stokes=records('PICES_STOKES')
   if rank==0:require({int(s['step']) for s in stokes}>=set(range(1,last+1)) and all(s['status']=='PASS' for s in stokes),'Stokes convergence')
   comp=records('PICES_COMPOSITION');require(len(comp)==last+1 and all(float(x['error'])<=1e-12 for x in comp),'composition lifecycle')
   projections=records('PICES_PROJECTION')
   if rank==0:require([int(x['step']) for x in projections if x['kind']=='absolute']==list(range(1,last+1)), 'absolute projection count')
   for x in projections:
    require(0<=float(x['residual'])<=float(x['tolerance']) and x['method']=='bounded_consistent_v1','projection')
   heat=records('PICES_EBA');qvis=records('PICES_VISCOUS');ta=records('PICES_TA')
   require(len(heat)==len(qvis)==len(ta)==last,'heat/TA count')
   for h,q,t,s in zip(heat,qvis,ta,steps):
    v=[float(h[k]) for k in ('storage','internal','adiabatic','viscous','phase_pressure','boundary_reaction')]
    err=abs(v[0]-sum(v[1:]));max_balance=max(max_balance,err);require(err<1e-12,'heat ledger')
    raw,applied=float(q['raw']),float(q['applied']);require(0<=applied<=raw*(1+1e-12) and raw>0 and applied==float(h['viscous']),'qvis ledger')
    require(q['state']=='accepted_stokes','source timing')
    if name=='cap_probe':require(applied<raw*1e-3,'inactive cap probe')
    if name=='uncapped':require(applied==raw,'uncapped source')
    require(t['calls']=='1' and float(t['time'])==float(s['time']) and float(t['mapping_error'])<1e-10,'TA arrival clock/mapping')
   if name=='multigrid' and rank==0:
    mg=records('PICES_MG_RESIDUAL');require(len(mg)>0 and all(x['status']=='PASS' and float(x['actual'])<float(x['tolerance']) for x in mg),'true multigrid residual')
   if rank==0:first_raw[name]=float(qvis[0]['raw'])
   for step in range(last+1):
    folder_step=folder/str(step)
    for f in folder_step.glob('*.gz'):
     if f.name.startswith('._'):continue
     lines=gzip.open(f,'rt').read().splitlines()
     require(all(math.isfinite(float(token)) for line in lines for token in line.split()),'finite output')
     if f.name.startswith(('comp_nd.','comp_el.')):
      expected_rows=729 if f.name.startswith('comp_nd.') else 512
      require(len(lines)==expected_rows+1 and int(lines[0].split()[1])==expected_rows,'composition shape')
      require(all(0<=float(v)<=1.000001 for line in lines[1:] for v in line.split()),'composition bounds')
     if f.name.startswith('tracer.'):
      require(len(lines)-1==int(lines[0].split()[1]),'particle count')
      for line in lines[1:]:
       v=list(map(float,line.split()));fl=int(v[3]);require(len(v)==5 and v[3]==fl and 0<=fl<25 and fl not in (18,19),'particle flavor');hist[step][fl]+=1
     if f.name.startswith('visc.'):
      vals=[float(x) for x in lines[1:]];visc_range[0]=min(visc_range[0],min(vals));visc_range[1]=max(visc_range[1],max(vals))
     if f.name.startswith('velo.'):
      data=[list(map(float,x.split())) for x in lines[2:]];require(len(data)==729,'nodal output count')
      age=start if step==0 else start-float(steps[step-1]['time'])*scalet
      expected=0 if name=='static_plate' else .01*(1+age)*scalev
      for v in data[8::9]:require(abs(v[0])<1e-5 and math.isclose(v[1],expected,rel_tol=2e-6,abs_tol=1e-5) and abs(v[2])<1e-5,'plate velocity at accepted age')
      if step==last:temps.extend(v[3] for v in data)
   if name=='present_day':require(abs(start-float(steps[-1]['time'])*scalet)<1e-6,'present-day termination')
  require(all(sum(h.values())==393216 and h[24]==hist[0][24] for h in hist),'particle/primordial conservation')
  require(.09999<=visc_range[0]<visc_range[1]<=100.001,'rheol7 variation/bounds')
  for step in range(1,last+1):
   for side in ('surf','botm'):boundary(d,step,side)
  temperatures[name]=temps
  report[name]=dict(steps=last,max_cfl=max_cfl,max_heat_balance=max_balance,viscosity_range=visc_range,primordial=hist[0][24])
 require(first_raw['coupled']==first_raw['uncapped']==first_raw['cap_probe'],'qvis must not alter initial Stokes viscosity')
 vectors=[]
 for name in ('coupled','multigrid'):
  v=[]
  for rank in range(12):
   lines=gzip.open(root/name/f'DATA/{rank}/0/velo.{rank}.0.gz','rt').read().splitlines()[2:]
   v.extend(float(x) for line in lines for x in line.split()[:3])
  vectors.append(v)
 relative=math.sqrt(sum((a-b)**2 for a,b in zip(*vectors))/sum(a*a for a in vectors[0]))
 require(relative<1e-4,'CG/MG initial velocity consistency')
 report['CG_MG_initial_velocity_relative_L2']=relative
 report['static_plate_temperature_RMS_K']=math.sqrt(sum((a-b)**2 for a,b in zip(temperatures['coupled'],temperatures['static_plate']))/len(temperatures['coupled']))
 require(report['static_plate_temperature_RMS_K']>1e-8,'inactive plate control')
 if not local:
  require((root/'launcher_exit_code.txt').read_text().strip()=='0','launcher')
  require((root/'p8b_complete.txt').read_text().strip()=='P8B_COMPLETE_PENDING_LOCAL_AUDIT','completion')
  require((root/'completed_cases.txt').read_text().strip()=='7','case count')
  require((root/'runtime_source_diff.txt').read_text()=='','dirty runtime')
  commit=(root/'solver_commit.txt').read_text().strip();require(re.fullmatch('[0-9a-f]{40}',commit) and commit==(root/'pices-p0-build.commit').read_text().strip(),'build identity')
  require('P0_BUILD_COMPLETE commit='+commit in (root/'pices-p0-build.log').read_text(),'build completion')
  hashed={}
  for line in (root/'input.sha256').read_text().splitlines():
   digest,file=line.split();require(file not in hashed and file.startswith('input/') and '..' not in Path(file).parts and sha(root/file)==digest,'input hash');hashed[file]=digest
  require(set(hashed)=={str(f.relative_to(root)) for f in (root/'input').rglob('*') if f.is_file()},'input inventory')
  mpi=(root/'mpi_path.txt').read_text().strip();library=(root/'binary_ldd.txt').read_text()
  require('/mpi/intel64/bin/mpiexec.hydra' in mpi and str(Path(mpi).parent.parent/'lib/libmpi.so.12') in library,'MPI provider match')
  require('Intel(R) MPI' in (root/'mpi_version.txt').read_text(),'MPI version')
  for f in ('submitted.lsf','runs_commit.txt','binary.sha256','platform.txt'):require((root/f).stat().st_size>0,'missing provenance '+f)
 return dict(status='PASS',cases=report,production_decision='HOLD_PENDING_P8C_AND_SCIENTIFIC_GATE')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('--local',action='store_true');p.add_argument('--summary',type=Path);a=p.parse_args();r=verify(a.root,a.local);s=json.dumps(r,indent=2)+'\n';print(s)
 if a.summary:a.summary.write_text(s)
