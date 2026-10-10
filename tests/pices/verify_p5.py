"""P5 integrity and paired diagnostics. Completion is not a production approval."""
import argparse,gzip,json,math,re,hashlib,subprocess
from pathlib import Path
from verify_p2 import require,fields,config
from verify_p4 import boundary

def temps(d,step,n):
 values=[]
 for rank in range(12):
  p=d/f'DATA/{rank}/{step}/velo.{rank}.{step}.gz'
  lines=gzip.decompress(p.read_bytes()).decode().splitlines();require(len(lines)==n**3+2,'grid length '+str(p))
  for i,l in enumerate(lines[2:],1):
   v=list(map(float,l.split()));require(len(v)==4 and all(map(math.isfinite,v)),'nonfinite field');require(v[3]>=0,'negative absolute temperature')
   if i%n==0:require(v[3]==300,'top boundary')
   if i%n==1:require(v[3]==3700,'bottom boundary')
   values.append(v[3])
 return values

def difference(a,b):return math.sqrt(sum((x-y)**2 for x,y in zip(a,b))/len(a))
def verify(root,local=False,partial=False,stage="P5"):
 root=Path(root);matrix=json.loads((root/'input/matrix.json').read_text());rows=matrix['cases'];results={};finals={};initials={};complete=[]
 expected={'P5':46,'P7':11,'P7_accuracy':10}.get(stage)
 require(stage in ('P5','P7','P7_accuracy') and len(rows)==expected and len({r['name'] for r in rows})==expected,'matrix size/names')
 for row in rows:
  name=row['name'];d=root/name
  if partial and not (d/'mpi_exit_code.txt').exists():continue
  require((d/'case.cfg').read_bytes()==(root/f'input/cases/{name}.cfg').read_bytes(),'actual cfg')
  which='assim' if row['scenario']=='assim' else 'transport';n=row['nodes'];steps=row['steps'];pic=row['method']=='pices'
  cfg=config(d/'case.cfg')
  for k,v in dict(energy_solver=row['method'],p5_case=row['scenario'],p5_prescribed_velocity='on' if row['prescribed'] else 'off',nodex=str(n),nodey=str(n),nodez=str(n),tracer=str(int(pic)),pices_checkpoint='off',filter_temp='off',maxstep=str(steps),restart='off').items():require(cfg.get(k)==v,'configuration '+name+' '+k)
  require(float(cfg['fixed_timestep'])==row['dt'] and float(cfg['p5_length_scale'])==row['length_scale'],'dt/length manifest')
  require(cfg['CBF_use_advection']==('off' if pic else 'on'),'CBF derivative policy')
  require((d/'refstate.txt').read_bytes()==(root/f'input/refstate_{which}_{n}.txt').read_bytes(),'actual refstate')
  for f in (root/f'input/forcing_{n}').iterdir():require((d/'forcing'/f.name).read_bytes()==f.read_bytes(),'actual forcing')
  require(int((d/'mpi_exit_code.txt').read_text()) in (0,8),'MPI exit '+name)
  err=(d/'solver.stderr').read_text();require(not re.search(r'PICES_ERROR|MPI_ABORT|Segmentation fault|No such file|Cannot open',err,re.I),'runtime error '+name)
  log=(d/'DATA/0/log').read_text();metric=[{k:float(v) for k,v in fields(l).items()} for l in log.splitlines() if l.startswith('P5_METRIC ')]
  require([int(v['step']) for v in metric]==list(range(steps+1)),'step completeness '+name)
  require(all(all(map(math.isfinite,v.values())) for v in metric),'finite metrics')
  require(all(math.isclose(v['time'],v['step']*row['dt'],rel_tol=2e-6,abs_tol=1e-14) for v in metric),'physical time')
  initial,final=metric[0],metric[-1]
  require(all(v['rss_max_KiB']>0 and v['volume']>0 for v in metric),'diagnostic geometry/memory')
  stokes=[]
  if not row['prescribed']:
   require(cfg.get('vlowstep')==('2000' if stage=='P7_accuracy' or (stage=='P7' and row['variant']=='rheol7') else '1000'),'coupled inner iteration budget')
   stokes=[fields(l) for l in log.splitlines() if l.startswith('PICES_STOKES ')]
   require([int(v['step']) for v in stokes]==list(range(steps+1)),'shared Stokes guard '+name)
   require(all(v['status']=='PASS' and float(v['inner_relative'])==1e-6 for v in stokes),'Stokes failure')
  heat={int(fields(l)['step']):fields(l) for l in log.splitlines() if l.startswith('PICES_EBA ')}
  adv={int(fields(l)['step']):fields(l) for l in log.splitlines() if l.startswith('PICES_STEP ')}
  ta={int(fields(l)['step']):fields(l) for l in log.splitlines() if l.startswith('PICES_TA ')}
  budget=None
  if pic:
   expected=12*(n-1)**3*row['particles'];require(all(v['particles']==expected for v in metric),'particle count')
   require(set(heat)==set(adv)==set(range(1,steps+1)),'PIC heat logs')
   require(all(0<=float(v['cfl'])<=.25 for v in adv.values()),'particle CFL')
   for rank in range(12):
    ls=(d/f'DATA/{rank}/log').read_text().splitlines()
    t=[fields(l) for l in ls if l.startswith('PICES_TA ')]
    require([int(v['step']) for v in t]==(list(range(1,steps+1)) if which=='assim' else []),'TA lifecycle')
    require(all(v['calls']=='1' for v in t),'TA duplicate')
   increments=sum(float(adv[i]['remap_energy'])+float(heat[i]['storage'])+(float(ta[i]['storage']) if i in ta else 0) for i in heat)
   budget=final['energy_proxy']-initial['energy_proxy']-increments
   # P5 has no phase entropy, so rho*Cp is a static heat capacity here.
   require(abs(budget)<1e-10*max(1,abs(initial['energy_proxy'])),'grid energy/remap/TA closure '+name)
  flux={side:boundary(d,steps,side,'pices' if pic else 'pg') for side in ['surf','botm']}
  if pic:
   reaction=float(heat[steps]['boundary_reaction'])/float(adv[steps]['dt'])*4*3400*6371000
   # Sharp plateaus have zero boundary conduction up to roundoff; use a
   # 1e-12 nondimensional absolute flux floor only for that explicit test.
   tolerance=1e-10*max(1,abs(reaction)+abs(flux['botm'])+abs(flux['surf']))
   if row['scenario']=='sharp':tolerance=max(tolerance,1e-12*4*3400*6371000)
   require(abs(flux['botm']-flux['surf']-reaction)<tolerance,'CBF heat boundary ledger')
   for step in [0,steps]:
    count=0
    for rank in range(12):
     with gzip.open(d/f'DATA/{rank}/{step}/tracer.{rank}.{step}.gz','rt') as f:
      header=f.readline().split();np=int(header[1]);require(int(header[0])==step and int(header[2])==4,'tracer header');num=0
      for line in f:
       v=list(map(float,line.split()));require(len(v)==4 and all(map(math.isfinite,v)) and 300+3400*v[3]>=0,'particle temperature');num+=1
      require(num==np,'tracer file count');count+=np
    require(count==expected,'global tracer count')
  # Verify initial and final primary fields independently of diagnostic text.
  initials[name]=temps(d,0,n);finals[name]=temps(d,steps,n)
  initial_error=0.
  if which=='transport':
   for rank in range(12):
    xyz=gzip.decompress((d/f'DATA/{rank}/coord.{rank}.gz').read_bytes()).decode().splitlines()[1:]
    require(len(xyz)==n**3,'coordinate count')
    for i,line in enumerate(xyz):
     th,ph,rad=map(float,line.split());center=.82 if row['scenario']=='cold' else .68
     x=rad*math.sin(th)*math.cos(ph);y=rad*math.sin(th)*math.sin(ph);z=rad*math.cos(th)
     anomaly=.15*math.sin(math.pi*(rad-.55)/.45)**2*math.exp(-((x-center)**2+y*y+z*z)/.12**2)
     expected_T=300+3400*((1-rad)/.45+(-1 if row['scenario']=='cold' else 1)*anomaly)
     if row['scenario']=='sharp':expected_T=3700 if rad<.77 else 300
     initial_error=max(initial_error,abs(initials[name][rank*n**3+i]-expected_T))
   require(initial_error<.02,'analytic initial field '+name) # printed coordinates/T have finite decimal precision

  results[name]=dict(initial=initial,final=final,CBF_final_W=flux,initial_analytic_error_K=initial_error,energy_proxy_change=final['energy_proxy']-initial['energy_proxy'],grid_ledger_residual=budget,
    stokes_iterations=[int(v['iterations']) for v in stokes],wall_seconds=float((d/'wall_seconds.txt').read_text()),
    particle_diagnostics=dict(max_cfl=max(float(v['cfl']) for v in adv.values()),
      remap_energy_by_step=[float(adv[i]['remap_energy']) for i in range(1,steps+1)],
      first_remap_fraction_of_initial_energy=float(adv[1]['remap_energy'])/initial['energy_proxy'],
      heat_storage_total=sum(float(v['storage']) for v in heat.values()),
      TA_storage_total=sum(float(v['storage']) for v in ta.values())) if pic else None)
  complete.append(name)
 pairs={}
 for row in rows:
  name=row['name']
  if row['method']!='pices' or name not in results:continue
  pg=name.replace('_pices_','_pg_')
  if pg not in results:continue
  require(math.isclose(results[name]['initial']['energy_proxy'],results[pg]['initial']['energy_proxy'],rel_tol=1e-13),'initial integral mismatch')
  require(initials[name]==initials[pg],'PG/PICES initial temperatures differ '+name)
  pairs[name]=dict(initial_temperatures_equal=True,final_T_rms_difference_K=difference(finals[name],finals[pg]),
    thermal_wall_ratio=results[name]['final']['heat_wall_max']/max(1e-20,results[pg]['final']['heat_wall_max']),
    rank_peak_memory_ratio=results[name]['final']['rss_max_KiB']/results[pg]['final']['rss_max_KiB'])
 refinements={}
 for scenario in ['cold','hot','assim']:
  for method in ['pg','pices']:
   key=f'{scenario}_{method}';names=[key+'_'+v for v in ['base','dt_half','dt_quarter']]
   if all(n in results for n in names):
    errors=[difference(finals[n],finals[names[2]]) for n in names[:2]]
    grids=[key+'_'+v for v in ['grid_coarse','base','grid_fine']]
    refinements[key]=dict(dt_difference_to_quarter_K=errors,dt_difference_decreased=errors[1]<=errors[0],
      grid_rotated_initial_L2=[results[n]['final']['rotated_initial_L2'] for n in grids] if scenario!='assim' and all(n in results for n in grids) else None,
      grid_final_diagnostics={n:{k:results[n]['final'][k] for k in ['energy_proxy','vrms','peak','wrong_sign_peak','anomaly_rms_radius']} for n in grids if n in results})
 sensitivities={}
 for scenario in ['cold','hot','assim']:
  base=f'{scenario}_pices_base'
  if base not in results:continue
  group={}
  for variant in ['np_low','np_high','length_half','length_double']:
   name=f'{scenario}_pices_{variant}'
   if name in results:
    require(initials[name]==initials[base],'sensitivity initial field mismatch '+name)
    b=results[base]['final'];v=results[name]['final']
    group[variant]=dict(final_T_rms_difference_from_base_K=difference(finals[name],finals[base]),
      rotated_initial_L2=v['rotated_initial_L2'],energy_proxy_difference_from_base=v['energy_proxy']-b['energy_proxy'],
      peak=v['peak'],wrong_sign_peak=v['wrong_sign_peak'],anomaly_rms_radius=v['anomaly_rms_radius'],
      heat_wall_ratio_to_base=v['heat_wall_max']/max(1e-20,b['heat_wall_max']))
  sensitivities[scenario]=group
 if not local:
  require((root/'launcher_exit_code.txt').read_text().strip()=='0','launcher failed')
  sc=(root/'solver_commit.txt').read_text().strip();require(re.fullmatch('[0-9a-f]{40}',sc)!=None,'solver commit')
  require((root/'runtime_source_diff.txt').read_text()=='','dirty runtime')
  require((root/'pices-p0-build.commit').read_text().strip()==sc,'stale build')
  require('P0_BUILD_COMPLETE commit='+sc in (root/'pices-p0-build.log').read_text(),'incomplete build')
  lines=(root/'input.sha256').read_text().splitlines();paths={str(p.relative_to(root)) for p in (root/'input').rglob('*') if p.is_file()}
  require({l.split()[1] for l in lines}==paths and len(lines)==len(paths),'input inventory')
  for l in lines:
   h,f=l.split();require(hashlib.sha256((root/f).read_bytes()).hexdigest()==h,'input hash '+f)
  for f in ['runs_commit.txt','binary.sha256','platform.txt','mpi_version.txt','submitted.lsf']:require((root/f).stat().st_size>0,'missing provenance '+f)
  require(re.fullmatch('[0-9a-f]{40}',(root/'runs_commit.txt').read_text().strip()) is not None,'runs commit')
  require(re.fullmatch('[0-9a-f]{64}',(root/'binary.sha256').read_text().split()[0]) is not None,'binary hash format')
 return dict(status='PARTIAL_PASS' if partial else 'PASS',production_decision='REQUIRES_REVIEW',cases_completed=len(complete),cases_expected=len(rows),provenance_checked=not local,cases=results,pairs=pairs,refinements=refinements,sensitivities=sensitivities,
  interpretation=['energy_proxy_change includes physical heating, boundary exchange, TA and remapping; it is not by itself conservation drift.','rotated_initial_L2 compares with a nodally sampled rotated initial field; diffusion is active, so this is a shape diagnostic, not an exact diffusion solution.','P5 completion does not authorize production switching; inspect oscillation, diffusion, refinement, dynamics, timing and memory together.'])
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('--local',action='store_true');p.add_argument('--partial',action='store_true');p.add_argument('--summary',type=Path);a=p.parse_args()
 try:r=verify(a.root,a.local,a.partial)
 except (ValueError,OSError,KeyError,AssertionError) as e:r=dict(status='FAIL',error=str(e))
 if a.summary:a.summary.write_text(json.dumps(r,indent=2)+'\n')
 print(json.dumps({k:v for k,v in r.items() if k not in ('cases','pairs','refinements','sensitivities')},indent=2));raise SystemExit(r['status']=='FAIL')
