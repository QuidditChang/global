"""Audit P8a composition, conduction/buoyancy controls and accepted restart."""
import argparse,gzip,hashlib,json,math,re,struct
from collections import Counter
from pathlib import Path
from verify_p2 import checkpoint,require,fields,sha

def verify(root,local=False):
 rows=json.loads((root/'input/matrix.json').read_text())['cases'];require([(x['name'],x['first'],x['last']) for x in rows]==[('continuous',1,6),('split',1,2),('restart',3,6),('conductivity_control',1,6),('buoyancy_control',1,6)],'matrix definition')
 summary={};states={};flavors={};field_count=0
 for row in rows:
  name=row['name'];d=root/name;first,last=row['first'],row['last'];initial=first-1
  require((d/'case.cfg').read_bytes()==(root/'input'/(name+'.cfg')).read_bytes(),'cfg mismatch')
  require((d/'refstate_EBA_PICES_P4.txt').read_bytes()==(root/'input/refstate_EBA_PICES_P4.txt').read_bytes(),'reference mismatch')
  for f in (root/'input/pices_p4_forcing').iterdir():
   if f.is_file() and not f.name.startswith('._'):require((d/'pices_p4_forcing'/f.name).read_bytes()==f.read_bytes(),'forcing mismatch')
  require((d/'mpi_exit_code.txt').read_text().strip() in ('0','8'),'MPI exit')
  require('PICES_ERROR' not in (d/'solver.stderr').read_text(),'PICES error')
  hist={s:Counter() for s in range(initial,last+1)};maxcfl=0;ledger=0
  for rank in range(12):
   folder=d/'DATA'/str(rank);log=(folder/'log').read_text().splitlines()
   records=[fields(x) for x in log if x.startswith('PICES_STEP ')]
   require([int(x['step']) for x in records]==list(range(first,last+1)),'missing/duplicate steps')
   for x in records:
    require(all(math.isfinite(float(v)) for k,v in x.items() if k not in ('derivative',)),'nonfinite step')
    maxcfl=max(maxcfl,float(x['cfl']));require(float(x['cfl'])<=.25,'CFL')
   comp=[fields(x) for x in log if x.startswith('PICES_COMPOSITION ')]
   require([int(x['step']) for x in comp]==list(range(initial,last+1)),'composition update count')
   require(all(float(x['error'])<=1e-12 for x in comp),'composition mismatch')
   for x in [fields(x) for x in log if x.startswith('PICES_PROJECTION ')]:
    require(math.isfinite(float(x['residual'])) and 0<=float(x['residual'])<=float(x['tolerance']),'projection convergence')
   stokes=[fields(x) for x in log if x.startswith('PICES_STOKES ')]
   if rank==0:require({int(x['step']) for x in stokes}>=set(range(first,last+1)) and all(x['status']=='PASS' for x in stokes),'Stokes')
   heats=[fields(x) for x in log if x.startswith('PICES_EBA ')]
   require(len(heats)==last-first+1,'heat ledger count')
   for x in heats:
    vals=[float(x[k]) for k in ['storage','internal','adiabatic','viscous','phase_pressure','boundary_reaction']]
    require(all(map(math.isfinite,vals)),'nonfinite heat');err=abs(vals[0]-sum(vals[1:]));require(err<1e-12*max(1,sum(map(abs,vals))),'heat imbalance');ledger=max(ledger,err)
   for step in range(initial,last+1):
    folder_step=folder/str(step)
    lines=gzip.decompress((folder_step/f'tracer.{rank}.{step}.gz').read_bytes()).decode().splitlines();header=lines[0].split()
    require(int(header[0])==step and int(header[1])==len(lines)-1 and int(header[2])==5,'tracer header')
    for line in lines[1:]:
     v=list(map(float,line.split()));require(len(v)==5 and all(map(math.isfinite,v)),'particle data')
     flavor=int(v[3]);require(v[3]==flavor and 0<=flavor<25 and flavor not in (18,19) and 300+3400*v[4]>=0,'flavor/Tp');hist[step][flavor]+=1
    for f in folder_step.glob('*.gz'):
     if f.name.startswith('._'):continue
     for token in gzip.decompress(f.read_bytes()).decode().split():
      try:v=float(token)
      except ValueError:continue
      require(math.isfinite(v),'nonfinite output '+str(f))
   for step in range(initial,last+1,2):
    path=folder/f'PICES_P4.chkpt.{rank}.{step}';meta,live=checkpoint(path);require(meta['schema']==3,'composition schema');states[name,rank,step]=live
    colors=[int(x[0]) for x in struct.iter_unpack('=d',live[11])];elements=[x[0] for x in struct.iter_unpack('=i',live[13])]
    counts=Counter(elements);pairs=Counter(zip(elements,colors))
    require(len(counts)==64,'empty checkpoint element')
    for i in range(24):
     values=[x[0] for x in struct.iter_unpack('=d',live[16+i])]
     require(len(values)==64,'composition elements')
     require(all(abs(v-pairs[e,i+1]/counts[e])<1e-12 for e,v in enumerate(values,1)),'checkpoint ratio reconstruction')
    if not local:require(meta['solver_commit']==(root/'solver_commit.txt').read_text().strip(),'checkpoint commit')
  require(all(sum(h.values())==98304 for h in hist.values()),'global count')
  require(len({h[24] for h in hist.values()})==1 and hist[initial][24]>0,'primordial conservation')
  flavors[name]=hist
  summary[name]=dict(max_cfl=maxcfl,max_heat_error=ledger,initial_flavors=dict(hist[initial]),final_flavors=dict(hist[last]))
 # Decode gzip before comparison: container timestamps are not solver state.
 for name,end,start in [('split',2,0),('restart',6,2)]:
  for rank in range(12):
   for step in range(start,end+1):
    sub=Path('DATA')/str(rank)/str(step)
    files={f.name for f in (root/'continuous'/sub).glob('*.gz') if not f.name.startswith('._')}
    require(files=={f.name for f in (root/name/sub).glob('*.gz') if not f.name.startswith('._')},'output inventory')
    for f in files:
     require(gzip.decompress((root/'continuous'/sub/f).read_bytes())==gzip.decompress((root/name/sub/f).read_bytes()),'restart field '+str(sub/f));field_count+=1
   for step in range(start,end+1,2):require(states['continuous',rank,step]==states[name,rank,step],'restart live state')
 def temperature(name):
  vals=[]
  for rank in range(12):
   lines=gzip.decompress((root/name/f'DATA/{rank}/6/velo.{rank}.6.gz').read_bytes()).decode().splitlines()[2:]
   vals.extend(float(x.split()[3]) for x in lines)
  return vals
 ref=temperature('continuous');deltas={}
 for name in ('conductivity_control','buoyancy_control'):
  other=temperature(name);delta=math.sqrt(sum((x-y)**2 for x,y in zip(ref,other))/len(ref));require(delta>1e-8,'inactive control '+name);deltas[name]=delta
 require(flavors['continuous'][0]!=flavors['continuous'][6],'inactive flavor reclassification')
 from verify_p4 import verify_p4
 cbf=verify_p4(root,True,6,2)
 # Every checkpoint set must have a matching collective completion manifest.
 for row in rows:
  name=row['name']
  for step in range(row['first']-1,row['last']+1,2):
   paths=[root/name/f'DATA/{rank}/PICES_P4.chkpt.{rank}.{step}' for rank in range(12)]
   collective=hashlib.sha256(b''.join(sha(Path(str(p)+'.pices.json')).encode()+b'\0' for p in paths)).hexdigest()
   for p in paths:require(json.loads(Path(str(p)+'.pices.manifest').read_text())==dict(magic='CITCOMS_EBA_PICES_COMPLETE',schema=1,mpi_size=12,metadata_set_sha256=collective),'manifest')
 if not local:
  require((root/'launcher_exit_code.txt').read_text().strip()=='0','launcher exit')
  require((root/'p8a_complete.txt').read_text().strip()=='P8A_COMPLETE_PENDING_LOCAL_AUDIT','completion')
  require((root/'runtime_source_diff.txt').read_text()=='','runtime source diff')
  commit=(root/'solver_commit.txt').read_text().strip()
  require(re.fullmatch('[a-f0-9]{40}',commit) is not None and commit==(root/'pices-p0-build.commit').read_text().strip(),'build commit')
  require('P0_BUILD_COMPLETE commit='+commit in (root/'pices-p0-build.log').read_text(),'incomplete build')
  require((root/'completed_cases.txt').read_text().strip()=='5','completed case count')
  hashed={}
  for line in (root/'input.sha256').read_text().splitlines():
   digest,file=line.split();require(file not in hashed and file.startswith('input/') and '..' not in Path(file).parts,'input hash path');hashed[file]=digest
   require(sha(root/file)==digest,'input hash')
  require(set(hashed)=={str(f.relative_to(root)) for f in (root/'input').rglob('*') if f.is_file()},'input inventory')
  mpi=(root/'mpi_path.txt').read_text().strip();library=(root/'binary_ldd.txt').read_text()
  require('/mpi/intel64/bin/mpiexec.hydra' in mpi and str(Path(mpi).parent.parent/'lib/libmpi.so.12') in library,'MPI provider match')
  require('Intel(R) MPI' in (root/'mpi_version.txt').read_text(),'MPI version')
  for f in ('submitted.lsf','runs_commit.txt','binary.sha256','platform.txt'):require((root/f).stat().st_size>0,'missing provenance '+f)
 return dict(status='PASS',TA_CBF=cbf,cases=summary,restart_outputs_compared=field_count,binary_live_state_equal=True,control_temperature_RMS_K=deltas,production_decision='HOLD_PENDING_P8B_P8C')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('--local',action='store_true');p.add_argument('--summary',type=Path);a=p.parse_args()
 try:r=verify(a.root,a.local)
 except (ValueError,OSError,KeyError,IndexError,struct.error) as e:r=dict(status='FAIL',error=str(e))
 text=json.dumps(r,indent=2)+'\n'
 if a.summary:a.summary.write_text(text)
 print(text);raise SystemExit(r['status']!='PASS')
