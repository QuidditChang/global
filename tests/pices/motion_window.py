"""Prepare/audit a frozen-field replay window. No writes to the source run."""
import argparse, hashlib, json, math, re, shutil, struct
from pathlib import Path

SOLVER_COMMIT='92f6fa44b52fec25fe283fa76c8be08d90b99499'

def fields(line):
 return dict(re.findall(r'(\w+)=([^\s]+)',line))

def parameters(text):
 # The C reader consumes key=value records regardless of section headers.
 result={}
 for line in text.splitlines():
  m=re.match(r'^\s*([\w]+)\s*=\s*(.*?)\s*$',line)
  if m: result[m[1]]=m[2]
 return result

def replace(text,key,value):
 pattern=r'(?m)^\s*'+re.escape(key)+r'\s*=.*$'
 if not re.search(pattern,text):raise ValueError('missing runtime parameter '+key)
 return re.sub(pattern,lambda m:key+'='+str(value),text)

def standalone_runtime(text):
 """Map Pyre's unset conductivity sentinels to absent standalone parameters.

 module/setProperties.c treats exactly -1 as unspecified. The standalone
 reader instead rejects any explicitly supplied value inconsistent with ks/k0.
 Keep explicit values (including invalid ones) so its validation still runs.
 """
 removed=[]
 for key in ('inputdiffusivity','reference_conductivity'):
  value=parameters(text).get(key)
  if value is not None and float(value)==-1.0:
   text=re.sub(r'(?m)^\s*'+key+r'\s*=.*$', '',text)
   removed.append(key)
 return text,removed

def restore_input_precision(text,original):
 """pid*.cfg is a %g display, not a lossless serialization (getProperty.h).
 Restore explicit numeric input only when its float/double C display agrees.
 Retain Pyre's expanded defaults and normalized booleans. Reject unexplained
 numeric differences rather than silently replacing transformed physics.
 """
 restored=[]
 for key,value in parameters(original).items():
  shown=parameters(text).get(key)
  if shown is None:continue
  try:
   actual=[float(x.strip()) for x in value.split(',')]
   rounded=[float(x.strip()) for x in shown.split(',')]
  except ValueError:continue
  if actual==rounded:continue
  if len(actual)!=len(rounded):raise ValueError('runtime vector length differs: '+key)
  for a,b in zip(actual,rounded):
   if not math.isfinite(a) or not math.isfinite(b):raise ValueError('nonfinite runtime input: '+key)
   f=struct.unpack('=f',struct.pack('=f',a))[0]
   if b not in (a,float(format(a,'.6g')),float(format(f,'.6g'))):
    raise ValueError('runtime differs from original beyond %g rounding: '+key)
  text=replace(text,key,value);restored.append(key)
 # setup_parser stops unquoted values at whitespace, including after commas.
 # Normalize every numeric record, not only values restored above. Preserve
 # each original numeric token verbatim: formatting via float would lose bits.
 for key,value in parameters(text).items():
  tokens=[v.strip() for v in value.split(',')]
  try:[float(v) for v in tokens]
  except ValueError:continue
  text=replace(text,key,','.join(tokens))
 return text,restored

def standalone_mesh(text):
 """Mirror module/setProperties.c Mesher_set_properties, not C defaults=1."""
 p=parameters(text);levels=int(p['levels'])
 if levels<1:raise ValueError('invalid multigrid levels')
 units={};local=[]
 for axis in 'xyz':
  nodes=int(p['node'+axis]);ranks=int(p['nproc'+axis]);denom=ranks*2**(levels-1)
  if ranks<1 or nodes<2 or (nodes-1)%denom:raise ValueError('nonintegral multigrid mesh: '+axis)
  key='mgunit'+axis;units[key]=(nodes-1)//denom;local.append((nodes-1)//ranks+1)
  if key in p and int(p[key])!=units[key]:raise ValueError('inconsistent '+key)
  text=replace(text,key,units[key]) if key in p else text+'\n'+key+'='+str(units[key])+'\n'
 return text,units,local

def checkpoint_state_bytes(p,local):
 # Schema layout in Pices_checkpoint.c; fail before MPI if cfg cannot read it.
 def enabled(key,default='0'):return p.get(key,default).lower() in ('1','on','true','yes')
 n=math.prod(local);e=math.prod(v-1 for v in local);surface=int(p['nodex'])*int(p['nodey'])
 p4=enabled('pices_p4');coupled=int(p.get('rheol','3'))==7 or enabled('file_vbcs') or int(p.get('qvis_mode','0'))!=0 or float(p.get('fixed_timestep','0'))<=0
 composition=enabled('tracer') and (enabled('chemical_buoyancy','1') or enabled('CDEPV'))
 return (12*n+(4*(3*n+(surface if enabled('lith_age') else 0))+4+64*e if p4 else 0)
         +(4*(8*e+n) if coupled else 0)+(8+(4*surface if enabled('tracer_reclassify_flavors') else 0) if composition else 0))

def prepare(source,root,template,build,first=91,last=100,ranks=384):
 if root.exists():raise ValueError('refusing to reuse '+str(root))
 source=source.resolve();root=root.resolve();build=build.resolve()
 provenance=json.loads((build/'replay_build.json').read_text())
 entry=provenance.get('entry','standalone')
 if entry not in ('standalone','producer_pyre_objects_v1'):raise ValueError('unknown replay entry')
 if provenance['solver_commit']!=SOLVER_COMMIT or not provenance['isolated_compile']:
  raise ValueError('rebuild isolated replay at the recorded solver commit')
 if (source/'runtime_source_diff.txt').exists() and (source/'runtime_source_diff.txt').stat().st_size:
  raise ValueError('source run had uncommitted runtime changes')
 if (source/'solver_commit.txt').read_text().strip()!=SOLVER_COMMIT:
  raise ValueError('unexpected source solver commit')
 exe=build/'PICESMotionReplay'
 if hashlib.sha256(exe.read_bytes()).hexdigest()!=provenance['binary_sha256']:
  raise ValueError('replay binary/build manifest mismatch')
 if entry=='producer_pyre_objects_v1':
  files=provenance.get('pyre_files',{})
  actual={str(p.relative_to(build/'PyreSite')) for p in (build/'PyreSite').rglob('*') if p.is_file()}
  if not files or actual!=files.keys():raise ValueError('Pyre package snapshot inventory mismatch')
  for name,digest in files.items():
   if hashlib.sha256((build/'PyreSite'/name).read_bytes()).hexdigest()!=digest:raise ValueError('Pyre package snapshot checksum mismatch')
 # Read the serialized Pyre inventory; adapt its unset sentinels below.
 runtime=list(source.glob('pid[0-9]*.cfg'))
 if len(runtime)!=1:raise ValueError('expected exactly one source pid*.cfg; found '+str(runtime))
 text=runtime[0].read_text();p=parameters(text);t=parameters(template.read_text());original=parameters((source/'input.cfg').read_text())
 allowed={'steps','walltime','name','datadir','datadir_old','restart','solution_cycles_init','refstate_file','coor_file'}
 differences={k:(t.get(k),original.get(k)) for k in t.keys()|original.keys() if t.get(k)!=original.get(k) and k not in allowed}
 if differences:raise ValueError('template changed physics: '+str(differences))
 text,restored_parameters=restore_input_precision(text,(source/'input.cfg').read_text())
 text,unset_parameters=standalone_runtime(text)
 text,mg_units,expected_mesh=standalone_mesh(text)
 state_bytes=checkpoint_state_bytes(parameters(text),expected_mesh)
 if int(p['nproc_surf'])*int(p['nprocx'])*int(p['nprocy'])*int(p['nprocz'])!=ranks:raise ValueError('rank layout mismatch')
 clocks={};local_mesh=None
 for rank in range(ranks):
  directory=source/'DATA'/str(rank);stem=directory/('global.chkpt.%d.%d'%(rank,first-1))
  for suffix in ('','.pices.json','.pices.state','.pices.manifest'):
   if Path(str(stem)+suffix).stat().st_size<=0:raise ValueError('empty checkpoint '+str(stem))
  meta=json.loads(Path(str(stem)+'.pices.json').read_text())
  if meta['solver_commit']!=SOLVER_COMMIT or meta['step']!=first-1:raise ValueError('checkpoint identity mismatch')
  if not {'motion_id','motion_selected'}<={v['name'] for v in meta['extraq']}:raise ValueError('missing motion IDs')
  if meta['local_mesh']!=expected_mesh:raise ValueError('checkpoint/config local mesh mismatch: '+str(stem))
  with stem.open('rb') as f:header=struct.unpack('=8i',f.read(32))
  if header!=tuple(expected_mesh+[int(p['nproc'+a]) for a in 'xyz']+[1,first-1]):raise ValueError('checkpoint header mesh/step mismatch: '+str(stem))
  if Path(str(stem)+'.pices.state').stat().st_size!=state_bytes:raise ValueError('checkpoint state size/config mismatch: '+str(stem))
  n=math.prod(meta['local_mesh']);local_mesh=meta['local_mesh'];expected=meta['time']
  for step in range(first,last+1):
   path=directory/('global.motion_velocity.%d.%d'%(rank,step))
   with path.open('rb') as f: head=struct.unpack('=4i2d',f.read(32))
   if head[:4]!=(0x504d5631,rank,step,n) or path.stat().st_size!=32+12*n:raise ValueError('invalid velocity dump '+str(path))
   clock,dt=head[4:]
   if clock!=expected or not math.isfinite(dt) or dt<=0:raise ValueError('discontinuous replay clock '+str(path))
   if rank==0:clocks[step]=[clock,dt]
   elif clocks[step]!=[clock,dt]:raise ValueError('rank clock mismatch')
   expected=struct.unpack('=f',struct.pack('=f',clock+dt))[0]
 print('checkpoint and velocity inventory verified',flush=True)
 refs={};dest={};ids=set();low_cells={s:[] for s in range(first,last+1)}
 for rank in range(ranks):
  with (source/'DATA'/str(rank)/('global.motion.%d'%rank)).open() as f:
   for line in f:
    if line.startswith('CELL '):
     _,step,element,count=line.split();step=int(step)-1
     if first<=step<=last:low_cells[step].append([rank,int(element),int(count)])
     continue
    if not line.startswith(('PATH ','HOST ')):continue
    v=line.split();step=int(v[1])
    if not first<=step<=last or v[2]!='2':continue
    key=(step,int(v[3]))
    if v[0]=='PATH':
     if key in refs:raise ValueError('duplicate source path')
     refs[key]=tuple(map(float,v[12:15]));ids.add(key[1])
    else:
     if key in dest or v[5]!='1':raise ValueError('duplicate/invalid source host')
     dest[key]=(rank,int(v[4]))
 if refs.keys()!=dest.keys() or not refs:raise ValueError('source trajectory pairing failed')
 source_rows=(source/'DATA/0/global.log').read_text().splitlines()
 hist={int(f['step']):list(map(int,f['histogram'].split(','))) for line in source_rows if line.startswith('PICES_COVERAGE_STATS ') for f in [fields(line)]}
 signatures={tuple(f[k] for k in ('count','hash_sum','hash_xor')) for line in source_rows if line.startswith('PICES_MOTION ') for f in [fields(line)]}
 if len(signatures)!=1 or any(s not in hist for s in range(first,last+1)):raise ValueError('incomplete source census')
 # Native Pyre reads the original sectioned input, not the lossy pid dump or
 # the standalone bridge. Only restart/output/short local input paths change.
 if entry=='producer_pyre_objects_v1':text=(source/'input.cfg').read_text()
 root.mkdir(parents=True)
 if entry=='producer_pyre_objects_v1':shutil.copytree(build/'PyreSite',root/'PyreSite')
 shutil.copy2(exe,root/'PICESMotionReplay');shutil.copy2(build/'replay_build.json',root)
 shutil.copy2(template,root/'submitted.cfg');shutil.copy2(source/'input.cfg',root/'source_input.cfg');shutil.copy2(runtime[0],root/'source_runtime.cfg')
 for name,key in [('refstate_EBA.txt','refstate_file'),('GLB.coor.global.dat','coor_file')]:
  path=Path(p[key]);path=path if path.is_absolute() else source/path
  shutil.copy2(path,root/name);text=replace(text,key,str(root/name))
 # Use a short relative mesh filename (the solver's coor_file buffer is 100 bytes).
 text=replace(text,'coor_file','GLB.coor.global.dat');text=replace(text,'refstate_file','refstate_EBA.txt')
 for name in ('Inca_plateau_15Ma_v5.xy','Nazca_ridge_15Ma_v2.xy','Juan_Fernandez_Ridge_15Ma_v2.xy'):
  shutil.copy2(source/name,root/name)
 (root/'tracked_ids.txt').write_text(''.join(str(uid)+'\n' for uid in sorted(ids)))
 with (root/'reference.tsv').open('w') as f:
  for (step,uid),xyz in sorted(refs.items()):f.write('%d %d %.17g %.17g %.17g %d %d\n'%((step,uid)+xyz+dest[step,uid]))
 for factor in (1,2,4):
  case=root/('sub%d'%factor);case.mkdir();(case/'DATA').mkdir()
  config=text
  for k,v in dict(datadir=str(case/'DATA/%RANK'),datadir_old=str(source/'DATA/%RANK'),restart='on' if entry=='producer_pyre_objects_v1' else 1,solution_cycles_init=first-1).items():config=replace(config,k,v)
  (case/'case.cfg').write_text(config)
 manifest=dict(source=str(source),entry=entry,first=first,last=last,ranks=ranks,solver_commit=SOLVER_COMMIT,tracked=len(ids),reference_rows=len(refs),standalone_unset_parameters=unset_parameters,restored_precision_parameters=restored_parameters,standalone_mg_units=mg_units,checkpoint_state_bytes=state_bytes,histograms=hist,low_cells=low_cells,signature=list(next(iter(signatures))),clocks=clocks,local_mesh=local_mesh,
  scope='Prescribed time sequence of saved velocity fields; each field frozen within its original timestep. No coupled physics evolution.')
 (root/'replay_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 print('prepared',root,'tracked IDs',len(ids),flush=True)

def read_positions(root,step,ranks):
 result={}
 for rank in range(ranks):
  path=root/'DATA'/str(rank)/('global.replay.%d.%d'%(rank,step))
  for line in path.open():
   if line.startswith('#'):continue
   uid,x,y,z,host=line.split();uid=int(uid);xyz=tuple(map(float,(x,y,z)))
   if uid in result or not all(math.isfinite(v) for v in xyz):raise ValueError('duplicate ID/nonfinite endpoint')
   result[uid]=(xyz,rank,int(host))
 return result

def audit(root,baseline_only=False):
 import numpy as np
 factors=(1,) if baseline_only else (1,2,4)
 m=json.loads((root/'replay_manifest.json').read_text());steps=range(m['first'],m['last']+1);ranks=m['ranks'];ids={int(s) for s in (root/'tracked_ids.txt').read_text().split()}
 references={s:{} for s in steps}
 for line in (root/'reference.tsv').open():
  s,uid,x,y,z,rank,host=line.split();references[int(s)][int(uid)]=(tuple(map(float,(x,y,z))),int(rank),int(host))
 report={'status':'PASS_REVIEW_REQUIRED','issues':[],'baseline_max_position_error_nd':0.,'baseline_host_mismatches':0,'steps':[],'scope':m['scope']}
 logs={}
 for factor in factors:
  case=root/('sub%d'%factor);lines=(case/'DATA/0/global.log').read_text().splitlines()
  if (case/'returncode.txt').read_text().strip()!='0':raise ValueError('replay failed: '+str(factor))
  rows=[fields(l) for l in lines if l.startswith('PICES_MOTION ')]
  expected=[(s,k) for s in steps for _ in range(factor) for k in (1,2)]
  if [(int(v['step']),int(v['stage'])) for v in rows]!=expected:raise ValueError('incomplete movement stages')
  if any([v[k] for k in ('count','hash_sum','hash_xor')]!=m['signature'] for v in rows):raise ValueError('particle signature mismatch')
  census=[fields(l) for l in lines if l.startswith('PICES_REPLAY ')]
  if [int(v['step']) for v in census]!=list(steps):raise ValueError('incomplete replay census')
  logs[factor]={int(v['step']):v for v in census}
 for step in steps:
  states=[read_positions(root/('sub%d'%factor),step,ranks) for factor in factors]
  if any(v.keys()!=ids for v in states):raise ValueError('tracked ID set changed')
  for uid,expected in references[step].items():
   actual=states[0][uid];err=math.dist(actual[0],expected[0]);report['baseline_max_position_error_nd']=max(err,report['baseline_max_position_error_nd']);report['baseline_host_mismatches']+=actual[1:]!=expected[1:]
  hist=[list(map(int,logs[f][step]['histogram'].split(','))) for f in factors]
  if hist[0]!=m['histograms'][str(step)]:report['issues'].append('baseline histogram differs at step '+str(step))
  if any(int(logs[f][step]['tracked'])!=len(ids) for f in factors):raise ValueError('tracked census mismatch')
  row={'step':step,'empty_counts':dict(zip(map(str,factors),[v[0] for v in hist])),'histograms':dict(zip(map(str,factors),hist)),'endpoint_error_to_4_m':{},'host_changes_to_4':{},'focus_cell_count':{}}
  for index,factor in enumerate(() if baseline_only else (1,2)):
   delta=np.array([math.dist(states[index][uid][0],states[2][uid][0])*6371000 for uid in ids]);row['endpoint_error_to_4_m'][str(factor)]=dict(rms=float(np.sqrt(np.mean(delta**2))),p99=float(np.quantile(delta,.99)),maximum=float(delta.max()))
   row['host_changes_to_4'][str(factor)]=sum(states[index][uid][1:]!=states[2][uid][1:] for uid in ids)
  empty_sets={}
  for factor in factors:
   low={}
   for rank in range(ranks):
    path=root/('sub%d'%factor)/'DATA'/str(rank)/('global.replay_cells.%d.%d'%(rank,step))
    for line in path.open():
     element,count=map(int,line.split())
     if (rank,element) in low or not 0<=count<=8:raise ValueError('invalid sparse census')
     low[rank,element]=count
   empty_sets[factor]={k for k,v in low.items() if v==0}
   if [sum(v==i for v in low.values()) for i in range(9)]!=hist[factors.index(factor)][:9]:raise ValueError('sparse census/global histogram disagreement')
   if factor==1 and {k:v for k,v in low.items() if v<=4}!={(r,e):n for r,e,n in m['low_cells'][str(step)]}:report['issues'].append('baseline low-cell map differs at step '+str(step))
   # Count files include every cell with <=8 particles. Absence means >=9.
   path=root/('sub%d'%factor)/'DATA/161'/('global.replay_cells.161.%d'%step)
   row['focus_cell_count'][str(factor)]=next((int(v[1]) for l in path.open() if (v:=l.split())[0]=='27151'),'>=9') if path.exists() else None
  if not baseline_only:row['empty_cell_symmetric_difference_to_4']={str(f):len(empty_sets[f]^empty_sets[4]) for f in (1,2)}
  report['steps'].append(row)
 if report['baseline_max_position_error_nd']>1e-12 or report['baseline_host_mismatches']:report['issues'].append('one-substep source reproduction failed')
 if report['issues']:report['status']='FAIL_BASELINE'
 report['interpretation']='Compare 1/2/4 convergence only after baseline passes. No automatic particle-addition or production-release decision. MPI decomposition is unchanged.'
 report['baseline_only']=baseline_only
 report_name='baseline_audit.json' if baseline_only else 'replay_audit.json'
 (root/report_name).write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
 if report['issues']:raise ValueError('baseline replay differs from original; see '+report_name)

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='mode',required=True)
 q=sub.add_parser('prepare');q.add_argument('--source',type=Path,required=True);q.add_argument('--root',type=Path,required=True);q.add_argument('--template',type=Path,required=True);q.add_argument('--build',type=Path,required=True)
 q=sub.add_parser('audit');q.add_argument('root',type=Path);q.add_argument('--baseline-only',action='store_true')
 a=p.parse_args()
 if a.mode=='prepare':prepare(a.source,a.root,a.template,a.build)
 else:audit(a.root,a.baseline_only)
