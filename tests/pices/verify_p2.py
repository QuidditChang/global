"""P2 same-partition checkpoint integrity and continuous/restart audit (stdlib)."""
import argparse,gzip,hashlib,json,math,re,struct
from pathlib import Path

def require(ok,why):
    if not ok:raise ValueError(why)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def config(p):
    out={}
    for line in p.read_text().splitlines():
        line=line.split('#',1)[0].strip()
        if line:
            k,v=map(str.strip,line.split('=',1));require(k not in out,'duplicate cfg key');out[k]=v
    return out
def fields(s):return dict(re.findall(r'(\w+)=([^\s]+)',s))
def checkpoint(path):
    meta=json.loads(Path(str(path)+'.pices.json').read_text());b=path.read_bytes();state=Path(str(path)+'.pices.state')
    require(meta['magic']=='CITCOMS_EBA_PICES' and meta['schema']==1 and meta['phase']=='accepted','checkpoint schema')
    require(meta['checkpoint_sha256']==sha(path) and meta['accepted_velocity_sha256']==sha(state),'checkpoint checksum')
    nx,ny,nz=meta['local_mesh'];nn=nx*ny*nz;ne=(nx-1)*(ny-1)*(nz-1);neq=nn*3
    h=struct.unpack_from('=8i3f',b);require(list(h[:3])==[nx,ny,nz] and list(h[3:6])==meta['decomposition'] and h[6]==1,'binary mesh')
    require(h[7]==meta['step'] and h[8]==meta['time'] and h[9]==meta['dt'],'binary clock')
    require(meta['total_timesteps']==meta['step']+1,'restart counter')
    require(meta['extraq']==([{'name':'flavor','slot':0}] if meta['flavors'] else [])+[{'name':'Tp','slot':int(meta['flavors']>0)}],'attribute registry')
    off=44;live=[]
    def sentinel():
        nonlocal off
        require(b[off:off+16]==bytes(16),'binary sentinel');off+=16
    def array(n,size=8,skip=0):
        nonlocal off
        v=b[off:off+n*size];require(len(v)==n*size,'truncated array');off+=n*size;live.append(v[skip*size:]);return v
    sentinel();array(nn+1,skip=1);array(nn+1,skip=1)
    sentinel();array(2,4);array(ne+1,skip=1);array(neq)
    sentinel();header=struct.unpack_from('=5i',b,off);off+=20
    require(header[0]==12 and header[1]==len(meta['extraq']) and header[2]==meta['flavors'] and header[4]==meta['particles'],'binary tracer header')
    for i in range(6+header[1]):
        v=array(header[4]+1,skip=1);require(all(math.isfinite(x[0]) for x in struct.iter_unpack('=d',v[8:])),'nonfinite checkpoint particle')
    array(header[4]+1,4,1);require(off==len(b),'unexpected checkpoint size')
    require(state.stat().st_size==nn*3*4,'velocity companion size');live.append(state.read_bytes())
    return meta,live

def verify(root,local=False,stage="P2"):
    require(stage in ("P2","P3"),"invalid audit stage")
    root=Path(root);summary={};states={}
    for name,first,last in [('continuous',1,4),('split',1,2),('restart',3,4)]:
        d=root/name;cfg=config(d/'case.cfg');dt=struct.unpack('f',struct.pack('f',float(cfg['fixed_timestep'])))[0]
        for k,v in dict(energy_solver='pices',pices_checkpoint='on',tracer='1',tracer_flavors='0',nodex='5',nodey='5',nodez='5',nproc_surf='12',nprocx='1',nprocy='1',nprocz='1',Q0=('1' if stage=='P3' else '0'),dissipation_number=('0.1' if stage=='P3' else '0'),CBF_frequency='0',tracers_per_element='128').items():require(cfg.get(k)==v,'unexpected '+name+' cfg '+k)
        if stage=='P3':
            for k,v in dict(pices_eba='on',kT_exponent='0.3',kd_lower_prefactor='4.8',phase_depth='0.15,0.25,0.35',phase_delta_rho='10,10,10',phase_delta_s='-0.02,-0.02,0.02',phase_clapeyron='0.02,0.02,-0.02',phase_width='0.03,0.03,0.03').items():require(cfg.get(k)==v,'unexpected P3 cfg '+k)
        require(cfg['restart']==('on' if name=='restart' else 'off'),'restart switch')
        require(int((d/'mpi_exit_code.txt').read_text()) in (0,8),'MPI failed')
        require('PICES_ERROR' not in (d/'solver.stderr').read_text(),'PICES runtime error')
        counts={};cfl=[];initial=first-1
        for rank in range(12):
            folder=d/'DATA'/str(rank);log=(folder/'log').read_text()
            require(log.count('PICES_INIT ')==(0 if name=='restart' else 1),'fresh Tp initialization during restart')
            require(log.count('PICES_RESTORE ')==(1 if name=='restart' else 0),'restore record')
            rows=[fields(x) for x in log.splitlines() if x.startswith('PICES_STEP ')];require([int(x['step']) for x in rows]==list(range(first,last+1)),'thermal steps')
            for row in rows:
                require(all(math.isfinite(float(row[k])) for k in ['time','dt','cfl','Tmin','Tmax','mismatch','subgrid_max','remap_energy','heat_energy']),'nonfinite diagnostic')
                require(math.isclose(float(row['time']),int(row['step'])*dt,rel_tol=2e-7),'time mismatch')
                cfl.append(float(row['cfl']));require(0<=cfl[-1]<=.25,'CFL')
            if stage=='P3':
                heat=[fields(x) for x in log.splitlines() if x.startswith('PICES_EBA ')]
                require([int(x['step']) for x in heat]==list(range(first,last+1)),'missing EBA heat ledger')
                for row in heat:
                    vals=[float(row[k]) for k in ['storage','internal','adiabatic','viscous','phase_pressure','boundary_reaction']]
                    require(all(map(math.isfinite,vals)),'nonfinite EBA heat ledger')
                    require(vals[1]>0 and vals[3]>0 and vals[2]!=0 and vals[4]!=0,'inactive EBA source')
                    require(abs(vals[0]-sum(vals[1:]))<1e-12*max(1,sum(map(abs,vals))),'EBA ledger imbalance')
            coverage=[fields(x) for x in log.splitlines() if x.startswith('PICES_COVERAGE ')]
            require(len(coverage)==last-first+1 and all(x['empty_elements']=='0' and x['zero_boundary_nodes']=='0' for x in coverage),'coverage')
            for step in range(initial,last+1):
                with gzip.open(folder/str(step)/f'velo.{rank}.{step}.gz','rt') as f:lines=f.read().splitlines()
                require(len(lines)==127,'node count')
                for n,line in enumerate(lines[2:],1):
                    values=list(map(float,line.split()));require(len(values)==4 and all(map(math.isfinite,values)) and values[3]>=0,'grid values')
                    if n%5==1:require(values[3]==3700,'bottom BC')
                    if n%5==0:require(values[3]==300,'top BC')
                with gzip.open(folder/str(step)/f'tracer.{rank}.{step}.gz','rt') as f:lines=f.read().splitlines()
                h=lines[0].split();require(int(h[0])==step and int(h[2])==4 and int(h[1])==len(lines)-1,'tracer header');counts[step]=counts.get(step,0)+int(h[1])
                for line in lines[1:]:
                    values=list(map(float,line.split()));require(len(values)==4 and all(map(math.isfinite,values)) and 300+3400*values[3]>=0,'particle values')
        require(all(n==98304 for n in counts.values()),'particle count')
        # The legacy CG accepts either increment or the configured divergence criterion.
        blocks=[]
        for line in (d/'DATA/0/log').read_text().splitlines():
            if not line.startswith('AhatP '):continue
            m=re.search(r'AhatP \((\d+)\).*?div/v=(\S+) dv/v=(\S+) and dp/p=(\S+) for step (\d+)',line);require(m is not None,'bad convergence record')
            n,div,dv,dp,step=m.groups();row=dict(iterations=int(n),div_v=float(div),dv_v=float(dv),dp_p=float(dp),step=int(step))
            if int(n)==0:blocks.append(row)
            else:require(bool(blocks),'missing Stokes start');blocks[-1]=row
        require(set(x['step'] for x in blocks)>=set(range(first,last+1)),'missing Stokes solve')
        for row in blocks:
            require(all(math.isfinite(row[k]) for k in ['div_v','dv_v','dp_p']),'nonfinite Stokes')
            require(row['dv_v']<float(cfg['accuracy']) or row['dp_p']<float(cfg['accuracy']) or row['div_v']<float(cfg.get('tole_compressibility','0')),'Stokes not converged')
        for step in range(initial,last+1,2):
            hashes=[];metas=[]
            for rank in range(12):
                path=d/'DATA'/str(rank)/f'PICES_{stage}.chkpt.{rank}.{step}';meta,live=checkpoint(path);states[(name,step,rank)]=live;metas.append(meta)
                require(meta['rank']==rank and meta['mpi_size']==12 and meta['step']==step,'checkpoint rank/step');hashes.append(sha(Path(str(path)+'.pices.json')).encode()+b'\0')
            collective=hashlib.sha256(b''.join(hashes)).hexdigest()
            for rank in range(12):
                path=d/'DATA'/str(rank)/f'PICES_{stage}.chkpt.{rank}.{step}.pices.manifest';m=json.loads(path.read_text());require(m==dict(magic='CITCOMS_EBA_PICES_COMPLETE',schema=1,mpi_size=12,metadata_set_sha256=collective),'incomplete/mixed checkpoint set')
                if not local:require(metas[rank]['solver_commit']==(root/'solver_commit.txt').read_text().strip(),'checkpoint build commit')
        summary[name]=dict(particles=counts,max_cfl=max(cfl),stokes=blocks)
    compared=0
    for name,steps in [('split',[0,1,2]),('restart',[2,3,4])]:
        for rank in range(12):
            for step in steps:
                for field in ['velo','tracer','visc']:
                    rel=Path('DATA')/str(rank)/str(step)/f'{field}.{rank}.{step}.gz'
                    require(gzip.decompress((root/'continuous'/rel).read_bytes())==gzip.decompress((root/name/rel).read_bytes()),'field mismatch '+str(rel));compared+=1
            for step in ([0,2] if name=='split' else [2,4]):require(states[('continuous',step,rank)]==states[(name,step,rank)],'binary live state mismatch')
    if not local:
        require((root/'launcher_exit_code.txt').read_text().strip()=='0','launcher failed')
        require((root/'runtime_source_diff.txt').read_text()=='','uncommitted runtime source')
        commit=(root/'solver_commit.txt').read_text().strip();require(re.fullmatch('[a-f0-9]{40}',commit) is not None,'solver commit')
        require(commit==(root/'pices-p0-build.commit').read_text().strip(),'stale build')
        require('P0_BUILD_COMPLETE commit='+commit in (root/'pices-p0-build.log').read_text(),'build incomplete')
        hashes=(root/'input.sha256').read_text().splitlines();expected={f'cmbhf_EBA_PICES_{stage}.cfg',f'cmbhf_EBA_PICES_{stage}_split.cfg',f'cmbhf_EBA_PICES_{stage}_restart.cfg',f'refstate_EBA_PICES_{stage}.txt'}
        require({line.split()[1] for line in hashes}==expected and len(hashes)==4,'input inventory')
        for line in hashes:
            digest,file=line.split();require(digest==sha(root/file),'input hash')
        for name,cfg in [('continuous',f'cmbhf_EBA_PICES_{stage}.cfg'),('split',f'cmbhf_EBA_PICES_{stage}_split.cfg'),('restart',f'cmbhf_EBA_PICES_{stage}_restart.cfg')]:require((root/name/'case.cfg').read_bytes()==(root/cfg).read_bytes(),'actual cfg differs from archived input')
        for name in ['continuous','split','restart']:
            require((root/name/f'refstate_EBA_PICES_{stage}.txt').read_bytes()==(root/f'refstate_EBA_PICES_{stage}.txt').read_bytes(),'actual reference state differs from archived input')
        for f in ['runs_commit.txt','binary.sha256','platform.txt','mpi_version.txt','submitted.lsf']:require((root/f).stat().st_size>0,'missing provenance '+f)
    return dict(status='PASS',scope=stage+' same-partition restart',provenance_checked=not local,cases=summary,decoded_outputs_compared=compared,binary_live_state_equal=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path);p.add_argument('--stage',choices=['P2','P3'],default='P2');p.add_argument('--local',action='store_true');p.add_argument('--summary',type=Path);a=p.parse_args()
    try:r=verify(a.root,a.local,a.stage)
    except (ValueError,OSError,KeyError,IndexError,struct.error,EOFError) as e:r=dict(status='FAIL',error=str(e))
    text=json.dumps(r,indent=2)+'\n';print(text)
    if a.summary:a.summary.write_text(text)
    raise SystemExit(r['status']!='PASS')
