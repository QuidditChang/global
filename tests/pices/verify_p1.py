"""Audit the P1 HPC smoke run (stdlib, Python >= 3.8). No restart/P2 claim."""
import argparse,gzip,hashlib,json,math,re,struct
from pathlib import Path

def require(ok,why):
    if not ok: raise ValueError(why)
def fields(line): return dict(re.findall(r'(\w+)=([^\s]+)',line))
def verify(root,local=False):
    root=Path(root)
    require(int((root/'mpi_exit_code.txt').read_text()) in (0,8),'MPI failed')
    cfg={}
    for line in (root/'cmbhf_EBA_PICES_P1.cfg').read_text().splitlines():
        line=line.split('#',1)[0].strip()
        if line:
            k,v=map(str.strip,line.split('=',1));require(k not in cfg,'duplicate cfg key');cfg[k]=v
    for k,v in dict(energy_solver='pices',tracer='1',tracer_flavors='0',chemical_buoyancy='off',tracer_reclassify_flavors='off',nodex='5',nodey='5',nodez='5',Q0='0',dissipation_number='0',CBF_frequency='0',output_q_surf_CBF='off',output_q_botm_CBF='off',lith_age='0',maxstep='2',maxtotstep='3',output_format='ascii-gz',tracers_per_element='128').items():
        require(cfg.get(k)==v,'unexpected cfg '+k)
    require(cfg.get('pices_test_no_diffusion','off')=='off','smoke requires diffusion')
    dt=struct.unpack('f',struct.pack('f',float(cfg['fixed_timestep'])))[0]
    require(dt>0,'nonpositive dt')
    ranks={p.name for p in (root/'DATA').iterdir() if p.is_dir()};require(ranks=={str(i) for i in range(12)},'rank directories missing')
    totals=[0,0,0];temps=[];mismatch=[];subgrid=[];cfl=[]
    for rank in range(12):
        directory=root/'DATA'/str(rank);log=(directory/'log').read_text()
        init=[fields(l) for l in log.splitlines() if l.startswith('PICES_INIT ')];require(len(init)==1,'missing/duplicate PICES init')
        require(init[0]['Tp_slot']=='0' and float(init[0]['kappa'])>0,'wrong Tp registration/diffusivity')
        rows=[fields(l) for l in log.splitlines() if l.startswith('PICES_STEP ')];require(len(rows)==2,'wrong PICES step count')
        audits=[fields(l) for l in log.splitlines() if l.startswith('TEMP_AUDIT stage=thermal_exit ')];require(len(audits)==2,'wrong thermal exit count')
        for step,rec in enumerate(rows,1):
            require(int(rec['step'])==step and int(rec['rank'])==rank,'wrong step/rank')
            for k in ('time','dt','cfl','Tmin','Tmax','mismatch','subgrid_max','remap_energy','heat_energy'):
                require(math.isfinite(float(rec[k])),'nonfinite '+k)
            require(math.isclose(float(rec['time']),step*dt,rel_tol=1e-6,abs_tol=1e-14),'wrong PICES clock')
            require(math.isclose(float(rec['dt']),dt,rel_tol=1e-6),'wrong dt')
            require(int(rec['substeps'])>=1 and 0<=float(rec['cfl'])<=.25,'unstable step')
            require(float(rec['Tmin'])>=-300/3400,'negative Kelvin')
            require(audits[step-1]['negative']=='0' and audits[step-1]['nonfinite']=='0','invalid thermal audit')
            mismatch.append(float(rec['mismatch']));subgrid.append(float(rec['subgrid_max']));cfl.append(float(rec['cfl']))
        require('PICES_ERROR' not in log and 'TEMP_AUDIT_NODE' not in log,'error in log')
        require(not list(directory.glob('*.chkpt.*')),'P1 must not publish restart checkpoint')
        require(not list(directory.glob('q.surf.*')) and not list(directory.glob('q.botm.*')),'legacy CBF output unexpected')
        for step in range(3):
            with gzip.open(directory/str(step)/('velo.%d.%d.gz'%(rank,step)),'rt') as f: lines=f.read().splitlines()
            h=lines[0].split();require(len(lines)==127 and int(h[0])==step and int(h[1])==125,'invalid node output')
            require(math.isclose(float(h[2]),step*dt,rel_tol=1e-5,abs_tol=1e-14),'wrong node clock')
            require(list(map(float,h[3:]))==[300,3700,3400],'wrong normalization')
            for i,line in enumerate(lines[2:],1):
                v=list(map(float,line.split()));require(len(v)==4 and all(map(math.isfinite,v)),'nonfinite node')
                require(v[3]>=0,'negative Kelvin');temps.append(v[3])
                if i%5==1:require(abs(v[3]-3700)<1e-5,'bottom boundary')
                if i%5==0:require(abs(v[3]-300)<1e-5,'top boundary')
            with gzip.open(directory/str(step)/('tracer.%d.%d.gz'%(rank,step)),'rt') as f:lines=f.read().splitlines()
            h=lines[0].split();count=int(h[1]);totals[step]+=count
            require(int(h[0])==step and int(h[2])==4 and len(lines)==count+1,'invalid Tp output')
            require(math.isclose(float(h[3]),step*dt,rel_tol=1e-5,abs_tol=1e-14),'wrong particle clock')
            if step: require(count==int(rows[step-1]['particles']),'particle count/log mismatch')
            else:require(count==int(init[0]['ntracers']),'particle init mismatch')
            for line in lines[1:]:
                v=list(map(float,line.split()));require(len(v)==4 and all(map(math.isfinite,v)),'invalid particle')
                require(300+3400*v[3]>=0,'negative particle Kelvin')
    require(totals==[98304]*3,'particle number not conserved')
    require(max(subgrid)>0,'subgrid correction never exercised')
    provenance={}
    if not local:
        for name in ('solver_commit.txt','runs_commit.txt','pices-p0-build.commit'):
            value=(root/name).read_text().strip();require(re.fullmatch('[0-9a-f]{40}',value)!=None,'bad commit '+name);provenance[name]=value
        require(provenance['solver_commit.txt']==provenance['pices-p0-build.commit'],'stale build')
        require((root/'launcher_exit_code.txt').read_text().strip()=='0','launcher failed')
        require((root/'runtime_source_diff.txt').read_text()=='','runtime source modified')
        for name in ('mpi_version.txt','platform.txt','submitted.lsf','pices-p0-build.log'):
            require((root/name).stat().st_size>0,'empty/missing '+name)
        require('P0_BUILD_COMPLETE commit='+provenance['solver_commit.txt'] in (root/'pices-p0-build.log').read_text(),'missing build success marker')
        require(re.fullmatch(r'[0-9a-f]{64}\s+\S[^\n]*\n?',(root/'binary.sha256').read_text())!=None,'invalid binary hash')
        require((root/'solver_status.txt').is_file() and (root/'runs_status.txt').is_file(),'missing git status')
        hashes=(root/'input.sha256').read_text().splitlines();require(len(hashes)==2,'wrong input hashes')
        for line,name in zip(hashes,('cmbhf_EBA_PICES_P1.cfg','refstate_EBA_PICES_P1.txt')):
            digest,record=line.split(maxsplit=1);require(record==name and digest==hashlib.sha256((root/name).read_bytes()).hexdigest(),'input hash mismatch')
    return dict(status='PASS',scope='P1 constant-physics smoke; not P2 restart/seam validation',provenance_checked=not local,particles=totals,temperature_K=[min(temps),max(temps)],max_cfl=max(cfl),max_grid_particle_mismatch=max(mismatch),max_subgrid_increment=max(subgrid),provenance=provenance)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path);p.add_argument('--local',action='store_true');p.add_argument('--summary',type=Path);a=p.parse_args()
    try: result=verify(a.root,a.local)
    except (ValueError,OSError,EOFError,IndexError,KeyError,struct.error) as e:result=dict(status='FAIL',error=str(e))
    text=json.dumps(result,indent=2)+'\n'
    if a.summary:a.summary.write_text(text)
    print(text);raise SystemExit(result['status']!='PASS')
