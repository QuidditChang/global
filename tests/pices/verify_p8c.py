"""P8c: accepted coupled state, forcing-knot crossing and exact restart audit."""
import argparse,gzip,hashlib,json,math,re,struct
from pathlib import Path
from verify_p2 import checkpoint,require,fields,sha,config
from verify_p4 import boundary

def verify(root,local=False):
    rows=json.loads((root/'input/matrix.json').read_text())['cases']
    require([(r['name'],r['first'],r['last']) for r in rows]==[('continuous',1,4),('split',1,2),('restart',3,4)],'matrix')
    states={};summary={};compared=0
    for row in rows:
        name=row['name'];d=root/name;first,last=row['first'],row['last'];cfg=config(d/'case.cfg')
        require((d/'case.cfg').read_bytes()==(root/'input'/(name+'.cfg')).read_bytes(),'cfg identity')
        for f in (root/'input').rglob('*'):
            if f.is_file() and (f.name=='refstate.txt' or 'forcing' in f.parts) and not f.name.startswith('._'):
                require((d/f.relative_to(root/'input')).read_bytes()==f.read_bytes(),'forcing/reference identity')
        for k,v in dict(p5_case='off',Solver='multigrid',rheol='7',qvis_mode='2',file_vbcs='1',fixed_timestep='0',pices_checkpoint='on',tracer_flavors='25').items():require(cfg[k]==v,'coupled cfg '+k)
        require((d/'mpi_exit_code.txt').read_text().strip() in ('0','8'),'MPI exit')
        require('PICES_ERROR' not in (d/'solver.stderr').read_text(),'runtime error')
        cfl=[];heat_error=0;count_by_step={};primordial={}
        for rank in range(12):
            folder=d/'DATA'/str(rank);log=(folder/'log').read_text().splitlines()
            require(sum(x.startswith('PICES_RESTORE ') for x in log)==int(name=='restart'),'restore lifecycle')
            records=[fields(x) for x in log if x.startswith('PICES_STEP ')]
            require([int(x['step']) for x in records]==list(range(first,last+1)),'step sequence')
            for x in records:
                for k in ('time','dt','cfl','Tmin','Tmax','mismatch'):require(math.isfinite(float(x[k])),'finite step')
                cfl.append(float(x['cfl']));require(0<=cfl[-1]<=.25,'CFL')
            ta=[fields(x) for x in log if x.startswith('PICES_TA ')]
            require([int(x['step']) for x in ta]==list(range(first,last+1)) and all(x['calls']=='1' for x in ta),'TA exactly once')
            for x,y in zip(ta,records):require(float(x['time'])==float(y['time']) and float(x['dt'])==float(y['dt']),'TA accepted clock')
            heat=[fields(x) for x in log if x.startswith('PICES_EBA ')]
            require(len(heat)==last-first+1,'heat records')
            for x in heat:
                v=[float(x[k]) for k in ('storage','internal','adiabatic','viscous','phase_pressure','boundary_reaction')]
                require(all(map(math.isfinite,v)),'finite heat');err=abs(v[0]-sum(v[1:]));heat_error=max(heat_error,err)
                require(err<1e-12*max(1,sum(map(abs,v))),'heat balance')
            for x in [fields(x) for x in log if x.startswith('PICES_PROJECTION ')]:require(0<=float(x['residual'])<=float(x['tolerance']),'projection')
            if rank==0:
                stokes=[fields(x) for x in log if x.startswith('PICES_STOKES ')]
                require({int(x['step']) for x in stokes}>=set(range(first,last+1)) and all(x['status']=='PASS' for x in stokes),'Stokes')
                mg=[fields(x) for x in log if x.startswith('PICES_MG_RESIDUAL ')]
                require(bool(mg) and all(x['status']=='PASS' and 0<=float(x['actual'])<=float(x['tolerance']) for x in mg),'MG true residual')
            for step in range(first-1,last+1,2):
                p=folder/f'PICES_P8c.chkpt.{rank}.{step}';m,live=checkpoint(p)
                require(m['schema']==4 and m['step']==step and m['rank']==rank,'schema/rank/step')
                states[name,rank,step]=live
                count_by_step[step]=count_by_step.get(step,0)+m['particles']
                colors=[int(x[0]) for x in struct.iter_unpack('=d',live[11])]
                require(all(0<=c<25 and c not in (18,19) for c in colors),'flavors')
                primordial[step]=primordial.get(step,0)+colors.count(24)
                if not local:require(m['solver_commit']==(root/'solver_commit.txt').read_text().strip(),'build identity')
        require(set(count_by_step.values())=={393216},'particle conservation')
        require(len(set(primordial.values()))==1 and next(iter(primordial.values()))>0,'primordial conservation')
        for step in range(first-1,last+1,2):
            paths=[d/f'DATA/{r}/PICES_P8c.chkpt.{r}.{step}' for r in range(12)]
            digest=hashlib.sha256(b''.join(sha(Path(str(p)+'.pices.json')).encode()+b'\0' for p in paths)).hexdigest()
            for p in paths:require(json.loads(Path(str(p)+'.pices.manifest').read_text())==dict(magic='CITCOMS_EBA_PICES_COMPLETE',schema=1,mpi_size=12,metadata_set_sha256=digest),'collective manifest')
        heat={int(fields(x)['step']):fields(x) for x in (d/'DATA/0/log').read_text().splitlines() if x.startswith('PICES_EBA ')}
        steps={int(fields(x)['step']):fields(x) for x in (d/'DATA/0/log').read_text().splitlines() if x.startswith('PICES_STEP ')}
        for step in range(first,last+1):
            surf=boundary(d,step,'surf');botm=boundary(d,step,'botm')
            reaction=float(heat[step]['boundary_reaction'])/float(steps[step]['dt'])*float(cfg['k0'])*(float(cfg['Tbottom'])-float(cfg['Ttop']))*float(cfg['radius'])
            require(abs(botm-surf-reaction)<1e-10*max(1,abs(botm)+abs(surf)+abs(reaction)),'CBF balance')
        summary[name]=dict(max_cfl=max(cfl),max_heat_error=heat_error,particles=count_by_step,primordial=primordial)
    for name,start,end in [('split',0,2),('restart',2,4)]:
        for rank in range(12):
            for step in range(start,end+1):
                if step>0:
                    for side in ('surf','botm'):
                        q=Path('DATA')/str(rank)/f'q.{side}.{rank}.{step}'
                        require((root/'continuous'/q).read_bytes()==(root/name/q).read_bytes(),'CBF restart equality')
                rel=Path('DATA')/str(rank)/str(step)
                files={f.name for f in (root/'continuous'/rel).glob('*.gz') if not f.name.startswith('._')}
                require(bool(files) and files=={f.name for f in (root/name/rel).glob('*.gz') if not f.name.startswith('._')},'output inventory')
                for f in files:
                    raw=gzip.decompress((root/'continuous'/rel/f).read_bytes())
                    require(all(math.isfinite(float(v)) for v in raw.split()),'finite output '+str(rel/f))
                    require(raw==gzip.decompress((root/name/rel/f).read_bytes()),'restart field '+str(rel/f));compared+=1
            for step in range(start,end+1,2):require(states['continuous',rank,step]==states[name,rank,step],'accepted binary state')
    # A step must land near 2 Ma and the final state lie below it.
    cfg=config(root/'continuous/case.cfg');log=(root/'continuous/DATA/0/log').read_text()
    ages=[float(fields(x)['age_Ma']) for x in log.splitlines() if x.startswith('PICES_TIMESTEP ')]
    scalet=6371000.**2/(4/(3300*1250))/(1e6*365.25*86400)
    final=checkpoint(root/'continuous/DATA/0/PICES_P8c.chkpt.0.4')[0]
    require(any(abs(a-2)<2e-6 for a in ages) and float(cfg['start_age'])-final['time']*scalet<2,'forcing-knot crossing')
    if not local:
        require((root/'launcher_exit_code.txt').read_text().strip()=='0','launcher')
        require((root/'completed_cases.txt').read_text().strip()=='3','completed cases')
        require((root/'runtime_source_diff.txt').read_text()=='','runtime source diff')
        commit=(root/'solver_commit.txt').read_text().strip()
        require(re.fullmatch('[a-f0-9]{40}',commit) and commit==(root/'pices-p0-build.commit').read_text().strip(),'build commit')
        require('P0_BUILD_COMPLETE commit='+commit in (root/'pices-p0-build.log').read_text(),'build completion')
        hashed={}
        for line in (root/'input.sha256').read_text().splitlines():
            digest,file=line.split();require(file.startswith('input/') and '..' not in Path(file).parts and file not in hashed,'hash path');hashed[file]=digest;require(sha(root/file)==digest,'input digest')
        require(set(hashed)=={str(f.relative_to(root)) for f in (root/'input').rglob('*') if f.is_file()},'input inventory')
        for f in ('submitted.lsf','runs_commit.txt','binary.sha256','platform.txt','mpi_version.txt'):require((root/f).stat().st_size>0,'provenance '+f)
    return dict(status='PASS',cases=summary,restart_outputs_compared=compared,accepted_binary_state_equal=True,production_release='HOLD: actual-data pilot and P7 accuracy investigation remain')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('--local',action='store_true');p.add_argument('--summary',type=Path);a=p.parse_args()
    try:r=verify(a.root,a.local)
    except (OSError,ValueError,KeyError,IndexError,struct.error) as e:r=dict(status='FAIL',error=str(e))
    text=json.dumps(r,indent=2)+'\n'
    if a.summary:a.summary.write_text(text)
    print(text);raise SystemExit(r['status']!='PASS')
