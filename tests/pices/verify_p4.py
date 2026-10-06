"""Independent P4 native-face integration, TA lifecycle and restart audit."""
import math
from pathlib import Path
from verify_p2 import require, config, fields

def boundary(root,step,side,expected_state="pices"):
    total=0.; area=0.; declared=[]
    for rank in range(12):
        p=root/'DATA'/str(rank)/f'q.{side}.{rank}.{step}'
        lines=p.read_text().splitlines();nodes={}
        state='state=PICES_heat_stage_average_before_TA derivative=material_heat_only advection=particles TA=excluded' if expected_state=='pices' else 'state=output_T_and_solver_Tdot'
        require(state in lines[3],'CBF stage semantics')
        scale=fields(lines[1]);length=float(scale['length_scale_m']);qscale=float(scale['k0_W_m_K'])*float(scale['deltaT_K'])/length
        hdr=fields(lines[2]);declared.append((float(hdr['global_heat_W']),float(hdr['global_area_m2'])))
        for line in lines:
            v=line.split()
            if v[0]=='N':
                q,rhs,mass=map(float,v[-3:]);require(all(map(math.isfinite,(q,rhs,mass))) and mass>0,'CBF finite positive area')
                expected=(1 if side=='surf' else -1)*rhs*qscale/mass
                require(math.isclose(q,expected,rel_tol=1e-12,abs_tol=1e-15),'CBF nodal sign/scale')
                nodes[int(v[2])]=q
            elif v[0]=='F':
                ns=list(map(int,v[3:7]));weights=list(map(float,v[7:11]));require(len(weights)==4 and min(weights)>0,'CBF face weights')
                total+=sum(nodes[n]*w for n,w in zip(ns,weights))*length**2;area+=sum(weights)*length**2
    for q,a in declared:
        require(math.isclose(total,q,rel_tol=1e-12,abs_tol=1e-4),'CBF independent face integral')
        require(math.isclose(area,a,rel_tol=1e-12),'CBF area integral')
    return total

def verify_p4(root,local=False,final_step=4,split_step=2):
    root=Path(root);max_balance=0.;max_delta=0.;max_mapping=0.;compared=0
    for name,first,last in [('continuous',1,final_step),('split',1,split_step),('restart',split_step+1,final_step)]:
        d=root/name;cfg=config(d/'case.cfg')
        for k,v in dict(pices_p4='on',lith_age='1',lith_age_asml='1',lith_age_time='1',temperature_bound_adj='0',output_q_surf_CBF='on',output_q_botm_CBF='on',CBF_use_advection='off').items():require(cfg.get(k)==v,'P4 cfg '+k)
        for rank in range(12):
            log=(d/'DATA'/str(rank)/'log').read_text().splitlines();ta=[fields(l) for l in log if l.startswith('PICES_TA ')]
            require([int(x['step']) for x in ta]==list(range(first,last+1)),'TA exactly once per accepted step, never at restart')
            for row in ta:
                require(row['calls']=='1','TA call count')
                vals=[float(row[k]) for k in ['time','dt','delta_max','mapping_error','storage']]
                require(all(map(math.isfinite,vals)) and vals[2]>0 and vals[3]>=0,'TA active finite ledger')
                require(math.isclose(vals[0],int(row['step'])*vals[1],rel_tol=max(2e-7,final_step*6e-8)),'TA accepted end time')
                max_delta=max(max_delta,vals[2]);max_mapping=max(max_mapping,vals[3])
        log=(d/'DATA/0/log').read_text().splitlines()
        heat={int(fields(l)['step']):fields(l) for l in log if l.startswith('PICES_EBA ')}
        steps={int(fields(l)['step']):fields(l) for l in log if l.startswith('PICES_STEP ')}
        scale=float(cfg['k0'])*(float(cfg['Tbottom'])-float(cfg['Ttop']))*float(cfg['radius'])
        for step in range(first,last+1):
            surf=boundary(d,step,'surf');botm=boundary(d,step,'botm')
            reaction=float(heat[step]['boundary_reaction'])/float(steps[step]['dt'])*scale
            residual=abs(botm-surf-reaction)/max(1.,abs(surf)+abs(botm)+abs(reaction));max_balance=max(max_balance,residual)
            require(residual<1e-10,'CBF material heat-stage boundary ledger')
        if name!='restart':require(not list((d/'DATA').glob('*/q.*.*.0')),'fabricated initial CBF')
        if not local:
            for f in (root/'pices_p4_forcing').iterdir():require((d/'pices_p4_forcing'/f.name).read_bytes()==f.read_bytes(),'actual forcing differs from archived input')
    for name,steps in [('split',range(1,split_step+1)),('restart',range(split_step,final_step+1))]:
        for rank in range(12):
            for step in steps:
                for side in ['surf','botm']:
                    rel=Path('DATA')/str(rank)/f'q.{side}.{rank}.{step}'
                    require((root/'continuous'/rel).read_bytes()==(root/name/rel).read_bytes(),'CBF restart mismatch '+str(rel));compared+=1
    return dict(status='PASS',cbf_outputs_compared=compared,max_relative_heat_boundary_residual=max_balance,max_TA_delta=max_delta,max_TA_mapping_error=max_mapping)
