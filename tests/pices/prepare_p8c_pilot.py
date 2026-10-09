"""Derive the real Pyre pilot without changing the target model's physics."""
import argparse,json,re
from pathlib import Path
TARGET='cmbhf_EBA_Q0_30_rheol7_scold1.0_LLSVPsV0.04_B24_0.3_TA'
def prepare(runs):
    original=(runs/(TARGET+'.cfg')).read_text();text=original
    changes={('CitcomS','steps'):'2',('CitcomS.job','name'):'cmbhf_EBA_PICES_P8c_pilot',('CitcomS.job','queue'):'medium',('CitcomS.controller','monitoringFrequency'):'1',('CitcomS.controller','profileMonitoringFrequency'):'1',('CitcomS.controller','checkpointFrequency'):'2',('CitcomS.solver','CBF_use_advection'):'off'}
    added=dict(energy_solver='pices',pices_projection='bounded_consistent',pices_eba='on',pices_p4='on',pices_checkpoint='on',pices_max_timestep_Ma='0.1')
    section='';lines=[];found=set()
    for line in original.splitlines(True):
        m=re.match(r'\[([^]]+)\]',line)
        if m:section=m[1]
        m=re.match(r'(\w+)\s*=\s*(.*?)\s*$',line)
        if m and (section,m[1]) in changes:
            key=(section,m[1]);line=m[1]+' = '+changes[key]+'\n';found.add(key)
        lines.append(line)
        if line.strip()=='[CitcomS.solver.output]':lines.append('output_format = ascii-gz\n')
        if line.strip()=='[CitcomS.solver.tsolver]':lines.extend(k+' = '+v+'\n' for k,v in added.items())
    assert found==set(changes)
    text=''.join(lines)
    (runs/'cmbhf_EBA_PICES_P8c_pilot.cfg').write_text(text)
    def parse(t):
        section='';out={}
        for line in t.splitlines():
            line=line.split('#',1)[0].strip()
            if line.startswith('['):section=line[1:-1]
            elif '=' in line:
                k,v=map(str.strip,line.split('=',1));assert (section,k) not in out;out[section,k]=v
        return out
    before,after=parse(original),parse(text)
    delta=[dict(section=s,key=k,before=before.get((s,k)),after=after.get((s,k))) for s,k in sorted(set(before)|set(after)) if before.get((s,k))!=after.get((s,k))]
    assert len(delta)==len(changes)+len(added)+1
    (runs/'PICES_P8C_PILOT_CONFIG_DIFF.json').write_text(json.dumps(dict(target=TARGET,solver_ranks=384,nodes_at_40_slots=10,mesh=[129,129,65],particles_per_element=27,initial_particles=12*128*128*64*27,changes=delta),indent=2)+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('runs',type=Path);prepare(p.parse_args().runs)
