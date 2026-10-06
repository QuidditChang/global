"""P7 targeted gates; synthetic inputs are not a full production configuration."""
import argparse,json,re,shutil
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('runs',type=Path);a=p.parse_args();r=a.runs;src=r/'pices_p6';out=r/'pices_p7';(out/'cases').mkdir(parents=True,exist_ok=True)
def edit(text,opts):
 for k,v in opts.items():
  pattern=r'^'+re.escape(k)+r'=.*$'
  text=re.sub(pattern,k+'='+str(v),text,flags=re.M) if re.search(pattern,text,re.M) else text+'\n'+k+'='+str(v)+'\n'
 return text
rows=[]
for method in ['pg','pices']:
 for variant,factor in [('base',1),('dt_half',2),('dt_quarter',4),('dt_eighth',8)]:
  row=dict(name=f'assim_{method}_{variant}',scenario='assim',method=method,variant=variant,nodes=9,particles=32 if method=='pices' else 0,steps=4*factor,dt=1e-7/factor,prescribed=False,length_scale=1)
  text=(src/'cases'/f'assim_{method}_base.cfg').read_text();text=edit(text,dict(maxstep=row['steps'],maxtotstep=row['steps']+1,fixed_timestep=format(row['dt'],'.17g'),storage_spacing=row['steps'],CBF_frequency=row['steps']))
  (out/'cases'/f"{row['name']}.cfg").write_text(text);rows.append(row)
# Exercise target rheology on the same synthetic mesh/forcing before enabling
# its checkpoint support or attempting the full-size production mesh.
for method in ['pg','pices']:
 row=dict(name=f'assim_{method}_rheol7',scenario='assim',method=method,variant='rheol7',nodes=9,particles=32 if method=='pices' else 0,steps=8,dt=1e-9,prescribed=False,length_scale=1)
 text=(src/'cases'/f'assim_{method}_base.cfg').read_text()
 text=edit(text,dict(topvbc=1,topvbxval=0,topvbyval=0,Solver='cgrad',levels=1,mgunitx=8,mgunity=8,mgunitz=8,vlowstep=2000,rheol=7,cold_scale=1.0,TDEPV='on',visc0='0.1,0.1,1.0,30',viscE='15,15,15,15',viscT='0.1,0.1,0.1,0.1',VMIN='on',VMAX='on',visc_min=.1,visc_max=100,maxstep=8,maxtotstep=9,fixed_timestep='1e-9',storage_spacing=8,CBF_frequency=8))
 (out/'cases'/f"{row['name']}.cfg").write_text(text);rows.append(row)
row=dict(name='sharp_pices_base',scenario='sharp',method='pices',variant='base',nodes=9,particles=32,steps=1,dt=1e-10,prescribed=True,length_scale=1);rows.append(row)
text=edit((src/'cases/cold_pices_base.cfg').read_text(),dict(p5_case='sharp',p5_omega=0,maxstep=1,maxtotstep=2,fixed_timestep='1e-10',storage_spacing=1,CBF_frequency=1))
(out/'cases/sharp_pices_base.cfg').write_text(text)
for kind in ['assim','transport']:shutil.copy(src/f'refstate_{kind}_9.txt',out)
shutil.copytree(src/'forcing_9',out/'forcing_9',dirs_exist_ok=True,ignore=shutil.ignore_patterns('._*'))
shutil.copytree(src/'restart',out/'restart',dirs_exist_ok=True,ignore=shutil.ignore_patterns('._*'))
for suffix,last in [('',64),('_split',32),('_restart',64)]:
 f=out/'restart'/f'cmbhf_EBA_PICES_P4{suffix}.cfg';opts=dict(maxstep=last,maxtotstep=last+1,fixed_timestep="2.5e-8")
 if suffix=='_restart':opts['solution_cycles_init']=32
 f.write_text(edit(f.read_text(),opts))
(out/'matrix.json').write_text(json.dumps(dict(schema=1,stage='P7',cases=rows,production_switch=False,restart_final=64,restart_split=32),indent=2)+'\n')
(out/'matrix.tsv').write_text(''.join(f"{x['name']}\t{x['nodes']}\t{x['steps']}\t{'assim' if x['scenario']=='assim' else 'transport'}\n" for x in rows))
print('11 targeted cases + 64/32/restart-to-64 endurance')
# Pin the user's exact production target without pretending these synthetic
# inputs reproduce its unsupported chemistry, heating and plate boundaries.
import hashlib
f=r/'cmbhf_EBA_Q0_30_rheol7_scold1.0_LLSVPsV0.04_B24_0.3_TA.cfg'
target=dict(target=f.name,target_sha256=hashlib.sha256(f.read_bytes()).hexdigest(),production_compatible=False,mpi_ranks=384,blockers=['chemical_buoyancy, 25 tracer flavors and flavor reclassification','kC_ratio=0.8 (PICES requires 1)','qvis_mode=2 (PICES requires 0)','file_vbcs=1 (PICES rejects prescribed plate files)','rheol7 restart (checkpoint currently requires uniform constant viscosity)'],note='P7 rheology probe matches cold_scale=1 and fixed top BC; uses synthetic P6 reference state/TA and single-level CG, not the full production setup.')
shutil.copyfile(f,out/'production_target.cfg') # Reference only; never in matrix.tsv.
(out/'production_target.json').write_text(json.dumps(target,indent=2)+'\n')
