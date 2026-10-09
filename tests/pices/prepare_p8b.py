"""Generate synthetic coupled P8b gates; these are not production inputs."""
import argparse,json,shutil,re
from pathlib import Path
from prepare_p8a import edit

def prepare(runs):
 out=runs/'pices_p8b';out.mkdir(exist_ok=True)
 text=(runs/'pices_p7/cases/assim_pices_rheol7.cfg').read_text()
 chemistry=(runs/'pices_p8a/continuous.cfg').read_text()
 keys=['tracer_flavors','tracer_reclassify_flavors','chemical_buoyancy','buoy_type','buoyancy_ratio','z_interface','ic_method_for_flavors','kC_primordial_flavor','kC_ratio']
 opts={k:re.search(r'^'+k+r'=(.*)$',chemistry,re.M)[1] for k in keys}
 opts.update(maxstep=4,maxtotstep=5,storage_spacing=1,CBF_frequency=1,start_age=2.13,fixed_timestep=0,pices_max_timestep_Ma=.05,file_vbcs=1,remove_rigid_rotation=0,vel_bound_file='forcing/bvel.',qvis_mode=2,qvis_cohesion_pa=1e7,qvis_friction_angle_rad=.085,tracers_per_element=64,output_optional='tracer,comp_nd,comp_el',datafile='PICES_P8b')
 text=edit(text,opts)
 shutil.copytree(runs/'pices_p7/forcing_9',out/'forcing',dirs_exist_ok=True,ignore=shutil.ignore_patterns('._*'))
 # Smooth change in time; constant spherical components are a boundary-reading
 # probe, not a reconstruction of Earth's plates.
 for age in range(4):
  for cap in range(12):
   (out/f'forcing/bvel.{age}.{cap}').write_text((f'0 {0.01*(1+age):.8f}\n')*81)
 lines=(runs/'pices_p7/refstate_assim_9.txt').read_text().splitlines()
 (out/'refstate.txt').write_text(''.join(f'{line} {1.35e11*(8-i)/8:.9g}\n' for i,line in enumerate(lines)))
 cases=[]
 for name,changes,last in [('coupled',{},4),('uncapped',dict(qvis_mode=0),4),('cap_probe',dict(qvis_cohesion_pa=1,qvis_friction_angle_rad=0),4),('static_plate',dict(file_vbcs=0),4),('present_day',dict(start_age=.03),1),('fixed_step',dict(fixed_timestep=3.769596e-8),4),('multigrid',dict(Solver='multigrid',levels=3,mgunitx=2,mgunity=2,mgunitz=2,see_convergence='on',vlowstep=50,maxstep=1,maxtotstep=2),1)]:
  (out/(name+'.cfg')).write_text(edit(text,changes));cases.append(dict(name=name,first=1,last=last))
 (out/'matrix.json').write_text(json.dumps(dict(stage='P8b',cases=cases),indent=2)+'\n')
 (out/'matrix.tsv').write_text(''.join(f"{c['name']}\t1\t{c['last']}\n" for c in cases))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('runs',type=Path);prepare(p.parse_args().runs)
