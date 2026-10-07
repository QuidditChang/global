"""Small production-composition gate using the established P4 forcing/physics."""
import argparse,re,shutil,json
from pathlib import Path

def edit(text,opts):
 for k,v in opts.items():
  pat=r'^'+re.escape(k)+r'=.*$'
  text=re.sub(pat,k+'='+str(v),text,flags=re.M) if re.search(pat,text,re.M) else text+'\n'+k+'='+str(v)+'\n'
 return text

def prepare(runs):
 src=runs/'pices_p7/restart';out=runs/'pices_p8a';out.mkdir(exist_ok=True)
 text=(src/'cmbhf_EBA_PICES_P4.cfg').read_text()
 target=(runs/'cmbhf_EBA_Q0_30_rheol7_scold1.0_LLSVPsV0.04_B24_0.3_TA.cfg').read_text()
 def target_value(k):return ','.join(v.strip() for v in re.search(r'^'+k+r'\s*=\s*([^\n#]+)',target,re.M)[1].split(','))
 common=dict(maxstep=6,maxtotstep=7,start_age=2.1,tracer_flavors=25,chemical_buoyancy='on',buoy_type=1,tracer_reclassify_flavors='on',kC_ratio=.8,kC_primordial_flavor=24,ic_method_for_flavors=0,buoyancy_ratio=target_value('buoyancy_ratio'),z_interface=target_value('z_interface'),output_optional='tracer,comp_nd,comp_el')
 text=edit(text,common)
 rows=[]
 for name,changes,first,last in [('continuous',{},1,6),('split',dict(maxstep=2,maxtotstep=3),1,2),('restart',dict(restart='on',solution_cycles_init=2,datadir_old='../split/DATA/%RANK',datafile_old='PICES_P4'),3,6),('conductivity_control',dict(kC_ratio=1),1,6),('buoyancy_control',dict(buoyancy_ratio=','.join(['0']*24)),1,6)]:
  (out/(name+'.cfg')).write_text(edit(text,changes));rows.append(dict(name=name,first=first,last=last))
 shutil.copy(src/'refstate_EBA_PICES_P4.txt',out)
 shutil.copytree(src/'pices_p4_forcing',out/'pices_p4_forcing',dirs_exist_ok=True,ignore=shutil.ignore_patterns('._*'))
 (out/'matrix.json').write_text(json.dumps(dict(stage='P8a',cases=rows),indent=2)+'\n')
 (out/'matrix.tsv').write_text(''.join(f"{r['name']}\t{r['first']}\t{r['last']}\n" for r in rows))
 print(out)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('runs',type=Path);prepare(p.parse_args().runs)
