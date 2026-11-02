"""12-rank RH integration: normalization, legacy regression and restart guards."""
import argparse, json, re, shutil, subprocess, struct
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('build',type=Path);p.add_argument('runs',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
a.output=a.output.resolve();a.output.mkdir(parents=True,exist_ok=True)
src=a.runs/'pices_p8a';base=(src/'continuous.cfg').read_text()
def run(name,changes,ok=True,token=None):
 d=a.output/name;d.mkdir();(d/'DATA').mkdir()
 shutil.copy(src/'refstate_EBA_PICES_P4.txt',d)
 shutil.copytree(src/'pices_p4_forcing',d/'pices_p4_forcing',ignore=shutil.ignore_patterns('._*'))
 text=base
 for k,v in changes.items():
  pat=r'^'+re.escape(k)+r'\s*=.*$'
  text=re.sub(pat,k+'='+str(v),text,flags=re.M) if re.search(pat,text,re.M) else text+'\n'+k+'='+str(v)+'\n'
 (d/'case.cfg').write_text(text)
 with (d/'stdout').open('w') as out,(d/'stderr').open('w') as err:
  r=subprocess.run(['mpiexec','--oversubscribe','-n','12',str(a.build.resolve()/'CitcomSFull'),'case.cfg'],cwd=d,stdout=out,stderr=err,timeout=180)
 if ok:assert r.returncode in (0,8),(name,r.returncode,(d/'stderr').read_text()[-3000:])
 else:assert r.returncode not in (0,8) and token in (d/'stderr').read_text(),name
 return d
def payload(path):
 data=bytearray(path.read_bytes())
 if path.name.endswith('.state'):return data
 nx,ny,nz,*_=struct.unpack_from('8i',data)
 nn=nx*ny*nz; ne=(nx-1)*(ny-1)*(nz-1)
 off=44+16
 for _ in range(2):data[off:off+8]=bytes(8);off+=8*(nn+1)
 off+=16+8
 data[off:off+8]=bytes(8);off+=8*(ne+1)+8*3*nn
 off+=16
 basic,extra,flavors,last,n=struct.unpack_from('5i',data,off);off+=20
 assert basic==12
 for _ in range(6+extra):data[off:off+8]=bytes(8);off+=8*(n+1)
 data[off:off+4]=bytes(4);off+=4*(n+1)+16
 nc=struct.unpack_from('i',data,off)[0];off+=4+16*nc
 for _ in range(nc):data[off:off+8]=bytes(8);off+=8*(ne+1)
 assert off==len(data),(off,len(data))
 return data
legacy=run('legacy',{})
one=run('one',{'RH':1})
full=run('twenty',{'RH':20})
for rank in range(12):
 for suffix in ['', '.pices.state']:
  filename='PICES_P4.chkpt.%d.6%s'%(rank,suffix)
  assert payload(legacy/'DATA'/str(rank)/filename)==payload(one/'DATA'/str(rank)/filename),('RH1 regression',rank,suffix)
row=next(l for l in (full/'DATA/0/log').read_text().splitlines() if l.startswith('RH_HEATING'))
r=dict(re.findall(r'(\w+)=([^ ]+)',row))
assert abs(float(r['Qprim'])/float(r['Qbg'])-20)<1e-12
assert abs(float(r['total_W'])/float(r['baseline_W'])-1)<1e-12
lines=(full/'DATA/0/log').read_text().splitlines()
def records(prefix):
 return {int(d['step']):d for l in lines if l.startswith(prefix) for d in [dict(re.findall(r'(\w+)=([^ ]+)',l))]}
heating=records('RH_HEATING ');heat=records('PICES_EBA ');steps=records('PICES_STEP ')
for step,h in heating.items():
 assert all(h[k]==r[k] for k in ('f0','Qbg','Qprim','baseline_W'))
 if step:
  # Mini-case dimensional scale: k0=4, DeltaT=3400 K, R=6371000 m.
  power=float(heat[step]['internal'])/float(steps[step]['dt'])*4*3400*6371000
  assert abs(power/float(h['total_W'])-1)<1e-10,(step,power,h)
split=run('split',{'RH':20,'maxstep':2,'maxtotstep':3})
restart={'RH':20,'restart':'on','solution_cycles_init':2,'datadir_old':str(split/'DATA/%RANK'),'datafile_old':'PICES_P4'}
resumed=run('resumed',restart)
for rank in range(12):
 for suffix in ['', '.pices.state']:
  filename='PICES_P4.chkpt.%d.6%s'%(rank,suffix)
  assert payload(full/'DATA'/str(rank)/filename)==payload(resumed/'DATA'/str(rank)/filename),('restart',rank,suffix)
run('wrong_RH',dict(restart,RH=10),False,'incompatible RH normalization')
# Corrupt authenticated state: must fail before loading coefficients.
f=split/'DATA/0/PICES_P4.chkpt.0.2.pices.state';f.write_bytes(f.read_bytes()[:-8])
run('truncated',restart,False,'checksum mismatch')
run('invalid',{'RH':-1},False,'RH must be')
print(json.dumps({'status':'PASS','initial':r,'checks':['RH1 legacy equality','RH20 power and contrast','continuous/restart byte equality','RH mismatch','truncated state','invalid RH']},indent=2))
