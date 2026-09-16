"""70 Ma radial-column validation using cfg, nodal Tref and compiled C kernel.
Usage: python3 tests/validate_rheol7_initial_profile.py WORKSPACE_ROOT
"""
from pathlib import Path
import sys, os, tempfile, subprocess, configparser, json
os.environ.setdefault('MPLCONFIGDIR',tempfile.gettempdir()+'/rheol7-mpl')
import numpy as np
from scipy.special import erfc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(sys.argv[1]).resolve();src=Path(__file__).resolve().parents[1]
run=root/'worktrees/runs-cmbhf_EBA';out=root/'output/rheol7_initial_validation';out.mkdir(exist_ok=True)
cfg=Path(sys.argv[2]).resolve() if len(sys.argv)>2 else run/'cmbhf_EBA.cfg'
c=configparser.ConfigParser(interpolation=None);c.read(cfg);p=c['CitcomS.solver.param'];k=c['CitcomS.solver.const'];v=c['CitcomS.solver.visc']
cold_scale=v.getfloat('cold_scale',fallback=0.5)
out=out/f'cold_scale_{cold_scale:g}';out.mkdir(exist_ok=True)
R=k.getfloat('radius');top=k.getfloat('Ttop');bottom=k.getfloat('Tbottom');eta0=k.getfloat('refvisc');delta=bottom-top
r=np.loadtxt(run/'GLB.coor.global.dat',skiprows=1,usecols=1);z=(1-r)*R/1000;z[0]=2891;z[-1]=0
tr=top+delta*np.loadtxt(run/'refstate_EBA.txt')[:,2]
assert p.getfloat('max_plate_age_Ma')==70
age=70.;kappa=k.getfloat('k0')/(k.getfloat('rho0')*k.getfloat('Cp0'));sec=age*1e6*365.25*86400
height=(z[0]-z)*1000;h=p.getfloat('bottom_tbl_thickness')*R
bot=(bottom-tr[0])*erfc(1.8214*height/h)
cold=(tr[-1]-top)*erfc(z*1000/(2*np.sqrt(kappa*sec)))
# Standard column: full lith_age_depth, no trench-dependent half-depth override.
cold[z>p.getfloat('lith_age_depth')*R/1000]=0
T=tr+bot-cold;T[0]=bottom;T[-1]=top
# Same nodal-to-Gauss linear radial interpolation as production.
gz=[];gt=[];gr=[];centers=[]
for i in range(len(z)-1):
 for xi in [-1/np.sqrt(3),1/np.sqrt(3)]:
  w=(1+xi)/2;gz.append((1-w)*z[i]+w*z[i+1]);gt.append((1-w)*T[i]+w*T[i+1]);gr.append((1-w)*tr[i]+w*tr[i+1]);centers.append((z[i]+z[i+1])/2)
gz,gt,gr,centers=map(np.array,(gz,gt,gr,centers))
with tempfile.TemporaryDirectory() as tmp:
 tmp=Path(tmp);(tmp/'main.c').write_text('#include <stdio.h>\n#include "Steinberger_nuref.h"\nint main(void){double z,r,t,e,s;while(scanf("%lf%lf%lf%lf%lf",&z,&r,&t,&e,&s)==5)printf("%.17g %.17g %.17g\\n",steinberger_Az(z),steinberger_nuref(z,r,e),steinberger_viscosity(z,r,t,e,s));return 0;}')
 subprocess.run(['cc','-std=c99','-Wall','-Wextra','-Werror','-I'+str(src/'lib'),str(tmp/'main.c'),'-lm','-o',str(tmp/'main')],check=True)
 def evaluate(z,r,t):
  data=''.join(f'{a:.17g} {b:.17g} {d:.17g} {eta0:.17g} {cold_scale:.17g}\n' for a,b,d in zip(z,r,t))
  return np.array([list(map(float,l.split())) for l in subprocess.check_output([str(tmp/'main')],input=data,text=True).splitlines()])
 node=evaluate(z,tr,T);gp=evaluate(gz,gr,gt)
 # Independent Kelvin formula, identity, signs and cold scaled-exponent tests.
 for zz,rr,tt,res in [(z,tr,T,node),(gz,gr,gt,gp)]:
  expected=res[:,1]*np.exp(np.where(tt<rr,cold_scale,1)*res[:,0]*(1/tt-1/rr))
  assert np.allclose(res[:,2],expected,rtol=2e-13)
  assert np.all(np.isfinite(res)) and np.all(res>0)
  identity=evaluate(zz,rr,rr);assert np.allclose(identity[:,1],identity[:,2],rtol=1e-14)
  probe=evaluate(zz,rr,rr-100)
  assert np.allclose(np.log(probe[:,2]/probe[:,1]),cold_scale*probe[:,0]*(1/(rr-100)-1/rr),atol=1e-13)
  hot=evaluate(zz,rr,rr+100);assert np.all(hot[:,2]<hot[:,1])
 low=v.getfloat('visc_min');high=v.getfloat('visc_max')
 def limited(raw,center):
  upper=np.where(1-center/(R/1000)>.89641,high,5*high)
  return np.maximum(low,np.minimum(upper,raw))
 nlim=limited(node[:,2],z);glim=limited(gp[:,2],centers)
 np.savetxt(out/'nodes_70Ma.txt',np.c_[z,tr,T,node,nlim],header='depth_km Tref_K Tinitial_K Az_K nuref eta_raw_nd eta_limited_nd',fmt='%.12e')
 np.savetxt(out/'gauss_70Ma.txt',np.c_[gz,gr,gt,gp,glim],header='depth_km Tref_K Tinitial_K Az_K nuref eta_raw_nd eta_limited_nd',fmt='%.12e')
 report={'cfg':str(cfg),'cold_scale':cold_scale,'age_Ma':age,'surface_assimilation_depth_km':p.getfloat('lith_age_depth')*R/1000,'bottom_TBL_1percent_height_km':h/1000,'kappa_m2_s':kappa,'surface_temperature_K':T[-1],'cmb_temperature_K':T[0],'surface_raw_Pa_s':node[-1,2]*eta0,'cmb_raw_Pa_s':node[0,2]*eta0,'surface_limited_Pa_s':nlim[-1]*eta0,'cmb_limited_Pa_s':nlim[0]*eta0,'gp_count':len(gz),'gp_min_limited':int(np.sum(gp[:,2]<low)),'gp_max_limited':int(np.sum(gp[:,2]>np.where(1-centers/(R/1000)>.89641,high,5*high))),'checks':'Kelvin formula, reference identity, cold scaled exponent, hot weakening, finite positive, C compilation passed','scope':'1D 70 Ma standard column with full lith_age_depth, not full 3D age map or MPI run; GP caps use element-center radius; nodal caps are illustrative.'}
 (out/'validation.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
 fig,ax=plt.subplots(1,3,figsize=(13,7),layout='constrained')
 ax[0].plot(tr,z,label='Tref');ax[0].plot(T,z,label='70 Ma + bottom TBL');ax[0].set_xlabel('Temperature (K)')
 for a in ax[1:]:
  a.semilogx(node[:,1]*eta0,z,label='nuref');a.semilogx(node[:,2]*eta0,z,label='Temperature-corrected');a.semilogx(nlim*eta0,z,ls='--',label='After viscosity limits')
  a.set_xlabel('Viscosity (Pa s)')
 for a in ax:a.set_ylim(2891,0);a.grid(alpha=.2);a.legend(fontsize=8);a.set_ylabel('Depth (km)')
 ax[2].set_ylim(200,0);ax[2].set_title('Upper 200 km');fig.suptitle(f'rheol=7: Az, cold exponent scale={cold_scale:g}, 70 Ma surface cooling')
 fig.savefig(out/'initial_viscosity_70Ma.png',dpi=150)
