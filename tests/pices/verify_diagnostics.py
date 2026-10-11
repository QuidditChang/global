"""Validate actual text/NPZ snapshots, volume integrals, closure and particles.
Run after a PICES_TEST_DIAGNOSTICS=1 manufactured case. No solver internals
are used to calculate the particle statistics or profile volume integrals.
"""
import argparse
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('case',type=Path);p.add_argument('--ranks',type=int,default=12);p.add_argument('--step',type=int,default=1);p.add_argument('--coupled',action='store_true');a=p.parse_args()
f=np.load(a.case/f'DATA/0/profiles_hist_{a.step}.npz')
assert f['pices_diagnostics_valid']==1
line=next(x for x in (a.case/'DATA/0/log').read_text().splitlines() if x.startswith(f'PICES_DIAGNOSTICS step={a.step} '))
values=dict(x.split('=',1) for x in line.split()[1:])
for key in f.files:
 if key.startswith('pices_') and key[6:] in values:
  assert float(f[key])==float(values[key[6:]]),key
v=lambda k:float(f['pices_'+k])
scale=max(abs(v('energy_before_nd')),abs(v('energy_after_nd')),1e-30)
assert abs(v('closure_residual_nd'))<2e-12*scale
assert abs(v('heat_residual_nd'))<2e-12*scale
radius=6371.
for q,total in [('qinternal','internal_input_nd'),('qadi','adiabatic_sink_nd'),('qvisc','viscous_input_nd'),('qphase','phase_sink_nd'),('qassim','ta_storage_nd')]:
 integral=np.dot(f[q+'_stats'][:,1],f[q+'_depth_volume_km3'])/radius**3*v('dt_nd')
 assert abs(integral-v(total))<2e-12*max(abs(v(total)),1), (q,integral,v(total))
expected_total=f['qinternal_stats'][:,1]+f['qvisc_stats'][:,1]-f['qadi_stats'][:,1]-f['qphase_stats'][:,1]+f['qassim_stats'][:,1]
assert np.allclose(f['qtotal_stats'][:,1],expected_total,rtol=1e-12,atol=1e-12)
assert np.allclose(f['qadi_stats'],f['qadi_base_stats'],rtol=0,atol=0)
if not a.coupled:
 particles=np.concatenate([np.loadtxt(a.case/f'particles.{r}.{a.step}.txt')[:,0] for r in range(a.ranks)])
 kelvin=float(f['T_surface_K'])+float(f['DeltaT'])*particles
 assert len(kelvin)==v('particle_count_global')
 for key,actual in [('particle_Tmin_K',kelvin.min()),('particle_Tmax_K',kelvin.max()),('particle_Tmean_count_K',kelvin.mean())]:
  assert abs(actual-v(key))<1e-10,(key,actual,v(key))
initial=np.load(a.case/'DATA/0/profiles_hist_0.npz')
assert initial['pices_diagnostics_valid']==0 and 'qassim_stats' not in initial.files
print('PASS: shared log/NPZ, source integrals, full-step closure, global particle statistics, initial validity')
