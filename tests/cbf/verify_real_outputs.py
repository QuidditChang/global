"""Summarize analytic conduction GRDs at 5/9/17-node local test resolutions.
Pass the three run directories in order. This is a sampled analytic field,
not a solved PDE benchmark. Also accepts --mpi-peer for exact field comparison.
"""
import argparse
from pathlib import Path
import numpy as np
from scipy.io import netcdf_file

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('runs',nargs='+',type=Path);p.add_argument('--mpi-peer',type=Path)
a=p.parse_args();previous={}
for run in a.runs:
    fields={}
    for boundary in ['cmbhf','eshf']:
        with netcdf_file(run/'PostProc/HF_CBF'/('%s_CBF_0.grd'%boundary),mmap=False) as f:
            fields[boundary]=(f.variables['z'][:].copy(),float(f.radius_m),float(f.native_integrated_heat_W))
    ri=fields['cmbhf'][1];ro=fields['eshf'][1]
    power=4*np.pi*4*3400*ri*ro/(ro-ri)
    for boundary,(z,r,total) in fields.items():
        exact=power/(4*np.pi*r*r)
        rms=float(np.sqrt(np.mean((z.astype(float)/exact-1)**2)))
        integral=abs(total/power-1)
        assert np.isfinite(z).all() and np.min(z)>0
        assert np.array_equal(z[:,0],z[:,-1])
        assert np.ptp(z[0])==0 and np.ptp(z[-1])==0
        if boundary in previous:assert rms<previous[boundary]
        previous[boundary]=rms
        print(run.name,boundary,'RMS_relative',rms,'integral_relative',integral)
        if a.mpi_peer and run==a.runs[0]:
            with netcdf_file(a.mpi_peer/'PostProc/HF_CBF'/('%s_CBF_0.grd'%boundary),mmap=False) as f:
                np.testing.assert_array_equal(z,f.variables['z'][:])
