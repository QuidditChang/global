"""Compile and exercise production CBF GRD output with a closed synthetic mesh.

Requires mpicc, a matching mpirun, nc-config, scipy and numpy. MPI needs local
socket permission. MPIRUN may specify an absolute matching launcher path.
"""
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile
import unittest
import numpy as np
from scipy.io import netcdf_file

ROOT=Path(__file__).resolve().parents[2]

class GridOutput(unittest.TestCase):
    def test_grid_and_mpi(self):
        with tempfile.TemporaryDirectory(prefix='cbf-grid-') as tmp:
            tmp=Path(tmp);exe=tmp/'fixture'
            flags=shlex.split(subprocess.check_output(['nc-config','--cflags'],text=True))
            libs=shlex.split(subprocess.check_output(['nc-config','--libs'],text=True))
            subprocess.run([os.environ.get('MPICC','mpicc'),'-std=gnu99','-w','-DUSE_CBF_NETCDF',
                '-I'+str(ROOT/'lib'),*flags,str(ROOT/'tests/cbf/grid_fixture.c'),*libs,'-lm','-o',str(exe)],check=True)
            results={}
            for ranks in [1,2]:
                wd=tmp/str(ranks);wd.mkdir()
                subprocess.run([os.environ.get('MPIRUN','mpirun'),'--oversubscribe','-n',str(ranks),str(exe)],cwd=wd,check=True,timeout=60)
                for boundary,mean,slope,total in [('eshf',3,.1,72e6),('cmbhf',7,-.2,42e6)]:
                    with netcdf_file(wd/'PostProc/HF_CBF'/('%s_CBF_50.grd'%boundary),mmap=False) as f:
                        z=f.variables['z'][:].copy()
                        lon,lat=np.meshgrid(np.deg2rad(f.variables['x'][:]),np.deg2rad(f.variables['y'][:]))
                        ray=np.array([np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)])
                        exact=mean+slope*np.sum(ray,axis=0)/np.max(abs(ray),axis=0)
                        self.assertEqual(z.shape,(361,721))
                        self.assertLess(np.max(abs(z-exact)),5e-7)
                        np.testing.assert_array_equal(z[:,0],z[:,-1])
                        self.assertEqual(np.ptp(z[0]),0);self.assertEqual(np.ptp(z[-1]),0)
                        self.assertAlmostEqual(float(f.native_integrated_heat_W),total,places=5)
                        self.assertEqual(f.method,b'CBF_GLL_Q1')
                        self.assertTrue(np.isfinite(f.time_seconds))
                        results[ranks,boundary]=z
            for boundary in ['eshf','cmbhf']:
                np.testing.assert_array_equal(results[1,boundary],results[2,boundary])

if __name__=='__main__':unittest.main()
