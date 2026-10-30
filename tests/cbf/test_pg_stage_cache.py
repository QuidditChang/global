"""Exercise the production cache lifecycle with mutable thermal/source state.
The physical Galerkin operator itself is covered by test_cbf_kernel.py.
"""
from pathlib import Path
import subprocess
import tempfile
import unittest
from test_cbf_kernel import function
ROOT=Path(__file__).resolve().parents[2]

class PGStageCache(unittest.TestCase):
    def test_frozen_endpoint_and_initial_unavailable(self):
        source=(ROOT/'lib/Process_buoyancy.c').read_text()
        parts='\n'.join(function(source,s) for s in (
            'void CBF_capture_pg_endpoint(', 'static void evaluate_heat_flux_CBF(',
            'void CBF_boundary_stats('))
        prefix=r'''
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <float.h>
#include "element_definitions.h"
#include "global_defs.h"
static int reduce(const void *a,void *b,int n,MPI_Datatype t,MPI_Op o,MPI_Comm c) {
 memcpy(b,a,n*(t==MPI_INT?sizeof(int):sizeof(double)));return 0;
}
#define MPI_Allreduce reduce
void parallel_process_termination(void) {abort();}
void CBF_heat_sources(struct All_variables *E,int m,double *a,double *v) {
 a[1]=2*E->T[m][1]; v[1]=3*E->T[m][1];
}
double CBF_element_thermal_residual(struct All_variables *E,int m,int e,double r[9]) {
 int a;for(a=1;a<=8;a++)r[a]=E->T[m][1]+E->Tdot[m][1]+E->heating_visc[m][e]-E->heating_adi[m][e];return 0;
}
static void heat_flux_CBF_boundary(struct All_variables *E,int top,double *slice[NCS],int write,double stats[4]) {
 int i;if(stats)for(i=0;i<4;i++)stats[i]=E->pg_cbf_residual[1][8];
}
'''
        main=r'''
int main(void) {
 struct All_variables *E=calloc(1,sizeof(*E)); double t[2]={0,2},dt[2]={0,5};
 double adi[2]={0,101},visc[2]={0,202},s[2][4];int i,j;
 E->parallel.me=1;E->sphere.caps_per_proc=1; E->lmesh.nel=E->lmesh.elz=1;
 E->mesh.toptbc=E->mesh.bottbc=1;E->T[1]=t;E->Tdot[1]=dt;
 E->heating_adi[1]=adi;E->heating_visc[1]=visc;
 CBF_boundary_stats(E,s);for(i=0;i<2;i++)for(j=0;j<4;j++)assert(isnan(s[i][j]));
 CBF_capture_pg_endpoint(E);assert(E->pg_cbf_valid);
 assert(E->heating_adi[1]==adi && E->heating_visc[1]==visc);
 assert(adi[1]==101 && visc[1]==202 && t[1]==2 && dt[1]==5);
 /* Subsequent TA and source/velocity updates must not enter the saved CBF. */
 t[1]=999;dt[1]=-100;adi[1]=500;visc[1]=-500;
 CBF_boundary_stats(E,s);assert(s[0][0]==9 && s[1][0]==9);
 E->pg_cbf_valid=0;CBF_boundary_stats(E,s);assert(isnan(s[0][0]));
 t[1]=3;dt[1]=4;CBF_capture_pg_endpoint(E);CBF_boundary_stats(E,s);assert(s[0][0]==10);
 free(E->pg_cbf_residual[1]);free(E);return 0;
}
'''
        with tempfile.TemporaryDirectory() as d:
            c=Path(d)/'cache.c';c.write_text(prefix+parts+main)
            exe=Path(d)/'cache'
            subprocess.run(['mpicc','-std=gnu99','-I'+str(ROOT/'lib'),str(c),'-lm','-o',str(exe)],check=True)
            subprocess.run([str(exe)],check=True)

if __name__=='__main__':unittest.main()
