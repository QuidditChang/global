"""Compile actual CFL functions and exercise serial and two-rank timestep bounds."""
from pathlib import Path
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
s = (ROOT / 'lib/Pices.c').read_text()
assert 'p->cfl_edge[e]=fmin(p->cfl_edge[e],sqrt(l2));' in s
assert 'pices_cfl' in (ROOT/'module/setProperties.c').read_text()
assert 'pices_cfl = prop.float("pices_cfl", default=0.24)' in (ROOT/'CitcomS/Components/Advection_diffusion/Advection_diffusion.py').read_text()
functions = s[s.index('static double particle_cfl_rate('):s.index('/* Read-only operator-stage minima.')]
stub = r'''
#include <mpi.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
struct All_variables {
 struct { int nel; } lmesh;
 struct { int node[9]; } ien[2][3];
 struct { struct { double V[4][17]; } cap[2]; } sphere;
 struct { double *cfl_edge,min_edge,cfl,max_timestep_Ma; int initialized; } pices;
 struct { float timestep; } advection;
 struct { double scalet; } data;
 struct { int lith_age,vbcs_file; double start_age; } control;
 struct { double elapsed_time; int solution_cycles; } monitor;
 FILE *fp;
};
static double reduce(struct All_variables *E,double x,MPI_Op op){double y;MPI_Allreduce(&x,&y,1,MPI_DOUBLE,op,MPI_COMM_WORLD);return y;}
static void pices_fail(struct All_variables *E,const char *message){fprintf(stderr,"%s\n",message);MPI_Abort(MPI_COMM_WORLD,72);exit(72);}
static int pices_time_finished(struct All_variables *E){return E->control.start_age-E->monitor.elapsed_time*E->data.scalet<=0;}
'''
main = r'''
int main(int argc,char **argv){
 MPI_Init(&argc,&argv);int rank,size;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
 struct All_variables E={0};double h[3]={0,1,10},legacy,cfl;E.fp=stdout;E.lmesh.nel=2;
 E.pices.cfl_edge=h;E.pices.min_edge=1;E.pices.cfl=.24;E.pices.initialized=1;
 E.pices.max_timestep_Ma=100;E.data.scalet=1;E.advection.timestep=5;
 for(int a=1;a<=8;a++){E.ien[1][1].node[a]=a;E.ien[1][2].node[a]=a+8;E.sphere.cap[1].V[1][a]=1;E.sphere.cap[1].V[1][a+8]=10;}
 if(rank==1){h[1]=4;for(int a=1;a<=8;a++)E.sphere.cap[1].V[1][a]=8;}
 if(argc>1 && !strcmp(argv[1],"nan"))E.sphere.cap[1].V[1][1]=NAN;
 if(argc>1 && !strcmp(argv[1],"length"))h[1]=0;
 double expected=size==1?1:2;
 assert(fabs(particle_cfl_rate(&E,&legacy)-expected)<1e-14);assert(legacy==10);
 pices_limit_timestep(&E);assert(E.advection.timestep<=.24/expected);
 assert(fabs(E.advection.timestep-.24/expected)<1e-7);
 cfl=checked_particle_cfl(&E,E.advection.timestep,&legacy);assert(cfl<=.24 && legacy>.24);
 E.pices.cfl=.4;E.advection.timestep=5;pices_limit_timestep(&E);
 assert(fabs(E.advection.timestep-.4/expected)<1e-7);
 if(argc>1 && !strcmp(argv[1],"oversize"))checked_particle_cfl(&E,.41/expected,&legacy);
 assert(checked_particle_cfl(&E,E.advection.timestep,&legacy)<=.4);
 /* Original candidate, Ma cap and age knot still take precedence. */
 E.advection.timestep=.01;pices_limit_timestep(&E);assert(fabs(E.advection.timestep-.01)<1e-8);
 E.pices.max_timestep_Ma=.02;E.advection.timestep=5;pices_limit_timestep(&E);assert(fabs(E.advection.timestep-.02)<1e-8);
 E.pices.max_timestep_Ma=100;E.control.lith_age=1;E.control.start_age=2.03;
 E.advection.timestep=5;pices_limit_timestep(&E);assert(fabs(E.advection.timestep-.03)<1e-7);
 /* Stagnant flow imposes no advective restriction. */
 memset(E.sphere.cap[1].V,0,sizeof(E.sphere.cap[1].V));E.control.lith_age=0;
 assert(particle_cfl_rate(&E,&legacy)==0);E.advection.timestep=5;pices_limit_timestep(&E);assert(E.advection.timestep==5);
 MPI_Finalize();return 0;
}
'''
with tempfile.TemporaryDirectory() as directory:
    root = Path(directory)
    (root/'test.c').write_text(stub+functions+main)
    subprocess.run(['mpicc','-std=gnu11',str(root/'test.c'),'-lm','-o',str(root/'test')],check=True)
    mpiexec = str(Path(shutil.which('mpicc')).resolve().with_name('mpiexec'))
    for ranks,mode in [(1,''),(2,''),(2,'nan'),(2,'length'),(2,'oversize')]:
        result=subprocess.run([mpiexec,'-n',str(ranks),str(root/'test'),mode],capture_output=True,text=True,timeout=40)
        assert (result.returncode==0)==(not mode), (mode,result.stdout,result.stderr)
    print('PASS: actual C CFL functions; serial/two-rank pairing, caps, zero flow and invalid/fixed-step guards')
