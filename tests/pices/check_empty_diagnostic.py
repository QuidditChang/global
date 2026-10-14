"""Replay empty coverage and stale flavor counters through actual C diagnostic."""
from pathlib import Path
import subprocess,tempfile,unittest
class EmptyDiagnostic(unittest.TestCase):
 def test_coverage_and_read_only_failure(self):
  source=(Path(__file__).resolve().parents[2]/'lib/Pices.c').read_text()
  body=source[source.index('static int *coverage_before;'):source.index('static void composition_check(')]
  stub=r'''
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#define MPI_LONG_LONG_INT 1
#define MPI_SUM 2
#define MPI_MAX 3
struct All_variables {
 struct {int on;} composition;
 struct {int nel,elx,ely,elz;} lmesh;
 struct {int me,world;} parallel;
 struct {int solution_cycles;double elapsed_time;} monitor;
 struct {double timestep;} advection;
 struct {double radius_km;} data;
 struct {int capid[2];double ro;} sphere;
 struct {int node[9];} ien[2][3];
 struct {double area;} eco[2][3];
 struct {FILE *fpt;int ntracers[2],ielement[2][3],nflavors,ntracer_flavor[2][2][3];} trace;
 double x[2][4][13],sx[2][4][13];
};
void empty_audit(struct All_variables *E,const int *c){}
void pices_coverage_probe(struct All_variables *E,int *c,int moved,int empty){}
static struct All_variables saved;
static double reduce(struct All_variables *E,double v,int op){return v;}
static void MPI_Allreduce(const void *a,void *b,int n,int t,int op,int world){memcpy(b,a,n*sizeof(long long));}
static void MPI_Barrier(int world){}
static void pices_fail(struct All_variables *E,const char *why){assert(!memcmp(E,&saved,sizeof(saved)));fprintf(stderr,"FAIL: %s\n",why);exit(72);}
'''.replace(' struct {int me,world;} parallel;',' FILE *fp;struct {int me,world;} parallel;')
  main=r'''
int main(int argc,char **argv){
 struct All_variables E={0};E.fp=stdout;E.trace.fpt=stderr;E.eco[1][1].area=E.eco[1][2].area=1;E.composition.on=1;E.lmesh.nel=2;E.lmesh.elx=2;E.lmesh.ely=E.lmesh.elz=1;
 E.parallel.me=0;E.monitor.solution_cycles=97;E.sphere.capid[1]=6;E.sphere.ro=1;E.data.radius_km=6371;
 E.trace.nflavors=1;E.trace.ntracers[1]=2;E.trace.ielement[1][1]=1;E.trace.ielement[1][2]=2;
 for(int e=1;e<=2;e++)for(int a=1;a<=8;a++){int n=a+4*(e-1);E.ien[1][e].node[a]=n;E.x[1][1][n]=1;E.sx[1][3][n]=.8;}
 capture_particle_coverage(&E);
 E.trace.ielement[1][1]=argc>1?1:2;E.trace.ntracer_flavor[1][0][2]=2;
 if(argc>1 && !strcmp(argv[1],"lost"))E.trace.ntracers[1]=1;
 saved=E;report_empty_elements(&E);assert(!memcmp(&E,&saved,sizeof(E)));
 report_empty_elements(&E); /* repeated observation must not increment streak */
 E.monitor.solution_cycles++;report_empty_elements(&E);
 E.trace.ielement[1][1]=1;E.trace.ntracer_flavor[1][0][1]=1;E.trace.ntracer_flavor[1][0][2]=1;
 E.monitor.solution_cycles++;report_empty_elements(&E);return 0;
}
'''
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);(p/'test.c').write_text(stub+body+main)
   subprocess.run(['cc','-std=gnu99',str(p/'test.c'),'-lm','-o',str(p/'test')],check=True)
   run=subprocess.run([str(p/'test')],capture_output=True,text=True)
   self.assertEqual(run.returncode,0,run.stderr)
   self.assertIn('empty=1 volume_nd=1 volume_fraction=0.5 longest_current=2',run.stdout)
   self.assertIn('event=recovered empty_steps=2 count=1',run.stderr)
   self.assertIn('empty=0 volume_nd=0',run.stdout)
   for mode in ('stale','lost'):
    run=subprocess.run([str(p/'test'),mode],capture_output=True,text=True)
    self.assertEqual(run.returncode,72,run.stderr)
    self.assertIn('particle conservation/host/count failure',run.stderr)
if __name__=='__main__':unittest.main()
