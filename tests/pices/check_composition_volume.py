"""Exercise the actual finite-support volume module with exact synthetic cells.

The two cells have volumes 1 and 4 and share one face. A second fixture removes
that connection to ensure insufficient support fails explicitly. MPI tests
share the fixture's nodes between ranks to exercise summed FE denominators.
No production source is patched or replaced.
"""
import csv
import json
import os
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
HARNESS = r'''
#include <assert.h>
#include <math.h>
#include <string.h>
#include "Pices_composition_volume.c"
void pices_fail(struct All_variables *E,const char *why) {
 fprintf(stderr,"rank%d %s\n",E->parallel.me,why);fflush(stderr);
 MPI_Abort(E->parallel.world,72);abort();
}
void get_global_shape_fn(struct All_variables *E,int e,struct Shape_function *n,
 struct Shape_function_dx *dx,struct Shape_function_dA *area,int pressure,
 int sphere,double rtf[4][9],int lev,int cap) {
 for(int g=1;g<=8;g++)area->vpt[g]=(e==1?1.:4.)/8;
}
void exchange_nodes(struct All_variables *E,double **p,int lev) {
 MPI_Allreduce(MPI_IN_PLACE,p[1]+1,E->lmesh.nno,MPI_DOUBLE,MPI_SUM,E->parallel.world);
}
static void legacy(struct All_variables *E) {
 for(int f=0;f<25;f++)memset(E->trace.ntracer_flavor[1][f],0,3*sizeof(int));
 for(int p=1;p<=4;p++)E->trace.ntracer_flavor[1][(int)E->trace.extraq[1][0][p]][E->trace.ielement[1][p]]++;
 for(int e=1;e<=2;e++) {
  int total=0;for(int f=0;f<25;f++)total+=E->trace.ntracer_flavor[1][f][e];
  if(total)for(int f=1;f<25;f++)E->composition.comp_el[1][f-1][e]=E->trace.ntracer_flavor[1][f][e]/(double)total;
 }
}
static void test_step(struct All_variables *E,int expected) {
 double old[24][3],oldnode[24][17],weights[5];
 for(int f=0;f<24;f++) {
  memcpy(old[f],E->composition.comp_el[1][f],3*sizeof(double));
  memcpy(oldnode[f],E->composition.comp_node[1][f],(E->lmesh.nno+1)*sizeof(double));
 }
 memcpy(weights,E->trace.extraq[1][E->pices.volume_slot],sizeof(weights));
 int success=pices_composition_volume_reconstruct(E);
 assert(success==expected);
 const struct Pices_volume_summary *s=pices_composition_volume_latest(E);
 if(success) {
  assert(s->species_error<1e-12 && s->node_error<1e-12);
  assert(s->fraction_error<1e-12 && s->partition_error<1e-12);
  for(int f=0;f<25;f++)assert(fabs(s->corrected[f]-s->target[f])<1e-11);
 }
 if(E->pices.volume_mode==1)for(int f=0;f<24;f++) {
  assert(!memcmp(old[f],E->composition.comp_el[1][f],3*sizeof(double)));
  assert(!memcmp(oldnode[f],E->composition.comp_node[1][f],(E->lmesh.nno+1)*sizeof(double)));
 }
 assert(!memcmp(weights,E->trace.extraq[1][E->pices.volume_slot],sizeof(weights)));
 pices_composition_volume_check(E);pices_composition_volume_commit(E);
}
int main(int argc,char **argv) {
 MPI_Init(&argc,&argv);struct All_variables *E=calloc(1,sizeof(*E));
 E->parallel.world=MPI_COMM_WORLD;MPI_Comm_rank(E->parallel.world,&E->parallel.me);MPI_Comm_size(E->parallel.world,&E->parallel.nproc);
 int disconnected=getenv("TEST_DISCONNECTED")!=NULL;
 E->pices.enabled=E->pices.eba=E->pices.checkpoint=1;
 E->pices.volume_mode=getenv("TEST_ACTIVE")?2:1;E->pices.volume_slot=-1;
 E->pices.volume_hops=2;E->pices.volume_iterations=100;E->pices.volume_tolerance=1e-12;
 E->trace.nflavors=25;E->trace.number_of_extra_quantities=2;
 E->sphere.caps_per_proc=1;E->composition.on=1;E->composition.ncomp=24;E->composition.ibuoy_type=1;
 E->lmesh.nel=2;E->lmesh.nno=disconnected?16:12;E->fp=stdout;
 E->exchange_node_d=exchange_nodes;
 pices_composition_volume_register(E);assert(E->pices.volume_slot==2 && E->trace.number_of_extra_quantities==3);
 E->trace.ntracers[1]=4;E->trace.ielement[1]=calloc(5,sizeof(int));
 for(int q=0;q<3;q++)E->trace.extraq[1][q]=calloc(5,sizeof(double));
 E->trace.ntracer_flavor[1]=calloc(25,sizeof(int*));
 for(int f=0;f<25;f++)E->trace.ntracer_flavor[1][f]=calloc(3,sizeof(int));
 E->composition.comp_el[1]=calloc(24,sizeof(double*));E->composition.comp_node[1]=calloc(24,sizeof(double*));
 for(int f=0;f<24;f++) {E->composition.comp_el[1][f]=calloc(3,sizeof(double));E->composition.comp_node[1][f]=calloc(E->lmesh.nno+1,sizeof(double));}
 E->ien[1]=calloc(3,sizeof(struct IEN));
 for(int d=1;d<=3;d++)E->x[1][d]=calloc(E->lmesh.nno+1,sizeof(double));
 for(int a=1;a<=8;a++) {
  E->ien[1][1].node[a]=a;E->ien[1][2].node[a]=a+(disconnected?8:4);
  for(int g=1;g<=8;g++)E->N.vpt[GNVINDEX(a,g)]=.125;
 }
 for(int n=1;n<=E->lmesh.nno;n++)E->x[1][1][n]=n<=4?0:n<=8?1:5;
 for(int p=1;p<=4;p++){E->trace.ielement[1][p]=p<=2?1:2;E->trace.extraq[1][0][p]=p<=2?24:0;}
 legacy(E);pices_composition_volume_initialize(E,0);
 assert(E->trace.extraq[1][2][1]==.5 && E->trace.extraq[1][2][4]==2);
 /* Verify the actual matrix-free kernel's symmetry and capacity invariant. */
 double x[2]={.3,2},y[2]={1.2,.7},kx[2],ky[2],one[2]={1,1},ko[2];
 kernel(E,x,kx);kernel(E,y,ky);kernel(E,one,ko);
 assert(fabs(x[0]*ky[0]+x[1]*ky[1]-y[0]*kx[0]-y[1]*kx[1])<1e-12);
 assert(fabs(ko[0]-1)<1e-12 && fabs(ko[1]-4)<1e-12);
 test_step(E,1);
 E->monitor.solution_cycles=1;E->monitor.elapsed_time=1;pices_composition_volume_begin(E);
 E->trace.ielement[1][2]=2;legacy(E);test_step(E,!disconnected);
 if(disconnected){MPI_Finalize();return 0;}
 E->monitor.solution_cycles=2;E->monitor.elapsed_time=2;pices_composition_volume_begin(E);
 E->trace.ielement[1][1]=2;legacy(E);test_step(E,1);
 assert(pices_composition_volume_latest(E)->empty_cells==E->parallel.nproc);
 E->monitor.solution_cycles=3;E->monitor.elapsed_time=3;pices_composition_volume_begin(E);
 pices_composition_volume_transition(E,1,3,0,5);E->trace.extraq[1][0][3]=5;legacy(E);test_step(E,1);
 assert(fabs(E->pices.volume_source[5]-2*E->parallel.nproc)<1e-12);
 E->monitor.solution_cycles=4;E->monitor.elapsed_time=4;pices_composition_volume_begin(E);
 pices_composition_volume_transition(E,1,3,5,0);E->trace.extraq[1][0][3]=0;
 E->trace.ielement[1][1]=E->trace.ielement[1][2]=1;legacy(E);test_step(E,1);
 pices_composition_volume_initialize(E,1);
 assert(E->pices.volume_source[5]==0 && E->pices.volume_source[24]==0);
 MPI_Finalize();return 0;
}
'''

CHAIN_MAIN = r'''
int main(int argc,char **argv) {
 MPI_Init(&argc,&argv);struct All_variables *E=calloc(1,sizeof(*E));
 E->parallel.world=MPI_COMM_WORLD;MPI_Comm_rank(E->parallel.world,&E->parallel.me);MPI_Comm_size(E->parallel.world,&E->parallel.nproc);
 E->pices.volume_hops=2;E->pices.volume_iterations=1000;E->pices.volume_tolerance=1e-12;
 E->sphere.caps_per_proc=1;E->composition.on=1;E->composition.ncomp=24;E->composition.ibuoy_type=1;
 E->lmesh.nel=4;E->lmesh.nno=20;E->fp=stdout;E->exchange_node_d=exchange_nodes;
 E->ien[1]=calloc(5,sizeof(struct IEN));
 for(int d=1;d<=3;d++)E->x[1][d]=calloc(21,sizeof(double));
 for(int e=1;e<=4;e++)for(int a=1;a<=8;a++)E->ien[1][e].node[a]=(e-1)*4+a;
 for(int a=1;a<=8;a++)for(int g=1;g<=8;g++)E->N.vpt[GNVINDEX(a,g)]=.125;
 for(int n=1;n<=20;n++)E->x[1][1][n]=(n-1)/4;
 geometry(E);struct Volume_state *v=&volume_state;struct Pices_volume_summary *s=&v->last;
 memset(s,0,sizeof(*s));s->initialized=s->feasible=1;s->shell=v->shell;
 int infeasible=getenv("TEST_CAPACITY_INFEASIBLE")!=NULL;
 if(infeasible)E->pices.volume_iterations=20;
 /* Cell 0 can reach surplus cell 1 but not cell 3 with two hops. The
  * infeasible case has capacity demand .5 > reachable surplus .3. */
 double q[4]={infeasible?.5:.9,2.3,infeasible?2.8:2.4,4.4};
 int flavors[4]={24,0,24,5};
 for(int e=0;e<4;e++)v->raw[e*25+flavors[e]]=q[e];
 for(int e=0;e<4;e++)s->target[flavors[e]]+=q[e];sums(E,s->target,25);
 int success=transport(E);
 if(infeasible) {
  assert(!success);assert(strstr(s->reason,"iteration limit"));
  assert(s->iterations==E->pices.volume_iterations && s->row_error>1e-4);
 } else {
  assert(success && s->iterations>1);assert(close_rows(E));assert(node_projection(E));
  assert(s->species_error<1e-12 && s->node_error<1e-12);
  assert(fabs(v->corrected[0*25+24]-.9)<1e-11);
  assert(fabs(v->corrected[0*25+0]-.1)<1e-11);
  assert(fabs(v->corrected[1*25+0]-2)<1e-11);
  assert(fabs(v->corrected[2*25+0]-.2)<1e-11);
  assert(fabs(v->corrected[2*25+5]-.4)<1e-11);
  assert(fabs(v->corrected[3*25+5]-4)<1e-11);
  assert(s->transfer_distance_rms_nd>0);
 }
 if(!E->parallel.me)printf("chain feasible=%d iterations=%d row_error=%.17g\n",success,s->iterations,s->row_error);
 MPI_Finalize();return 0;
}
'''


def main():
    compiler = os.environ.get('MPICC', '/usr/local/bin/mpicc')
    mpiexec = os.environ.get('MPIEXEC', '/usr/local/bin/mpiexec')
    with tempfile.TemporaryDirectory(prefix='pices-volume-module-') as temporary:
        root = Path(temporary)
        source = root / 'harness.c'
        source.write_text(HARNESS)
        exe = root / 'harness'
        subprocess.run([compiler, '-std=gnu99', '-O1', '-I' + str(ROOT / 'lib'),
                        str(source), '-lm', '-o', str(exe)], check=True)
        chain_source = root / 'chain.c'
        chain_source.write_text(HARNESS.split('int main(int argc,char **argv)')[0].replace(
            '(e==1?1.:4.)/8', 'e/8.0') + CHAIN_MAIN)
        chain_exe = root / 'chain'
        subprocess.run([compiler, '-std=gnu99', '-O1', '-I' + str(ROOT / 'lib'),
                        str(chain_source), '-lm', '-o', str(chain_exe)], check=True)
        evidence = []
        for ranks in (1, 2):
            for mode, disconnected in (('shadow', False), ('active', False),
                                       ('shadow', True), ('active', True)):
                case = root / ('%s_%d_%s' % (mode, ranks, disconnected))
                case.mkdir()
                env = dict(os.environ)
                if mode == 'active': env['TEST_ACTIVE'] = '1'
                else: env.pop('TEST_ACTIVE', None)
                if disconnected: env['TEST_DISCONNECTED'] = '1'
                else: env.pop('TEST_DISCONNECTED', None)
                run = subprocess.run([mpiexec, '--oversubscribe', '-n', str(ranks), str(exe)],
                                     cwd=case, env=env, capture_output=True, text=True, timeout=60)
                expected_failure = mode == 'active' and disconnected
                assert (run.returncode != 0) == expected_failure, (mode, ranks, disconnected, run.stderr)
                with (case / 'Diagnostic/VolumeConserve/steps.csv').open() as stream:
                    steps = list(csv.DictReader(stream))
                if disconnected:
                    failures = [r for r in steps if int(r['step']) == 1]
                    assert failures and all(r['feasible'] == '0' for r in failures)
                    assert 'unreachable' in failures[-1]['reason']
                else:
                    accepted = [r for r in steps if r['status'] == 'accepted']
                    assert [int(r['step']) for r in accepted] == list(range(5))
                    assert all(r['feasible'] == '1' for r in accepted)
                    assert all(float(r['species_error']) < 1e-12 for r in accepted)
                evidence.append(dict(ranks=ranks, mode=mode, disconnected=disconnected,
                                     expected_active_abort=expected_failure, status='PASS'))
            for infeasible in (False, True):
                env = dict(os.environ)
                if infeasible: env['TEST_CAPACITY_INFEASIBLE'] = '1'
                else: env.pop('TEST_CAPACITY_INFEASIBLE', None)
                run = subprocess.run([mpiexec, '--oversubscribe', '-n', str(ranks), str(chain_exe)],
                                     cwd=root, env=env, capture_output=True, text=True, timeout=60)
                assert run.returncode == 0, (ranks, infeasible, run.stderr)
                evidence.append(dict(ranks=ranks, constrained_chain=True,
                                     capacity_infeasible=infeasible, status='PASS', output=run.stdout.strip()))
        print(json.dumps(dict(status='PASS', cases=evidence), indent=2))


if __name__ == '__main__':
    main()
