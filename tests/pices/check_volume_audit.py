"""Exercise actual audit hooks with unequal cells, labels, empties and MPI.

The geometry stub supplies exact quadrature for two volumes 1 and 4. Both
cells share their eight nodes, deliberately exercising owned-element nodal
weights. Each MPI rank owns a copy of this two-cell patch. No transport
algorithm is mocked into the audit: particle positions/labels are explicit.
"""
from pathlib import Path
import csv
import json
import os
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
HARNESS = r'''
#include <assert.h>
#include <math.h>
#include <string.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "pices_volume_audit.h"
void pices_fail(struct All_variables *E,const char *reason) {
 fprintf(stderr,"rank %d: %s\n",E->parallel.me,reason);MPI_Abort(E->parallel.world,72);abort();
}
void get_global_shape_fn(struct All_variables *E,int e,struct Shape_function *n,
 struct Shape_function_dx *dx,struct Shape_function_dA *area,int pressure,
 int sphere,double rtf[4][9],int lev,int cap) {
 for(int g=1;g<=8;g++)area->vpt[g]=(e==1?1.0:4.0)/8;
}
static void recount(struct All_variables *E) {
 for(int f=0;f<25;f++)memset(E->trace.ntracer_flavor[1][f],0,3*sizeof(int));
 for(int p=1;p<=4;p++)E->trace.ntracer_flavor[1][(int)E->trace.extraq[1][0][p]][E->trace.ielement[1][p]]++;
}
static void reconstruct(struct All_variables *E) {
 recount(E);
 for(int e=1;e<=2;e++) {
  int total=0;for(int f=0;f<25;f++)total+=E->trace.ntracer_flavor[1][f][e];
  if(total)for(int f=1;f<25;f++)E->composition.comp_el[1][f-1][e]=E->trace.ntracer_flavor[1][f][e]/(double)total;
 }
 for(int f=1;f<25;f++)for(int n=1;n<=8;n++)E->composition.comp_node[1][f-1][n]=
  (E->composition.comp_el[1][f-1][1]+4*E->composition.comp_el[1][f-1][2])/5;
 /* Isolated precision fixture: raw overshoots remain untouched by audit. */
 if(getenv("AUDIT_NODE_BOUND_DELTA"))for(int n=1;n<=8;n++) {
  for(int f=0;f<24;f++)E->composition.comp_node[1][f][n]=0;
  E->composition.comp_node[1][23][n]=1+atof(getenv("AUDIT_NODE_BOUND_DELTA"));
 }
}
/* Audit calls must not change the actual particle/element/node arrays. */
static void readonly_call(struct All_variables *E,void (*fn)(struct All_variables *)) {
 double old_comp[24][3],old_node[24][9],old_label[5];int old_host[5],old_count[25][3];
 for(int f=0;f<24;f++) {
  memcpy(old_comp[f],E->composition.comp_el[1][f],sizeof old_comp[f]);
  memcpy(old_node[f],E->composition.comp_node[1][f],sizeof old_node[f]);
 }
 for(int f=0;f<25;f++)memcpy(old_count[f],E->trace.ntracer_flavor[1][f],sizeof old_count[f]);
 memcpy(old_label,E->trace.extraq[1][0],sizeof old_label);memcpy(old_host,E->trace.ielement[1],sizeof old_host);
 fn(E);
 for(int f=0;f<24;f++) {
  assert(!memcmp(old_comp[f],E->composition.comp_el[1][f],sizeof old_comp[f]));
  assert(!memcmp(old_node[f],E->composition.comp_node[1][f],sizeof old_node[f]));
 }
 for(int f=0;f<25;f++)assert(!memcmp(old_count[f],E->trace.ntracer_flavor[1][f],sizeof old_count[f]));
 assert(!memcmp(old_label,E->trace.extraq[1][0],sizeof old_label));assert(!memcmp(old_host,E->trace.ielement[1],sizeof old_host));
}
int main(int argc,char **argv) {
 struct All_variables *E=calloc(1,sizeof(*E));MPI_Init(&argc,&argv);
 E->parallel.world=MPI_COMM_WORLD;MPI_Comm_rank(MPI_COMM_WORLD,&E->parallel.me);MPI_Comm_size(MPI_COMM_WORLD,&E->parallel.nproc);
 E->fp=stdout;E->pices.enabled=E->pices.eba=E->composition.on=1;E->composition.ibuoy_type=1;
 E->trace.nflavors=25;E->composition.ncomp=24;E->sphere.caps_per_proc=1;
 E->lmesh.nel=E->lmesh.elz=E->mesh.elz=2;E->lmesh.nno=8;E->lmesh.nzs=1;
 E->control.start_age=100;E->data.scalet=1;E->data.radius_km=6371;E->sphere.ro=1;
 E->trace.ntracers[1]=4;E->trace.ielement[1]=calloc(5,sizeof(int));E->trace.extraq[1][0]=calloc(5,sizeof(double));
 E->trace.ntracer_flavor[1]=calloc(25,sizeof(int*));
 for(int f=0;f<25;f++)E->trace.ntracer_flavor[1][f]=calloc(3,sizeof(int));
 E->composition.comp_el[1]=calloc(24,sizeof(double*));E->composition.comp_node[1]=calloc(24,sizeof(double*));
 for(int f=0;f<24;f++) { E->composition.comp_el[1][f]=calloc(3,sizeof(double));E->composition.comp_node[1][f]=calloc(9,sizeof(double)); }
 E->ien[1]=calloc(3,sizeof(struct IEN));E->sx[1][3]=calloc(9,sizeof(double));
 for(int n=1;n<=8;n++) {
  E->ien[1][1].node[n]=E->ien[1][2].node[n]=n;E->sx[1][3][n]=.8;
  for(int g=1;g<=8;g++)E->N.vpt[GNVINDEX(n,g)]=.125;
 }
 for(int p=1;p<=4;p++){E->trace.ielement[1][p]=p<=2?1:2;E->trace.extraq[1][0][p]=p<=2?24:0;}
 reconstruct(E);readonly_call(E,pices_volume_audit_initialize);readonly_call(E,pices_volume_audit_commit);
 for(int step=1;step<=4;step++) {
  E->monitor.solution_cycles=step;readonly_call(E,pices_volume_audit_begin);
  if(step==1)E->trace.ielement[1][2]=2;
  if(step==2)E->trace.ielement[1][1]=2;
  if(step==4)E->trace.ielement[1][3]=1;
  readonly_call(E,pices_volume_audit_moved);
  if(step==1) {pices_volume_audit_transition(E,1,2,0,5);E->trace.extraq[1][0][3]=5;}
  reconstruct(E);E->monitor.elapsed_time=step;readonly_call(E,pices_volume_audit_rebuilt);
  if(step<4)readonly_call(E,pices_volume_audit_commit);
 }
 MPI_Finalize();return 0;
}
'''

def rows(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def check(path, ranks):
    metadata=json.loads((path/'metadata.json').read_text())
    assert metadata['baseline_step']==0 and metadata['radial_layers']==2
    data=rows(path/'globals.csv')
    accepted={(int(r['step']),int(r['flavor'])):r for r in data if r['status']=='accepted'}
    assert len(accepted)==4*25 and not any(k[0]==4 for k in accepted)
    expected=[1,7/3,3,3]
    for step,c24 in enumerate(expected):
        r=accepted[step,24]
        assert abs(float(r['volume_after_nd'])/ranks-c24)<1e-12
        assert int(r['count_before'])==int(r['count_moved'])==int(r['count_after'])==2*ranks
        assert float(r['source_nd'])==0
        for flavor in range(25):
            rr=accepted[step,flavor]
            assert abs(float(rr['residual_nd']))<1e-12*ranks
            assert abs(float(rr['node_volume_nd'])-float(rr['volume_after_nd']))<1e-12*ranks
        assert abs(sum(float(accepted[step,f]['volume_after_nd']) for f in range(25))-5*ranks)<1e-12*ranks
    changed=[r for r in rows(path/'transitions.csv') if r['status']=='accepted']
    assert len(changed)==1 and changed[0]['old_flavor']=='0' and changed[0]['new_flavor']=='5'
    assert int(changed[0]['count'])==ranks and abs(float(changed[0]['volume_nd'])-4*ranks/3)<1e-12*ranks
    bins={(int(r['step']),r['bin']):r for r in rows(path/'bins.csv') if r['kind']=='occupancy' and r['status']=='accepted'}
    assert float(bins[2,'new_empty']['elements'])==ranks
    assert float(bins[3,'persistent_empty']['elements'])==ranks
    assert abs(float(bins[2,'new_empty']['C24_after_nd'])-ranks)<1e-12
    pending={(int(r['step']),r['bin']):r for r in rows(path/'bins.csv') if r['kind']=='occupancy' and r['status']=='pending'}
    assert float(pending[4,'recovered']['elements'])==ranks
    for step in range(4):
        rad=[r for r in rows(path/'radial.csv') if int(r['step'])==step and r['status']=='accepted']
        assert len(rad)==2
        assert sum(float(r['particles']) for r in rad)==4*ranks
        assert abs(sum(float(r['C24_after_nd']) for r in rad)-expected[step]*ranks)<1e-12*ranks
    return accepted


def main():
    compiler=shutil.which('mpicc');launcher=shutil.which('mpiexec')
    if not compiler or not launcher:
        raise SystemExit('MPI compiler/runtime required for volume audit unit harness')
    with tempfile.TemporaryDirectory(prefix='citcom-volume-audit-') as directory:
        tmp=Path(directory);(tmp/'test.c').write_text(HARNESS)
        subprocess.run([compiler,'-std=gnu99','-w','-I'+str(ROOT/'lib'),str(tmp/'test.c'),str(ROOT/'lib/Pices_volume_audit.c'),'-lm','-o',str(tmp/'test')],check=True)
        results=[]
        for n in (1,2):
            work=tmp/str(n);work.mkdir()
            env=dict(os.environ,PICES_VOLUME_AUDIT='1',PICES_VOLUME_AUDIT_INTERVAL='1',OMPI_ALLOW_RUN_AS_ROOT='1',OMPI_ALLOW_RUN_AS_ROOT_CONFIRM='1')
            run=subprocess.run([launcher,'-n',str(n),str(tmp/'test')],cwd=work,env=env,text=True,capture_output=True,timeout=60)
            assert run.returncode==0,(run.stdout,run.stderr)
            assert 'VOLUME_AUDIT_ANOMALY' not in run.stderr,run.stderr
            results.append(check(work/'Diagnostic/VolumeAudit',n))
        # Float node projection roundoff is distinct from double ledger closure.
        for name, delta, anomalous in [('node_roundoff', 5e-8, False), ('node_violation', 5e-6, True)]:
            work=tmp/name;work.mkdir()
            env=dict(os.environ,PICES_VOLUME_AUDIT='1',PICES_VOLUME_AUDIT_INTERVAL='1',AUDIT_NODE_BOUND_DELTA=str(delta))
            run=subprocess.run([str(tmp/'test')],cwd=work,env=env,text=True,capture_output=True,timeout=60)
            assert run.returncode==0,run.stderr
            assert ('VOLUME_AUDIT_ANOMALY' in run.stderr)==anomalous,run.stderr
            first=next(r for r in rows(work/'Diagnostic/VolumeAudit/globals.csv') if r['status']=='accepted' and r['flavor']=='24')
            assert abs(float(first['node_volume_nd'])/float(first['shell_volume_nd'])-(1+delta))<1e-14
            meta=json.loads((work/'Diagnostic/VolumeAudit/metadata.json').read_text())
            assert meta['element_fraction_tolerance']==1e-10 and meta['node_fraction_tolerance']==1e-7
        # Off switch must create no diagnostic directory, even with audit hooks.
        work=tmp/'off';work.mkdir()
        run=subprocess.run([str(tmp/'test')],cwd=work,env=dict(os.environ,PICES_VOLUME_AUDIT='0'),text=True,capture_output=True,timeout=60)
        assert run.returncode==0 and not (work/'Diagnostic').exists(),run.stderr
        for key in results[0]:
            for field in ('volume_before_nd','volume_moved_nd','volume_after_nd','node_volume_nd','source_nd','residual_nd'):
                assert abs(float(results[1][key][field])/2-float(results[0][key][field]))<1e-12
    print('PASS: unequal-cell drift, immutable C24, conversion budget, empty/recovered cells, shared-node integration, one/two MPI ranks, pending failure exclusion, float-node/double-ledger precision and observational/off behavior')

if __name__=='__main__':
    main()
