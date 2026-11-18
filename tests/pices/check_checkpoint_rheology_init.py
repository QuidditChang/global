"""Exercise production initialization/validation with dirty reserved flags.

The allocator supplies nonzero FREEZE/weak-region bytes, reproducing reused
malloc storage. Both input frontends use citcom_init and neither reads these
reserved switches. MPI rank queries and unrelated phase helpers are stubbed;
the actual constructor, coupled selector and validation function are compiled.
"""
import os
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
source = (ROOT / "lib/Pices.c").read_text()
validate = source[source.index("void pices_validate("):
                  source.index("/* Split the existing PG latent term")]
source = (ROOT / "lib/Pices_checkpoint.c").read_text()
coupled = source[source.index("int pices_checkpoint_coupled("):
                 source.index("static int schema(")]

HARNESS = r'''
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "global_defs.h"
#include "pices.h"
#include "phase_change.h"
#include "Steinberger_nuref.h"

static int rank_query(MPI_Comm world, int *rank) {*rank=16;return 0;}
static int size_query(MPI_Comm world, int *size) {*size=384;return 0;}
int get_process_identifier(void) {return 123;}
static void *dirty_allocation(size_t size) {
    struct All_variables *E=calloc(1,size);
    assert(E && size==sizeof(*E));
    E->viscosity.FREEZE=0x5a5a5a5a;
    E->viscosity.weak_blobs=0x5a5a5a5a;
    E->viscosity.weak_zones=0x5a5a5a5a;
    return E;
}
#define malloc dirty_allocation
#define MPI_Comm_rank rank_query
#define MPI_Comm_size size_query
#include "Citcom_init.c"
#undef MPI_Comm_size
#undef MPI_Comm_rank
#undef malloc

void pices_fail(struct All_variables *E,const char *why) {
    fprintf(stderr,"PICES_ERROR rank=%d step=%d: %s\n",
            E->parallel.me,E->monitor.solution_cycles,why);
    exit(72);
}
int phase_change_count(const struct All_variables *E) {return 3;}
double phase_change_entropy(const struct All_variables *E,int index) {return 0;}
COUPLED
VALIDATE

int main(int argc,char **argv) {
    MPI_Comm world=MPI_COMM_WORLD;
    struct All_variables *E=citcom_init(&world);
    const char *reject=argv[3];
    E->pices.enabled=E->pices.eba=E->pices.p4=E->pices.checkpoint=1;
    E->pices.cfl=.5;E->pices.max_timestep_Ma=.1;E->pices.max_substeps=10000;
    E->viscosity.RHEOL=atoi(argv[1]);E->viscosity.TDEPV=1;
    E->viscosity.FROM_SYSTEM=E->viscosity.update_allowed=1;
    E->viscosity.num_mat=4;
    E->control.pPv_enabled=atoi(argv[2]);
    E->control.pPv_k_ratio=E->control.pPv_eta_ratio=1;
    E->control.tracer=E->control.eba_formulation=1;
    E->control.lith_age=E->control.lith_age_asml=E->control.lith_age_time=1;
    E->control.kC_primordial_flavor=24;
    E->control.vbcs_file=1;E->data.timedir=1;
    E->composition.on=E->composition.ibuoy_type=1;E->composition.ncomp=24;
    E->trace.nflavors=25;E->trace.reclassify_flavors=1;E->trace.itperel=1;
    E->sphere.caps=12;E->sphere.caps_per_proc=1;
    E->mesh.topvbc=E->mesh.toptbc=E->mesh.bottbc=1;
    E->advection.ADVECTION=1;E->advection.fixed_timestep=0;
    strcpy(E->output.format,"ascii-gz");
    if(!strcmp(reject,"update_allowed"))E->viscosity.update_allowed=0;
    if(!strcmp(reject,"FROM_SYSTEM"))E->viscosity.FROM_SYSTEM=0;
    if(!strcmp(reject,"RHEOL"))E->viscosity.RHEOL=3;
#define REJECT_FLAG(flag) if(!strcmp(reject,#flag))E->viscosity.flag=1
    REJECT_FLAG(SDEPV);REJECT_FLAG(PDEPV);REJECT_FLAG(CDEPV);
    REJECT_FLAG(FREEZE);REJECT_FLAG(channel);REJECT_FLAG(wedge);
    REJECT_FLAG(weak_blobs);REJECT_FLAG(weak_zones);
    pices_validate(E);
    free(E);
    return 0;
}
'''

with tempfile.TemporaryDirectory(prefix="checkpoint-rheology-init-") as tmp:
    path = Path(tmp)
    c = path / "check.c"
    c.write_text(HARNESS.replace("COUPLED\n", coupled)
                 .replace("VALIDATE\n", validate))
    exe = path / "check"
    subprocess.run([os.environ.get("MPICC", "mpicc"), "-std=gnu99",
                    "-I" + str(ROOT / "lib"), str(c), "-lm", "-o", str(exe)],
                   check=True)
    for mode in (1, 7, 71, 72):
        for pPv in (0, 1):
            run = subprocess.run([str(exe), str(mode), str(pPv), "none"],
                                 capture_output=True, text=True)
            assert run.returncode == 0, (mode, pPv, run.stderr)
    for flag in ("update_allowed", "FROM_SYSTEM", "RHEOL", "SDEPV", "PDEPV",
                 "CDEPV", "FREEZE", "channel", "wedge", "weak_blobs", "weak_zones"):
        run = subprocess.run([str(exe), "7", "0", flag],
                             capture_output=True, text=True)
        assert run.returncode == 72, (flag, run.returncode, run.stderr)
        assert "checkpoint requires supported Newtonian rheology" in run.stderr
        assert flag + "=" in run.stderr, (flag, run.stderr)
print("PASS: dirty-allocation initialization for rheol 1/7/71/72 with pPv off/on; "
      "all 11 unsupported rheology conditions still rejected with flag diagnostics.")
