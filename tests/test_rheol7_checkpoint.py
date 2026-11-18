"""Check actual PICES physics fingerprints distinguish the new rheology family.

The production fingerprint and coupled-mode selector are compiled unchanged.
Forcing-file hashing is stubbed because this test varies only the rheology;
the full MPI smoke test exercises actual forcing/checkpoint I/O.
"""
from pathlib import Path
import subprocess
import tempfile

root=Path(__file__).resolve().parents[1]
source=(root/'lib/Pices_checkpoint.c').read_text()
coupled=source[source.index('int pices_checkpoint_coupled('):source.index('static int schema(')]
fingerprint=source[source.index('static void fingerprint('):source.index('static void metadata_policy(')]
tag_start=fingerprint.index(' /* The new baseline 7')
tag_end=fingerprint.index(' if(strcmp(E->control.pices_stokes_inner_policy,',tag_start)
untagged=(fingerprint[:tag_start]+fingerprint[tag_end:]).replace('static void fingerprint(', 'static void untagged_fingerprint(')
prefix=r'''
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "global_defs.h"
#include "pices_sha256.h"
#include "Steinberger_nuref.h"
static struct All_variables state;
static void forcing_hash(struct All_variables *E,PicesSHA *s) {
 const char forcing[]="fixed-forcing-for-rheology-isolation";
 pices_sha_add(s,forcing,sizeof(forcing));
}
'''
main=r'''
static void defaults(struct All_variables *E) {
 steinberger_defaults(E->viscosity.rheol7_layer_factors,E->viscosity.rheol7_poly_a,
                      E->viscosity.rheol7_poly_b,E->viscosity.rheol7_poly_c);
}
int main(void) {
 struct All_variables *E=&state;
 char original[65],changed[65],untagged[65];
 E->advection.fixed_timestep=1.;E->pices.eba=1;E->control.RH=1;
 strcpy(E->control.pices_stokes_inner_policy,"current_rhs");
 E->control.pices_stokes_inner_rtol=1.e-6;
 E->lmesh.nno=E->lmesh.noz=1;
 for(int d=1;d<=3;d++) {
  E->x[1][d]=calloc(2,sizeof(*E->x[1][d]));
  E->sphere.cap[1].TB[d]=calloc(2,sizeof(*E->sphere.cap[1].TB[d]));
 }
 E->refstate.rho=calloc(2,sizeof(double));
 E->refstate.thermal_expansivity=calloc(2,sizeof(double));
 E->refstate.gravity=calloc(2,sizeof(double));
 E->refstate.heat_capacity=calloc(2,sizeof(double));
 E->mat[1]=calloc(1,sizeof(*E->mat[1]));
 defaults(E);
 E->viscosity.RHEOL=3;assert(!pices_checkpoint_coupled(E));
 E->viscosity.RHEOL=7;assert(pices_checkpoint_coupled(E));
 fingerprint(E,original,0);untagged_fingerprint(E,untagged,0);
 assert(strcmp(original,untagged)); /* historical7 cannot silently become baseline7 */
 E->viscosity.rheol7_layer_factors[0]=9.;E->viscosity.rheol7_poly_a[0]=9.;
 fingerprint(E,changed,0);assert(!strcmp(original,changed)); /* inactive coefficients */
 E->viscosity.RHEOL=71;assert(pices_checkpoint_coupled(E));defaults(E);
 fingerprint(E,original,0);
 E->viscosity.rheol7_layer_factors[2]=1.000000000000001;
 fingerprint(E,changed,0);assert(strcmp(original,changed)); /* preserve double precision */
 defaults(E);E->viscosity.rheol7_poly_c[0]=9.;
 fingerprint(E,changed,0);assert(!strcmp(original,changed));
 E->viscosity.RHEOL=72;assert(pices_checkpoint_coupled(E));defaults(E);
 fingerprint(E,original,0);
 double *arrays[]={E->viscosity.rheol7_poly_a,E->viscosity.rheol7_poly_b,E->viscosity.rheol7_poly_c};
 for(int j=0;j<3;j++) {
  defaults(E);arrays[j][1]+=.001;
  fingerprint(E,changed,0);assert(strcmp(original,changed));
 }
 defaults(E);E->viscosity.rheol7_layer_factors[0]=9.;
 fingerprint(E,changed,0);assert(!strcmp(original,changed));
 E->viscosity.cold_scale=.5;fingerprint(E,changed,0);assert(strcmp(original,changed));
 return 0;
}
'''
with tempfile.TemporaryDirectory() as temporary:
    path=Path(temporary)
    (path/'fingerprint.c').write_text(prefix+coupled+fingerprint+untagged+main)
    subprocess.run(['mpicc','-std=gnu99','-I'+str(root/'lib'),str(path/'fingerprint.c'),'-lm','-o',str(path/'fingerprint')],check=True)
    subprocess.run([str(path/'fingerprint')],check=True)
print('PASS: all-family PICES coupled selection; baseline version tag; active double coefficients '
      'and cold_scale alter fingerprints; inactive coefficients do not.')
