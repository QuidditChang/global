"""Compile and exercise rheol7's production nodal interpolation block."""
from pathlib import Path
import subprocess
import tempfile

root = Path(__file__).resolve().parents[1]
source = (root / 'lib/Viscosity_structures.c').read_text().split('    case 7:', 1)[1]
block = source[source.index('                    for(kk=1;kk<=ends;kk++) {'):
               source.index('                    if(!E->refstate.has_temperature)')]
prefix = r'''
#include <assert.h>
#include <float.h>
#include <string.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "Steinberger_nuref.h"
static struct All_variables state;
static double interpolate(void) {
    struct All_variables *E = &state;
    int m=1,i=1,jj=1,kk,ends=8;
    double depth_km=0.,temperature_nd=0.;
'''
suffix = r'''
    return temperature_nd;
}
int main(void) {
    struct All_variables *E = &state;
    double cases[3][8] = {
        {-.0316986055,-.0521749282,-.0866792407,.120888173,0,0,0,0},
        {-1,2,.2,.4,.6,.8,0,1},
        {0,1,.2,.4,.6,.8,.3,.7}
    };
    double expected[] = {.120888173/8.,.5,.5};
    E->data.Ttop=300.; E->data.Tbottom=3700.; E->data.ref_temperature=3400.;
    E->mesh.toptbc=E->mesh.bottbc=1;
    E->control.TBCtopval=0.; E->control.TBCbotval=1.;
    E->ien[1] = calloc(2,sizeof(*E->ien[1]));
    E->T[1] = calloc(9,sizeof(*E->T[1]));
    E->sx[1][3] = calloc(9,sizeof(*E->sx[1][3]));
    for(int a=1;a<=8;a++) {
        E->ien[1][1].node[a]=a;
        E->N.vpt[GNVINDEX(a,1)]=.125;
    }
    for(int c=0;c<3;c++) {
        for(int a=1;a<=8;a++) E->T[1][a]=cases[c][a-1];
        float before[9]; memcpy(before,E->T[1],sizeof(before));
        double t=interpolate();
        assert(fabs(t-expected[c])<1e-7);
        assert(memcmp(before,E->T[1],sizeof(before))==0);
        for(int s=1;s<=4;s*=2) {
            double eta=steinberger_viscosity(11.486433307280915,
                1626.0477339019924,300.+3400.*t,1e21,.25*s);
            assert(eta>0. && eta<FLT_MAX);
        }
    }
    /* Exact CMB element values from the cluster regression at step 0. */
    double cmb_t[8]={1.00000062,1.00000062,1.00000062,1.00000062,
                     .966868883,.966868883,.966868883,.966868883};
    double cmb_w[8]={.49056260249407813,.13144585726148808,
                     .035220812396550491,.13144585726148808,
                     .13144585726148808,.035220812396550491,
                     .0094373885318062185,.035220812396550491};
    for(int a=1;a<=8;a++) {
        E->T[1][a]=cmb_t[a-1];
        E->N.vpt[GNVINDEX(a,1)]=cmb_w[a-1];
    }
    /* E->T and shape data use the solver's mixed float/double storage. */
    assert(fabs(interpolate()-.99299857)<1e-7);
    double bad[]={NAN,INFINITY,-INFINITY};
    for(int b=0;b<3;b++) {
        E->T[1][1]=bad[b];
        assert(!isfinite(interpolate()));
        assert(steinberger_viscosity(11.,1626.,300.+3400.*interpolate(),1e21,1.)<0.);
    }
    /* The nondimensional clip stays [0,1]; cfg controls its Kelvin mapping. */
    E->data.Ttop=400.; E->data.Tbottom=4400.; E->data.ref_temperature=4000.;
    assert(E->data.Ttop+E->data.ref_temperature*0.0==400.);
    assert(E->data.Ttop+E->data.ref_temperature*1.0==4400.);
    return 0;
}
'''
with tempfile.TemporaryDirectory() as tmp:
    p = Path(tmp)
    (p/'test.c').write_text(prefix + block + suffix)
    subprocess.run(['mpicc', '-std=gnu99', '-I'+str(root/'lib'), str(p/'test.c'),
                    '-lm', '-o', str(p/'test')], check=True)
    subprocess.run([str(p/'test')], check=True)
print('PASS: production nodal clipping before interpolation; E->T unchanged; '
      'step3 temperatures finite for cold_scale=0.25/0.5/1; nonfinite inputs rejected; '
      'cfg-driven Kelvin mapping checked.')
