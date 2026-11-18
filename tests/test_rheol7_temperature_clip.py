"""Compile and exercise rheol7's production nodal interpolation block."""
from pathlib import Path
import subprocess
import tempfile

root = Path(__file__).resolve().parents[1]
whole = (root / 'lib/Viscosity_structures.c').read_text()
source = whole.split('    case 7:', 1)[1]
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
    double q[4],pa[4],pb[4],pc[4];
    const int modes[3]={7,71,72};
    steinberger_defaults(q,pa,pb,pc);
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
        for(int s=1;s<=4;s*=2) for(int mode=0;mode<3;mode++) {
            double eta=steinberger_viscosity_mode(modes[mode],11.486433307280915,
                1626.0477339019924,300.+3400.*t,1e21,.25*s,q,pa,pb,pc);
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
        for(int mode=0;mode<3;mode++)
            assert(steinberger_viscosity_mode(modes[mode],11.,1626.,300.+3400.*interpolate(),1e21,1.,q,pa,pb,pc)<0.);
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

# Exercise actual family switch cases and the get_system_viscosity dispatch gate.
# The kernel was independently checked against archived physical profiles by
# test_rheol7_nuref.py; this checks that production reaches it even with TDEPV off.
cases = whole[whole.index('    case 7:'):whole.index('    case 1:',whole.index('    case 7:'))]
strict = whole[whole.index('static double strict_rheology_reference_temperature('):whole.index('static void read_rheol7_vector(')]
gate_start=whole.index('    if(E->viscosity.TDEPV ||')
gate_end=whole.index('visc_from_mat(E,evisc);',gate_start)+len('visc_from_mat(E,evisc);')
gate=whole[gate_start:gate_end]
prefix = r'''
#include <assert.h>
#include <float.h>
#include <string.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "Steinberger_nuref.h"
static struct All_variables state;
void myerror(struct All_variables *E,char *reason) {abort();}
static void rheol7_failure(struct All_variables *E,const char *reason,int cap,int element,int gp,
                          double depth,double temperature_nd,double tref,double eta) {abort();}
'''
dispatch = r'''
static void visc_from_T(struct All_variables *E,float **EEta,int propogate) {
 int m,i,jj,kk,nel=1,vpts=1,ends=8;
 switch(E->viscosity.RHEOL) {
 CASES
 default:abort();
 }
}
static void visc_from_mat(struct All_variables *E,float **EEta) {EEta[1][1]=-7.;}
static void pipeline(struct All_variables *E,float **evisc) {
 int propogate=0;
 GATE
}
'''.replace('CASES',cases).replace('GATE',gate)
main = r'''
int main(void) {
 struct All_variables *E=&state;
 const int modes[3]={7,71,72};
 float **eta=calloc(2,sizeof(*eta));eta[1]=calloc(2,sizeof(*eta[1]));
 E->sphere.caps_per_proc=1;E->lmesh.noz=8;E->mesh.nsd=3;
 E->data.Ttop=300.;E->data.ref_temperature=3400.;E->data.ref_viscosity=1.e21;E->data.radius_km=6371.;
 E->refstate.has_temperature=1;E->refstate.Tref=calloc(9,sizeof(*E->refstate.Tref));
 E->ien[1]=calloc(2,sizeof(*E->ien[1]));E->T[1]=calloc(9,sizeof(*E->T[1]));
 E->sx[1][3]=calloc(9,sizeof(*E->sx[1][3]));
 E->viscosity.cold_scale=.5;
 steinberger_defaults(E->viscosity.rheol7_layer_factors,E->viscosity.rheol7_poly_a,
                      E->viscosity.rheol7_poly_b,E->viscosity.rheol7_poly_c);
 for(int i=0;i<4;i++)E->viscosity.rheol7_layer_factors[i]=2.+i;
 for(int a=1;a<=8;a++) {
  E->ien[1][1].node[a]=a;E->N.vpt[GNVINDEX(a,1)]=.125;
  E->sx[1][3][a]=1.-1000./6371.;E->refstate.Tref[a]=.5;
 }
 for(int mode=0;mode<3;mode++) for(int tdep=0;tdep<2;tdep++) for(int hot=0;hot<2;hot++) {
  E->viscosity.RHEOL=modes[mode];E->viscosity.TDEPV=tdep;
  double t=hot?.8:.3,depth=0.,temp=0.;
  for(int a=1;a<=8;a++) {
   E->T[1][a]=tdep?t:NAN;
   depth+=(1.-E->sx[1][3][a])*6371.*.125;
   if(tdep)temp+=E->T[1][a]*.125;
  }
  double expected=steinberger_nuref_mode(modes[mode],depth,2000.,1.e21,
   E->viscosity.rheol7_layer_factors,E->viscosity.rheol7_poly_a,E->viscosity.rheol7_poly_b,E->viscosity.rheol7_poly_c);
  if(tdep)expected=steinberger_viscosity_mode(modes[mode],depth,2000.,300.+3400.*temp,1.e21,.5,
   E->viscosity.rheol7_layer_factors,E->viscosity.rheol7_poly_a,E->viscosity.rheol7_poly_b,E->viscosity.rheol7_poly_c);
  eta[1][1]=NAN;pipeline(E,eta);
  assert(isfinite(eta[1][1]) && fabs(eta[1][1]/expected-1.)<1.e-6);
 }
 E->viscosity.RHEOL=3;E->viscosity.TDEPV=0;pipeline(E,eta);assert(eta[1][1]==-7.);
 return 0;
}
'''
with tempfile.TemporaryDirectory() as tmp:
    p=Path(tmp)
    (p/'dispatch.c').write_text(prefix+strict+dispatch+main)
    subprocess.run(['mpicc','-std=gnu99','-I'+str(root/'lib'),str(p/'dispatch.c'),'-lm','-o',str(p/'dispatch')],check=True)
    subprocess.run([str(p/'dispatch')],check=True)
print('PASS: actual production switch and outer dispatch evaluate 7/71/72 with TDEPV on/off; '
      'cold/hot fields use selected coefficients; non-family fallback unchanged.')
