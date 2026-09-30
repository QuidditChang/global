/* Test driver only: prescribed temperatures/velocity; production PICES unchanged. */
#include <math.h>
#include <string.h>
#include "pices.h"
static void p1_test_velocity(struct All_variables *E)
{
    int n;
    const char *mode=getenv("PICES_TEST_MODE");
    double omega=mode && !strncmp(mode,"rotation",8) ? 1.0 : 0.0;
    for(n=1;n<=E->lmesh.nno;n++) {
        E->sphere.cap[1].V[1][n]=E->sphere.cap[1].V[3][n]=0;
        E->sphere.cap[1].V[2][n]=omega*E->sx[1][3][n]*sin(E->sx[1][1][n]);
        if(mode && (!strcmp(mode,"rotation_x") || !strcmp(mode,"rotation_y"))) {
            double th=E->sx[1][1][n],ph=E->sx[1][2][n],r=E->sx[1][3][n];
            int xaxis=!strcmp(mode,"rotation_x");
            E->sphere.cap[1].V[1][n]=r*(xaxis ? -sin(ph):cos(ph));
            E->sphere.cap[1].V[2][n]=-r*cos(th)*(xaxis ? cos(ph):sin(ph));
        }
    }
}
static void p1_test_initial(struct All_variables *E)
{
    const char *mode=getenv("PICES_TEST_MODE");
    int n,p,a,nodes[9];double w[9];
    FILE *f;char path[512];
    double lo=E->sphere.ri,hi=E->sphere.ro;
    for(n=1;n<=E->lmesh.nno;n++) {
        double r=E->sx[1][3][n];
        E->T[1][n]=.5;
        if(mode && !strcmp(mode,"diffusion"))
            E->T[1][n]+=.1*sin(M_PI*(r-lo)/(hi-lo))/r;
        if(mode && !strcmp(mode,"rotation"))
            E->T[1][n]+=.1*sin(M_PI*(r-lo)/(hi-lo))*sin(E->sx[1][1][n])*cos(E->sx[1][2][n]);
        if(mode && !strcmp(mode,"nonlinear_steady")) {
            double offset=E->data.Ttop/E->data.ref_temperature;
            double exponent=1-E->control.kT_exponent;
            double f=(1/lo-1/r)/(1/lo-1/hi);
            E->T[1][n]=pow((1-f)*pow(1+offset,exponent)+f*pow(offset,exponent),1/exponent)-offset;
        }
        for(a=1;a<=3;a++) E->sphere.cap[1].TB[a][n]=
            mode && !strcmp(mode,"nonlinear_steady") ? E->T[1][n] : .5;
    }
    for(p=1;p<=E->trace.ntracers[1];p++) {
        double t=0;
        tracer_temperature_weights(E,1,p,nodes,w);
        for(a=1;a<=8;a++)t+=w[a]*E->T[1][nodes[a]];
        E->trace.extraq[1][E->pices.slot][p]=getenv("PICES_TEST_IDS") ? .25+(E->parallel.me*E->trace.ntracers[1]+p)*1e-6 : t;
    }
    if(mode && !strcmp(mode,"empty")) { E->trace.ntracers[1]=0; E->trace.ilast_tracer_count=0; }
    /* Local FE matrices for independent symmetry/PSD and analytic tests. */
    snprintf(path,sizeof(path),"matrix.%d.txt",E->parallel.me);f=fopen(path,"w");
    for(n=1;n<=E->lmesh.nel;n++) {
        int b;
        for(a=0;a<8;a++)for(b=0;b<8;b++)fprintf(f,"%.17g%c",E->pices.K[n*64+a*8+b],b==7?'\n':' ');
    }
    fclose(f);
    p1_test_velocity(E);
}

static void p1_test_snapshot(struct All_variables *E)
{
    FILE *f; char path[128]; int n,p,a,nodes[9];double w[9];
    if(getenv("PICES_TEST_IDS") && E->monitor.solution_cycles!=0 && E->monitor.solution_cycles!=E->advection.max_timesteps && !E->control.restart) return;
    snprintf(path,sizeof(path),"nodes.%d.%d.txt",E->parallel.me,E->monitor.solution_cycles);
    f=fopen(path,"w");
    for(n=1;n<=E->lmesh.nno;n++)
        fprintf(f,"%d %.17g %.17g %.17g %.17g %d %.17g\n",n,E->x[1][1][n],E->x[1][2][n],E->x[1][3][n],(double)E->T[1][n],!!(E->node[1][n]&(TBX|TBY|TBZ)),E->pices.mass[n]);
    fclose(f);
    snprintf(path,sizeof(path),"particles.%d.%d.txt",E->parallel.me,E->monitor.solution_cycles);
    f=fopen(path,"w");
    for(p=1;p<=E->trace.ntracers[1];p++) {
        tracer_temperature_weights(E,1,p,nodes,w);
        fprintf(f,"%.17g",E->trace.extraq[1][E->pices.slot][p]);
        for(a=1;a<=8;a++) fprintf(f," %d %.17g",nodes[a],w[a]);
        fprintf(f," %.17g %.17g %.17g\n",E->trace.basicq[1][3][p],E->trace.basicq[1][4][p],E->trace.basicq[1][5][p]);
    }
    fclose(f);
    snprintf(path,sizeof(path),"elements.%d.txt",E->parallel.me);f=fopen(path,"w");
    for(n=1;n<=E->lmesh.nel;n++) {
        for(a=1;a<=8;a++)fprintf(f,"%d ",E->ien[1][n].node[a]);
        fprintf(f,"%.17g\n",E->pices.length[n]);
    }
    fclose(f);
}

static void p1_test_restart_bcs(struct All_variables *E)
{
    int a,n;
    for(a=1;a<=3;a++)for(n=1;n<=E->lmesh.nno;n++)E->sphere.cap[1].TB[a][n]=.5;
}
