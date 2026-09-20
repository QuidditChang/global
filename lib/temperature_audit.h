#ifndef CITCOM_TEMPERATURE_AUDIT_H
#define CITCOM_TEMPERATURE_AUDIT_H
#include <float.h>
#include <math.h>

/* Read-only, rank-local records. Shared nodes may appear on several ranks;
 * counts are NOT a global volume fraction. Log every invalid node so that
 * checkpoint, predictor, corrector and assimilation can be compared by ID. */
static void audit_temperature(struct All_variables *E, const char *stage,
                              int attempt, int pass)
{
    int m,i,negative=0,nonfinite=0,mincap=0,minnode=0;
    double lo=DBL_MAX,hi=-DBL_MAX,t;
    if(!E->control.temperature_audit) return;
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(i=1;i<=E->lmesh.nno;i++) {
            t=E->data.Ttop+E->data.ref_temperature*E->T[m][i];
            if(isfinite(t)) {
                if(t<lo) {lo=t;mincap=m;minnode=i;}
                if(t>hi) hi=t;
                if(t<0.0) negative++;
            } else nonfinite++;
            if(!isfinite(t) || t<0.0)
                fprintf(E->fp,"TEMP_AUDIT_NODE stage=%s step=%d attempt=%d pass=%d rank=%d cap=%d node=%d theta=%.17g phi=%.17g radius=%.17g T_K=%.17g T_nd=%.17g Tdot=%.17g flags=%u\n",
                    stage,E->monitor.solution_cycles,attempt,pass,E->parallel.me,
                    E->sphere.capid[m],i,E->sx[m][1][i],E->sx[m][2][i],
                    E->sx[m][3][i],t,E->T[m][i],E->Tdot[m][i],E->node[m][i]);
        }
    fprintf(E->fp,"TEMP_AUDIT stage=%s step=%d attempt=%d pass=%d rank=%d time_nd=%.17g dt_nd=%.17g min_K=%.17g max_K=%.17g min_cap_local=%d min_node=%d negative=%d nonfinite=%d\n",
        stage,E->monitor.solution_cycles,attempt,pass,E->parallel.me,
        E->monitor.elapsed_time,E->advection.timestep,lo,hi,mincap,minnode,
        negative,nonfinite);
    fflush(E->fp);
}
#endif
