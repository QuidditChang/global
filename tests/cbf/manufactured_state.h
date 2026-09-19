/* Analytic sampled steady conduction on the actual CitcomS mesh.
 * This tests extraction/geometry/MPI; it is not a discrete PDE solve. */
static void cbf_manufactured_state(struct All_variables *E)
{
    int m,n,d;
    double B=E->sphere.ri*E->sphere.ro/(E->sphere.ro-E->sphere.ri);
    E->control.disptn_number=0;E->control.Q0=0;
    E->control.kd_upper_prefactor=E->data.ks;
    E->control.kd_lower_prefactor=E->data.ks;
    E->control.kd_upper_linear=E->control.kd_upper_quadratic=0;
    E->control.kd_lower_linear=E->control.kd_lower_quadratic=0;
    E->control.kT_exponent=0;E->control.kC_ratio=1;
    for(d=0;d<PHASE_TRANSITIONS;++d)E->control.phase[d].entropy_jump=0;
    for(m=1;m<=E->sphere.caps_per_proc;++m)
        for(n=1;n<=E->lmesh.nno;++n) {
            E->T[m][n]=B*(1/E->sx[m][3][n]-1/E->sphere.ro);
            E->Tdot[m][n]=0;
            for(d=1;d<=3;++d)E->sphere.cap[m].V[d][n]=0;
        }
}
