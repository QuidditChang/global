/* Test-only suffix compiled with a copy of Pices.c. Never linked into production.
 * Evaluate real P/Q and FE energy before any time advancement. */
void pices_projection_probe(struct All_variables *E)
{
    int n,p,a,e,nodes[9],nn=E->lmesh.nno;
    double w[9],row[17],*raw=array(E,nn+1),*bc=array(E,nn+1),*local_mass=array(E,nn+1);
    double energies[3],error=0;char path[128];FILE *f;
    const char *profile=getenv("PICES_PROBE_PROFILE");
    if(sizeof(double)!=8)pices_fail(E,"probe requires 8-byte double");
    for(n=1;n<=nn;n++) {
        for(a=1;a<=3;a++)E->sphere.cap[1].V[a][n]=0;
        if(profile && !strcmp(profile,"constant"))E->T[1][n]=.5;
        if(profile && !strcmp(profile,"cartesian_linear"))E->T[1][n]=.5+.2*E->x[1][1][n];
    }
    if(profile && strcmp(profile,"keep")) {
        for(p=1;p<=E->trace.ntracers[1];p++) {
            double t=0;tracer_temperature_weights(E,1,p,nodes,w);
            for(a=1;a<=8;a++)t+=w[a]*E->T[1][nodes[a]];
            E->trace.extraq[1][E->pices.slot][p]=t;
        }
    }
    assemble(E,E->T[1]);
    project(E,E->trace.extraq[1][E->pices.slot],raw,0);
    project(E,E->trace.extraq[1][E->pices.slot],bc,1);
    energies[0]=energy(E,E->T[1]);energies[1]=energy(E,raw);energies[2]=energy(E,bc);
    for(e=1;e<=E->lmesh.nel;e++)for(a=1;a<=8;a++)
        local_mass[E->ien[1][e].node[a]]+=E->pices.emass[e*8+a-1];
    snprintf(path,sizeof(path),"probe.%d.nodes.bin",E->parallel.me);f=fopen(path,"wb");
    if(!f)pices_fail(E,"cannot write probe nodes");
    for(n=1;n<=nn;n++) {
        for(a=0;a<3;a++)row[a]=E->x[1][a+1][n];
        row[3]=E->T[1][n];row[4]=raw[n];row[5]=bc[n];row[6]=fixed(E,n);row[7]=local_mass[n];
        if(fwrite(row,sizeof(double),8,f)!=8)pices_fail(E,"short probe node write");
    }
    fclose(f);
    snprintf(path,sizeof(path),"probe.%d.particles.bin",E->parallel.me);f=fopen(path,"wb");
    if(!f)pices_fail(E,"cannot write probe particles");
    for(p=1;p<=E->trace.ntracers[1];p++) {
        tracer_temperature_weights(E,1,p,nodes,w);row[0]=E->trace.extraq[1][E->pices.slot][p];
        for(a=1;a<=8;a++){row[2*a-1]=nodes[a];row[2*a]=w[a];}
        if(fwrite(row,sizeof(double),17,f)!=17)pices_fail(E,"short probe particle write");
    }
    fclose(f);
    for(n=1;n<=nn;n++)error=fmax(error,fabs(bc[n]-E->T[1][n]));
    error=reduce(E,error,MPI_MAX);
    if(E->parallel.me==0) {
        f=fopen("probe.energy.txt","w");if(!f)pices_fail(E,"cannot write probe energy");
        fprintf(f,"%.17g %.17g %.17g %.17g\n",energies[0],energies[1],energies[2],error);fclose(f);
    }
    free(raw);free(bc);free(local_mass);
}
