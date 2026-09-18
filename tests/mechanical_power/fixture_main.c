int main(int argc,char **argv) {
    struct All_variables state={0},*E=&state;
    struct Npz_writer writer;
    int i,j,k,nz,mode=atoi(argv[2]);
    double raw[9],avg=0;
    MPI_Init(&argc,&argv);
    E->parallel.world=MPI_COMM_WORLD;
    E->parallel.nprocz=1; E->parallel.me_loc[3]=0;
    E->sphere.caps_per_proc=1;
    E->mesh.nsd=3; E->mesh.levmax=0; E->mesh.noz=2;
    E->lmesh.nno=8; E->lmesh.noz=2; E->lmesh.nel=1;
    E->lmesh.npno=1; E->lmesh.elz=1; E->lmesh.neq=24;
    E->lmesh.NEL[0]=1; E->lmesh.NPNO[0]=1; E->lmesh.NEQ[0]=24;
    E->control.eba_formulation=1; E->control.Atemp=4;
    E->control.disptn_number=mode==3 ? 0 : 2;
    /* Nonzero inv_gruneisen must not introduce ALA beta into EBA. */
    E->control.inv_gruneisen=7;
    E->monitor.solution_cycles=12; E->monitor.elapsed_time=0.125;
    E->fp=stdout; E->solver.exchange_id_d=exchange;
    E->ien[1]=calloc(2,sizeof(struct IEN)); E->IEN[0][1]=E->ien[1];
    E->id[1]=calloc(9,sizeof(struct ID)); E->ID[0][1]=E->id[1];
    E->elt_del[0][1]=calloc(2,sizeof(struct EG));
    E->node[1]=calloc(9,sizeof(unsigned int));
    E->U[1]=calloc(25,sizeof(double)); E->P[1]=calloc(2,sizeof(double));
    E->buoyancy[1]=calloc(9,sizeof(double));
    E->T[1]=calloc(9,sizeof(double)); E->EVi[1]=calloc(9,sizeof(float));
    E->eco[1]=calloc(2,sizeof(*E->eco[1])); E->eco[1][1].area=1;
    E->refstate.rho=calloc(3,sizeof(double));
    E->refstate.gravity=calloc(3,sizeof(double));
    E->refstate.thermal_expansivity=calloc(3,sizeof(double));
    E->boundary.nel=1; E->boundary.element[1]=calloc(2,sizeof(int));
    E->boundary.element[1][1]=1; E->P[1][1]=0.9;
    E->control.tracer=1; E->composition.ichemical_buoyancy=1;
    E->composition.ncomp=1;
    E->composition.buoyancy_ratio=calloc(1,sizeof(double));
    E->composition.buoyancy_ratio[0]=0.3;
    E->composition.comp_node[1]=calloc(1,sizeof(double *));
    E->composition.comp_node[1][0]=calloc(9,sizeof(double));
    for(nz=1;nz<=2;nz++) {
        E->refstate.rho[nz]=1.2; E->refstate.gravity[nz]=1.1;
        E->refstate.thermal_expansivity[nz]=0.8;
    }
    for(k=0;k<3;k++) {
        E->control.phase[k].Ra=0.2*(k+1);
        E->phase_B[k][1]=calloc(9,sizeof(float));
    }
    for(i=1;i<=8;i++) {
        E->ien[1][1].node[i]=i; E->node[1][i]=VBX|VBY|VBZ;
        E->T[1][i]=0.1*i; E->composition.comp_node[1][0][i]=0.03*i;
        for(j=1;j<=3;j++) {
            int eq=(i-1)*3+j-1;
            E->id[1][i].doff[j]=eq; E->U[1][eq]=0.2+eq*0.01;
            E->elt_del[0][1][1].g[eq][0]=0.02*(eq+1);
        }
        E->N.vpt[GNVINDEX(i,i)]=1; E->EVi[1][i]=2;
        raw[i]=4*1.2*0.8*E->T[1][i]-4*0.3*0.03*i;
        for(k=0;k<3;k++) {
            E->phase_B[k][1][i]=0.05*(i+k);
            raw[i]-=E->control.phase[k].Ra*(E->phase_B[k][1][i]-0.25f);
        }
        raw[i]*=1.1; avg+=raw[i]/8;
    }
    for(i=1;i<=8;i++) E->buoyancy[1][i]=raw[i]-avg;
    if(mode==1) E->node[1][1]&=~VBX;
    if(mode==4) E->P[1][1]+=1000;
    if(mode==5) {
        for(i=0;i<24;i++) E->U[1][i]=0;
    }
    if(mode!=2) write_eba_mechanical_power(E);
    /* Verify output uses the solve's metadata rather than the later state. */
    E->monitor.solution_cycles=13; E->monitor.elapsed_time=0.25;
    if(npz_open(&writer,argv[1]) || write_eba_power_npz(&writer,E)
       || npz_close(&writer)) return 2;
    MPI_Finalize(); return 0;
}
