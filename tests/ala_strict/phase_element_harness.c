/* Fixture for unmodified production geometry, phase buoyancy and thermal
 * residual functions inserted ahead of this file by test_phase_element.py.
 * Conductivity is held constant; no MPI mesh/solver or boundary flux is run. */
void audit_element(double rin, double rout, double theta_half, double phi_half,
                   const double *reference, const double *parameters,
                   double temperature_shift, double temperature_rate,
                   double radial_velocity, double conductivity, int mask,
                   double *out)
{
    static struct All_variables state;
    struct All_variables *E=&state;
    struct IEN ien[2];
    struct Shape_function GN, PG;
    struct COORD eco[2];
    struct Shape_function_dx GNx;
    struct Shape_function_dA volume;
    struct SOURCES sources;
    double xyz[4][9], spherical[4][9], T[9], Tdot[9], buoyancy[9];
    double *fields[NCS]={0}, *rates[NCS]={0}, *buoy[NCS]={0};
    double rho[3], gravity[3], cp[3], tref[3], heating[2]={0};
    double phase_heating[2]={0}, rtf[4][9], residual[9], base[9];
    float fraction[PHASE_TRANSITIONS][9], boundary[PHASE_TRANSITIONS][5];
    float VV[4][9]={{0}};
    unsigned int flags[9]={0};
    int j,n,d,p,i,side;
    memset(E,0,sizeof(*E)); memset(&sources,0,sizeof(sources));
    E->mesh.nsd=3; E->mesh.dof=3; E->mesh.levmax=0;
    E->lmesh.nno=8; E->lmesh.noz=2; E->lmesh.nox=2; E->lmesh.noy=2;
    E->lmesh.nel=1; E->lmesh.elz=1;
    E->sphere.caps_per_proc=1; E->sphere.ro=1.;
    E->ien[1]=ien; E->IEN[0][1]=ien; E->node[1]=flags;
    E->eco[1]=eco;
    eco[1].size[1]=(rin+rout)*theta_half;
    eco[1].size[2]=(rin+rout)*sin(1.1)*phi_half;
    eco[1].size[3]=rout-rin;
    E->T[1]=T; fields[1]=T; rates[1]=Tdot; buoy[1]=buoyancy;
    E->refstate.rho=rho; E->refstate.gravity=gravity;
    E->refstate.heat_capacity=cp; E->refstate.Tref=tref;
    E->heating_phase[1]=phase_heating;
    E->heating_adi[1]=heating; E->heating_visc[1]=heating;
    E->control.surface_temp=parameters[21];
    E->control.disptn_number=1.;
    for(side=1;side<=2;side++) {
        rho[side]=reference[4*(side-1)];
        gravity[side]=reference[4*(side-1)+1];
        cp[side]=reference[4*(side-1)+2];
        tref[side]=reference[4*(side-1)+3];
    }
    for(d=1;d<=3;d++) { E->X[0][1][d]=xyz[d]; E->sx[1][d]=spherical[d]; }
    for(j=1;j<=8;j++) {
        double r,theta,phi;
        side=bb[2][j]; n=2*((j-1)%4)+side;
        ien[1].node[j]=n;
        r=side==1 ? rin : rout;
        theta=1.1+(2*bb[0][j]-3)*theta_half;
        phi=.4+(2*bb[1][j]-3)*phi_half;
        xyz[1][n]=r*sin(theta)*cos(phi);
        xyz[2][n]=r*sin(theta)*sin(phi); xyz[3][n]=r*cos(theta);
        spherical[1][n]=theta; spherical[2][n]=phi; spherical[3][n]=r;
        T[n]=tref[side]+temperature_shift; Tdot[n]=temperature_rate;
        VV[3][j]=radial_velocity; buoyancy[n]=0.;
    }
    for(p=0;p<PHASE_TRANSITIONS;p++) {
        E->control.phase[p].depth=parameters[7*p];
        E->control.phase[p].clapeyron=parameters[7*p+1];
        E->control.phase[p].inv_width=1./parameters[7*p+2];
        E->control.phase[p].transT=parameters[7*p+3];
        E->control.phase[p].entropy_jump=0.;
        E->control.phase[p].density_jump=parameters[7*p+5];
        E->control.phase[p].Ra=parameters[7*p+6];
        E->phase_B[p][1]=fraction[p]; E->phase_boundary[p][1]=boundary[p];
    }
    construct_shape_functions(E);
    get_global_shape_fn(E,1,&GN,&GNx,&volume,0,1,rtf,0,1);
    pg_shape_fn(E,1,&PG,&GNx,VV,rtf,conductivity,1);
    element_residual(E,1,PG,GNx,volume,VV,fields,rates,sources,
                     base,rtf,conductivity,NULL,NULL,1);
    for(p=0;p<PHASE_TRANSITIONS;p++)
        if(mask & (1<<p)) E->control.phase[p].entropy_jump=parameters[7*p+4];
    element_residual(E,1,PG,GNx,volume,VV,fields,rates,sources,
                     residual,rtf,conductivity,NULL,NULL,1);
    phase_change_apply(E,buoy);
    for(i=1;i<=8;i++) {
        double x[4]={0}, temp=0., density=0., rho_g=0., d_rg=0., gradT=0.;
        for(j=1;j<=8;j++) {
            double shape=E->N.vpt[GNVINDEX(j,i)];
            n=ien[1].node[j]; side=(n-1)%2+1;
            for(d=1;d<=3;d++) x[d]+=shape*xyz[d][n];
            temp+=shape*T[n]; density+=shape*rho[side];
            rho_g+=shape*rho[side]*gravity[side];
            d_rg+=GNx.vpt[GNVXINDEX(2,j,i)]*rho[side]*gravity[side];
            gradT+=GNx.vpt[GNVXINDEX(2,j,i)]*T[n];
        }
        out[8*(i-1)]=rtf[3][i];
        out[8*(i-1)+1]=sqrt(x[1]*x[1]+x[2]*x[2]+x[3]*x[3]);
        out[8*(i-1)+2]=temp; out[8*(i-1)+3]=density;
        out[8*(i-1)+4]=rho_g; out[8*(i-1)+5]=d_rg;
        out[8*(i-1)+6]=gradT; out[8*(i-1)+7]=volume.vpt[i];
        out[64+i-1]=residual[i]-base[i];
        out[72+i-1]=buoyancy[i];
    }
    out[80]=phase_heating[1];
    for(p=0;p<PHASE_TRANSITIONS;p++) {
        for(n=1;n<=8;n++) {
            out[81+16*p+n-1]=fraction[p][n];
            out[81+16*p+8+n-1]=phase_change_reference_fraction(E,p,1,n);
        }
    }
}
