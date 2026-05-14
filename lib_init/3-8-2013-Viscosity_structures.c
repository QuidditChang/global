/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 *<LicenseText>
 *
 * CitcomS by Louis Moresi, Shijie Zhong, Lijie Han, Eh Tan,
 * Clint Conrad, Michael Gurnis, and Eun-seo Choi.
 * Copyright (C) 1994-2005, California Institute of Technology.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *</LicenseText>
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
/* Functions relating to the determination of viscosity field either
   as a function of the run, as an initial condition or as specified from
   a previous file */


#include <math.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "parsing.h"


void myerror(struct All_variables *,char *);

static void apply_low_visc_wedge_channel(struct All_variables *E, float **evisc);
static void low_viscosity_channel_factor(struct All_variables *E, float *F);
static void low_viscosity_wedge_factor(struct All_variables *E, float *F);
void parallel_process_termination();
//void weak_plate_boundary_from_velocities();


void viscosity_system_input(struct All_variables *E)
{
    int m=E->parallel.me;
    int i;

    /* default values .... */
    for(i=0;i<40;i++) {
        E->viscosity.N0[i]=1.0;
        E->viscosity.T[i] = 0.0;
        E->viscosity.Z[i] = 0.0;
        E->viscosity.E[i] = 0.0;

	E->viscosity.pdepv_a[i] = 1.e20; /* \sigma_y = min(a + b * (1-r),y) */
	E->viscosity.pdepv_b[i] = 0.0;
	E->viscosity.pdepv_y[i] = 1.e20;


    }
    for(i=0;i<10;i++)
      E->viscosity.cdepv_ff[i] = 1.0; /* flavor factors for CDEPV */


    /* read in information */
    input_boolean("VISC_UPDATE",&(E->viscosity.update_allowed),"on",m);
    input_int("rheol",&(E->viscosity.RHEOL),"3",m);
    input_int("num_mat",&(E->viscosity.num_mat),"1",m);
    input_float_vector("visc0",E->viscosity.num_mat,(E->viscosity.N0),m);

    input_boolean("TDEPV",&(E->viscosity.TDEPV),"on",m);
    if (E->viscosity.TDEPV) {
        input_float_vector("viscT",E->viscosity.num_mat,(E->viscosity.T),m);
        input_float_vector("viscE",E->viscosity.num_mat,(E->viscosity.E),m);
        input_float_vector("viscZ",E->viscosity.num_mat,(E->viscosity.Z),m);
    }


    E->viscosity.sdepv_misfit = 1.0;
    input_boolean("SDEPV",&(E->viscosity.SDEPV),"off",m);
    if (E->viscosity.SDEPV) {
      E->viscosity.sdepv_visited = 0;
      input_float_vector("sdepv_expt",E->viscosity.num_mat,(E->viscosity.sdepv_expt),m);
    }


    input_boolean("PDEPV",&(E->viscosity.PDEPV),"off",m); /* plasticity addition by TWB */
    if (E->viscosity.PDEPV) {
      E->viscosity.pdepv_visited = 0;

      input_boolean("pdepv_eff",&(E->viscosity.pdepv_eff),"on",m);
      input_float_vector("pdepv_a",E->viscosity.num_mat,(E->viscosity.pdepv_a),m);
      input_float_vector("pdepv_b",E->viscosity.num_mat,(E->viscosity.pdepv_b),m);
      input_float_vector("pdepv_y",E->viscosity.num_mat,(E->viscosity.pdepv_y),m);

      input_float("pdepv_offset",&(E->viscosity.pdepv_offset),"0.0",m);
    }
    if(E->viscosity.PDEPV || E->viscosity.SDEPV)
      input_float("sdepv_misfit",&(E->viscosity.sdepv_misfit),"0.001",m);


    input_boolean("CDEPV",&(E->viscosity.CDEPV),"off",m);
    for(i=0;i<10;i++)
      E->viscosity.cdepv_ff[i] = 1.0; /* flavor factors for CDEPV */
    if(E->viscosity.CDEPV){
      /* compositional viscosity */
      if(E->control.tracer < 1){
        fprintf(stderr,"error: CDEPV requires tracers, but tracer is off\n");
        parallel_process_termination();
      }
      if(E->trace.nflavors > 10)
        myerror(E,"error: too many flavors for CDEPV");
      /* read in flavor factors */
      input_float_vector("cdepv_ff",E->trace.nflavors,
                         (E->viscosity.cdepv_ff),m);
      /* and take the log because we're using a geometric avg */
      for(i=0;i<E->trace.nflavors;i++)
        E->viscosity.cdepv_ff[i] = log(E->viscosity.cdepv_ff[i]);
    }


    input_boolean("low_visc_channel",&(E->viscosity.channel),"off",m);
    input_boolean("low_visc_wedge",&(E->viscosity.wedge),"off",m);

    input_float("lv_min_radius",&(E->viscosity.lv_min_radius),"0.9764",m);
    input_float("lv_max_radius",&(E->viscosity.lv_max_radius),"0.9921",m);
    input_float("lv_channel_thickness",&(E->viscosity.lv_channel_thickness),"0.0047",m);
    input_float("lv_reduction",&(E->viscosity.lv_reduction),"0.5",m);

    input_boolean("VMAX",&(E->viscosity.MAX),"off",m);
    if (E->viscosity.MAX)
        input_float("visc_max",&(E->viscosity.max_value),"1e22,1,nomax",m);

    input_boolean("VMIN",&(E->viscosity.MIN),"off",m);
    if (E->viscosity.MIN)
        input_float("visc_min",&(E->viscosity.min_value),"1e20",m);

    return;
}


void viscosity_input(struct All_variables *E)
{
    int m = E->parallel.me;

    input_string("Viscosity",E->viscosity.STRUCTURE,"system",m);
    input_int ("visc_smooth_method",&(E->viscosity.smooth_cycles),"0",m);

    if ( strcmp(E->viscosity.STRUCTURE,"system") == 0)
        E->viscosity.FROM_SYSTEM = 1;
    else
        E->viscosity.FROM_SYSTEM = 0;

    if (E->viscosity.FROM_SYSTEM)
        viscosity_system_input(E);

    return;
}



/* ============================================ */

void get_system_viscosity(E,propogate,evisc,visc)
     struct All_variables *E;
     int propogate;
     float **evisc,**visc;
{
    void visc_from_mat();
    void visc_from_T();
    void visc_from_S();

    void visc_from_P();
    void visc_from_C();

    void apply_viscosity_smoother();
    void visc_from_gint_to_nodes();



    int i,j,m;
    float temp1,temp2,*vvvis;
    double *TG;

    const int vpts = vpoints[E->mesh.nsd];

    if(E->viscosity.TDEPV)
        visc_from_T(E,evisc,propogate);
    else
        visc_from_mat(E,evisc);

    if(E->viscosity.CDEPV) {	/* compositional prefactor */
      //fprintf(E->fp,"before visc_from_C!\n");
      visc_from_C(E,evisc);
    }

    if(E->viscosity.SDEPV)
      visc_from_S(E,evisc,propogate);

    if(E->viscosity.PDEPV)	/* "plasticity" */
      visc_from_P(E,evisc);


    /* i think this should me placed differently i.e.  before the
       stress dependence but I won't change it because it's by
       someone else

       TWB
    */
    if(E->viscosity.channel || E->viscosity.wedge)
        apply_low_visc_wedge_channel(E, evisc);


    /* min/max cut-off */

    if(E->viscosity.MAX) {
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(i=1;i<=E->lmesh.nel;i++)
                for(j=1;j<=vpts;j++)
                    if(evisc[m][(i-1)*vpts + j] > E->viscosity.max_value)
                        evisc[m][(i-1)*vpts + j] = E->viscosity.max_value;
    }

    if(E->viscosity.MIN) {
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(i=1;i<=E->lmesh.nel;i++)
                for(j=1;j<=vpts;j++)
                    if(evisc[m][(i-1)*vpts + j] < E->viscosity.min_value)
                        evisc[m][(i-1)*vpts + j] = E->viscosity.min_value;
    }

    if (E->control.verbose)  {
      fprintf(E->fp_out,"output_evisc \n");
      for(m=1;m<=E->sphere.caps_per_proc;m++) {
        fprintf(E->fp_out,"output_evisc for cap %d\n",E->sphere.capid[m]);
      for(i=1;i<=E->lmesh.nel;i++)
          fprintf(E->fp_out,"%d %d %f %f\n",i,E->mat[m][i],evisc[m][(i-1)*vpts+1],evisc[m][(i-1)*vpts+7]);
      }
      fflush(E->fp_out);
    }

    /* interpolate from gauss quadrature points to node points for output */
    visc_from_gint_to_nodes(E,evisc,visc,E->mesh.levmax);
    fprintf(E->fp,"after visc_from_gint_to_nodes\n");

    return;
}



void initial_viscosity(struct All_variables *E)
{
    if (E->viscosity.FROM_SYSTEM)
        get_system_viscosity(E,1,E->EVI[E->mesh.levmax],E->VI[E->mesh.levmax]);

    return;
}


void visc_from_mat(E,EEta)
     struct All_variables *E;
     float **EEta;
{

    int i,m,jj;

    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(i=1;i<=E->lmesh.nel;i++)
            for(jj=1;jj<=vpoints[E->mesh.nsd];jj++)
                EEta[m][ (i-1)*vpoints[E->mesh.nsd]+jj ]=E->viscosity.N0[E->mat[m][i]-1];

    return;
}

void visc_from_T(E,EEta,propogate)
     struct All_variables *E;
     float **EEta;
     int propogate;
{
    int m,i,j,k,l,z,jj,kk,imark,iii,jjj,kkk,intage,node,nodeg,nnn,ttt,count,bd_i,bd_j,bd_k;
    int theta_start,fi_start,rheo_trick,cap,flag,flag1,flag2,*PB3[4],cut_ocean_plate;
    float zero,e_6,one,eta0,Tave,depth,temp,tempa,temp1,temp2,TT[9],lon_inc,lon_inc1,lon_dcr;
    float zzz,zz[9],lon_temp,r_temp,r_min,lon_max,guide_dep,dip,dist1,d_wdg;
    float visc1, visc2, tempa_exp,*PB1[4],*PB2[4];
    float phi, theta, rad, fi, r, fi_trench, the_trench,v_trench;
    float find_age_in_MY();
    float T_max, T_min,age,age1,age2,viscE,viscT;
    char input_s[1000],pb_1[255],pb_2[255],bounds[255];
    double Tmaxd();
    const int vpts = vpoints[E->mesh.nsd];
    const int ends = enodes[E->mesh.nsd];
    const int nel = E->lmesh.nel;
    const int dims=E->mesh.nsd;
    FILE *fp,*fp1,*fp2,*fp3,*fp4;

    void weak_plate_boundary_from_velocities();

    e_6 = 1.e-6;
    one = 1.0;
    zero = 0.0;
    imark = 0;

    switch (E->viscosity.RHEOL)   {
    case 1:			/* eta = N_0 exp( E * (T_0 - T))  */
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(i=1;i<=nel;i++)   {
                l = E->mat[m][i];

                if(E->control.mat_control==0)
                    tempa = E->viscosity.N0[l-1];
                else if(E->control.mat_control==1)
                    tempa = E->viscosity.N0[l-1]*E->VIP[m][i];

                for(kk=1;kk<=ends;kk++) {
                    TT[kk] = E->T[m][E->ien[m][i].node[kk]];
                }

                for(jj=1;jj<=vpts;jj++) {
                    temp=0.0;
                    for(kk=1;kk<=ends;kk++)   {
                        temp += TT[kk] * E->N.vpt[GNVINDEX(kk,jj)];
                    }

                    EEta[m][ (i-1)*vpts + jj ] = tempa*
                        exp( E->viscosity.E[l-1] * (E->viscosity.T[l-1] - temp));

                }
            }
        break;

    case 2:			/* eta = N_0 exp(-T/T_0) */
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(i=1;i<=nel;i++)   {
                l = E->mat[m][i];

                if(E->control.mat_control==0)
                    tempa = E->viscosity.N0[l-1];
                else if(E->control.mat_control==1)
                    tempa = E->viscosity.N0[l-1]*E->VIP[m][i];

                for(kk=1;kk<=ends;kk++) {
                    TT[kk] = E->T[m][E->ien[m][i].node[kk]];
                }

                for(jj=1;jj<=vpts;jj++) {
                    temp=0.0;
                    for(kk=1;kk<=ends;kk++)   {
                        temp += TT[kk] * E->N.vpt[GNVINDEX(kk,jj)];
                    }

                    EEta[m][ (i-1)*vpts + jj ] = tempa*
                        exp( -temp / E->viscosity.T[l-1]);

                }
            }
        break;

    case 3:			/* eta = N_0 exp(E/(T+T_0) - E/(1+T_0)) */
	//fprintf(E->fp,"emax=%d\n",E->mesh.nel);
	//fprintf(E->fp,"nel=%d\n",E->mesh.elx*E->mesh.elz*E->mesh.ely);
	fp=fopen("filter_limit","r");
        if(!fp) {
              T_max=Tmaxd(E,E->T);
              T_min=Tmind(E,E->T);
              fp=fopen("filter_limit","w");
              fprintf(fp,"%f %f\n",T_max,T_min);
              fclose(fp);
        }
        else {
              fgets(input_s,1000,fp);
              sscanf(input_s,"%f %f",&T_max,&T_min);
              fclose(fp);
        }
        //one = T_max;
	one = E->control.lith_age_mantle_temp;

	//temporaray change
        T_max = one; T_min = 0.0;

        age=find_age_in_MY(E);
        intage = age;
        age1=1.0*age;
        age2=age1+1.0;

        fp2=fopen("rheo.dat","r");
        if(fp2==NULL) 
	    rheo_trick = 0;
	else {
            fscanf(fp2,"%d",&(rheo_trick));
            fclose(fp2);
	}

        fp3=fopen("velo_trench","r");
        if(fp3==NULL)
            v_trench = 2.0;
        else {
            fscanf(fp3,"%f",&(v_trench));
            fclose(fp3);
        }     
        if(E->parallel.me==0) {
	    fprintf(E->fp,"v_trench=%f\n",v_trench);
	}

	for(m=1;m<=E->sphere.caps_per_proc;m++)
          for(kkk=1;kkk<=E->lmesh.noy;kkk++)
            for(jjj=1;jjj<=E->lmesh.nox;jjj++){
                if(E->sphere.caps > 1 && rheo_trick != 0) {
                    flag = 0;
                    flag1 = 0;
                    }

		flag2 = 0;
                for(iii=E->lmesh.noz;iii>=1;iii--)  {
                    if(kkk == 1)
                         bd_k = 1;
                    else
                         bd_k = 2;
                    if(jjj == 1)
                         bd_j = 1;
                    else
                         bd_j = 2;
                    if(iii == 1)
                         bd_i = 0;
                    else
                         bd_i = 1;

                    i = (iii-bd_i) + (jjj-bd_j)*(E->lmesh.noz-1)
                                + (kkk-bd_k)*(E->lmesh.noz-1)*(E->lmesh.nox-1);

                    node = iii + (jjj-1)*E->lmesh.noz
                             + (kkk-1)*E->lmesh.noz*E->lmesh.nox;

		    nodeg=E->lmesh.nxs-1+jjj+(E->lmesh.nys+kkk-2)*E->mesh.nox;

                    theta = E->sx[m][1][node];
                    phi = E->sx[m][2][node];
                    r = E->sx[m][3][node];

		    l = E->mat[m][i];
                    tempa = E->viscosity.N0[l-1];

		    /* creat a smooth radial visc. profile */
                    if(r>0.993){
                        tempa = E->viscosity.N0[0] + (E->viscosity.N0[1]-E->viscosity.N0[0])*(cos((r-0.993)/0.007*3.14159)+1)/2;
                        //tempa = E->viscosity.N0[0] + (E->viscosity.N0[1]-E->viscosity.N0[0])*(cos((r-0.99)/0.01*3.14159)+1)/2;
                        viscE = E->viscosity.E[0] + (E->viscosity.E[1]-E->viscosity.E[0])*(cos((r-0.993)/0.007*3.14159)+1)/2;
                        viscT = E->viscosity.T[0] + (E->viscosity.T[1]-E->viscosity.T[0])*(cos((r-0.993)/0.007*3.14159)+1)/2;
                    }
                    else if(r<0.97 && r>0.95) {
		    //else if(r<0.95 && r>0.92) {
                        tempa = E->viscosity.N0[1] + (E->viscosity.N0[2]-E->viscosity.N0[1])*(cos((r-0.95)/0.02*3.14159)+1)/2;
                        viscE = E->viscosity.E[1] + (E->viscosity.E[2]-E->viscosity.E[1])*(cos((r-0.95)/0.02*3.14159)+1)/2;
                        viscT = E->viscosity.T[1] + (E->viscosity.T[2]-E->viscosity.T[1])*(cos((r-0.95)/0.02*3.14159)+1)/2;
                    }
                    else if(r<0.91 && r>0.85) {
                        tempa = E->viscosity.N0[2] + (E->viscosity.N0[3]-E->viscosity.N0[2])*(cos((r-0.85)/0.06*3.14159)+1)/2;
                        viscE = E->viscosity.E[2] + (E->viscosity.E[3]-E->viscosity.E[2])*(cos((r-0.85)/0.06*3.14159)+1)/2;
                        viscT = E->viscosity.T[2] + (E->viscosity.T[3]-E->viscosity.T[2])*(cos((r-0.85)/0.06*3.14159)+1)/2;
                    }
                    else {
                        tempa = E->viscosity.N0[l-1];
                        viscE = E->viscosity.E[l-1];
                        viscT = E->viscosity.T[l-1];
                    }
                    // end of smooth


                    j = 0;

                    for(kk=1;kk<=ends;kk++) {
                        TT[kk] = E->T[m][E->ien[m][i].node[kk]];
                        zz[kk] = (1.-E->sx[m][3][E->ien[m][i].node[kk]]);
                    }

                    for(jj=1;jj<=vpts;jj++) {

                        temp=0.0;
                        zzz=0.0;
                        for(kk=1;kk<=ends;kk++)   {
                            TT[kk]=max(TT[kk],zero);
                            //temp += min(TT[kk],one) * E->N.vpt[GNVINDEX(kk,jj)];
                            temp += min(TT[kk],T_max) * E->N.vpt[GNVINDEX(kk,jj)];
                            zzz += zz[kk] * E->N.vpt[GNVINDEX(kk,jj)];
                        }

                        if(E->control.mat_control==0) {
				EEta[m][ (i-1)*vpts + jj ] = tempa*exp(viscE/(temp+viscT) - viscE/(one+viscT) );
				temp1 =500*tempa; // *(1.1-cos((r-0.9)/0.1*1.57));
                                if(EEta[m][(i-1)*vpts+jj]>temp1 && r>0.92 && r<0.993)
				   EEta[m][(i-1)*vpts+jj]= temp1; 

                        }


			if(E->control.mat_control==1) {
                                //Lijun addition: maintain the same slab stiffness for all visc_um values
                                if(r > 0.91 && r < 1.0)
                                  //if(rheo_trick == 4 && (E->control.TBCbotval-E->T[m][node]) > 0.5*(T_max-T_min) || rheo_trick == 5 && (E->control.TBCbotval-E->T[m][node]) > 0.5*(T_max-T_min))
                                  if((E->control.TBCbotval-E->T[m][node]) > 0.8*(E->control.TBCbotval-T_min))
                                        tempa = 1.0*T_max/E->control.TBCbotval;
                                //End of Lijun's addition
                                EEta[m][ (i-1)*vpts + jj ] = tempa*E->VIP[m][i]*exp( E->viscosity.E[l-1]/((temp-T_min)/(one-T_min)*one+E->viscosity.T[l-1]) - E->viscosity.E[l-1]/(one +E->viscosity.T[l-1]) );
                        }

			if(rheo_trick == 1) { //a strong slab in the upper mantle as a stress guide
                                temp1 = 100.0;
                                temp2 = 0.1;
                                if(flag > 0 && E->sphere.caps > 1) {
                                        fi = PB1[2][flag]+2*3.1415926;
                                        lon_inc = 0.7;
                                        if (E->sx[m][2][node]>fi-0.1 && E->sx[m][2][node]<fi+lon_inc && r>0.91) {
                                            if (E->sx[m][2][node]>fi && E->sx[m][2][node]<fi+0.06 && r>0.99)
                                                EEta[m][ (i-1)*vpts + jj ] = temp2;
                                            if (E->sx[m][2][node]>fi-0.1 && temp<0.9*E->control.TBCbotval && r<0.975)
                                                if (E->sx[m][2][node]<fi-0.02 && flag2<=1*vpts || E->sx[m][2][node]>fi-0.02 && E->sx[m][2][node]<fi+0.03 && flag2<=2*vpts || E->sx[m][2][node]>fi+0.03 && flag2<=3*vpts) {
                                                    flag2 += 1;
                                                    EEta[m][ (i-0)*vpts + jj ] = temp2;
                                                }
                                            if (EEta[m][ (i+1)*vpts + jj ] == temp2 && EEta[m][ (i)*vpts + jj ] != temp2) {
                                                EEta[m][ (i-1)*vpts + jj ] = temp1;
                                                EEta[m][ (i-2)*vpts + jj ] = temp1;
                                            }
                                        }
                                }
                        }

			else if(rheo_trick == 2) { //lateral viscosity variation underneath continent
                                temp1 = 100.0;
                                temp2 = 0.1;
                                fi = PB1[2][flag]+2*3.1415926;
                                if(r>0.97 && E->sx[m][2][node]>=fi && E->sx[m][2][node]<=fi+PB1[3][flag] && age >50) {
                                        EEta[m][ (i)*vpts + jj ] = temp1;
                                }
                        }
			
			else if(rheo_trick >= 3) { // stress guide + lateral viscosity variation
                                temp1 = 100.0;
                                temp2 = 0.1;
                                if(flag > 0 && E->sphere.caps > 1) {
                                        fi = PB1[2][flag]+2*3.1415926-0.06;
                                        lon_inc = 1.5;
                                        //lon_inc = 0.5;
                                        lon_inc1 = lon_inc;
                                        lon_dcr = 0.15;
                                        guide_dep = 0.91;


					if (E->sx[m][2][node]>fi-lon_dcr && E->sx[m][2][node]<fi+0.5) {
					    if(r > 0.98)
                                                EEta[m][ (i-1)*vpts + jj ] = temp1;	
					    if(E->sphere.cap[m].V[3][node] > 500.0 && E->T[m][node] <= 1.0 && r > 0.95)
					            E->T[m][node] += (1.0-E->control.TBCbotval)*(E->sphere.cap[m].V[3][node]/10000.0)*(0.001); //melt ~ vertical velo. * timestep
                                        }

			
					if (E->sx[m][2][node]>fi && E->sx[m][2][node]<fi+lon_inc1+0.03 && r>0.99){
                                                EEta[m][ (i-1)*vpts + jj ] = temp2;
                                        } //end of weak layer definition

                                } //end of caps > 1

				else if(E->sphere.caps == 1) {
					    temp2 = 7.0; 
					    //if(v_trench>2.0) temp2 = 7.0;
					    //if(age<20.0 && age>=15.0) 
						//v_trench += (20.0-age)/(20.0-15.0)*(-6.0-v_trench);
					    //if(age<15.0) v_trench = -6.0;
					    fi = 0.5+v_trench*(3.14159/(11.1*180.0))*(200.0-age); //v_trench = 2.0cm/yr
					    //fi = 1.5+v_trench*(3.14159/(11.1*180.0))*(200.0-age); //v_trench = 2.0cm/yr

					    //fi=4.4+v_trench*(3.14159/(11.1*180.0))*(40.0-age);

					    //if(age<191.0) 
						//fi -= 2.0*(3.14159/(11.1*180.0))*(191.0-age);
					    //if(E->parallel.me==0)
					    //	fprintf(E->fp,"fi=%f\n",fi);
                                            temp1 = phi-((fi+0.0)+temp2*(1.0-r));
                                            dist1 = 0.08;
                                            if(temp1>-dist1 && temp1<0.0 && r>0.985) { // && theta<=lat_max+0.015 && theta>=lat_min+0.005) {
						temp2 = 0.01+(1.0-(1*cos((pow(temp1/dist1,0.3)*1.57))));
                                                if(EEta[m][(i-1)*vpts+jj]>temp2) { // && theta<=lat_max-0.005 && theta>=lat_min+0.0) {
                                                  if(r>0.99)
                                                     EEta[m][(i-1)*vpts+jj]=temp2*1*cos((r-0.99)/0.01*1.57);
                                                  else {
                                                     temp1 = EEta[m][(i-1)*vpts+jj];
                                                     temp2 = temp2*1 + (temp1-temp2*1)*(0.99-r)/(0.99-0.98);
                                                     EEta[m][(i-1)*vpts+jj] = temp2;
                                                  }
                                                }
                                            } //end of weak top

					    temp2 = 5.0;
					    if(v_trench<-2.0) temp2 = 7.0;
                                            temp1 = phi-((fi-0.0)+temp2*(1.0-r));
                                            dist1 = 0.03;
                                            if(temp1>0.0-dist1 && temp1<0.04 && r>0.985 && phi>fi-0.03){
                                                if(temp1<0.0) temp2 = (cos(pow(temp1/dist1,3)*3.1416)+1)/2;
                                                if(temp1>=0.0) temp2 = (cos(pow(temp1/0.05,3)*3.1416)+1)/2;
						//if(temp1>=0.0) temp2 = (cos(pow(temp1/0.03,1)*3.1416)+1)/2;
                                                temp2 = 0.001 + (1.0 - temp2);
                                                if(EEta[m][(i-1)*vpts+jj]>temp2) { // && theta<=lat_max-0.005 && theta>=lat_min+0.0) {
                                                  if(temp1>0.0) // && r>0.99)
                                                      EEta[m][(i-1)*vpts+jj]=temp2;
                                                  if(phi>fi && phi<fi+0.03 && r>0.99 && phi-((fi-0.0)+1*(1.0-r))>0.0) {
                                                     dist1 = 1.0 - cos(pow((1.0-r)/0.01,2)*1.57);
                                                     if(EEta[m][(i-1)*vpts+jj]>temp2) EEta[m][(i-1)*vpts+jj] = dist1;
                                                  }
                                                }
                                            }  //end of weak subduction zone

					    /* add a weak mantle wedge atop slab */
                                            temp1 = phi-fi;
                                            dist1 = 0.1;
                                            dip = 1.0; d_wdg = 0.975;
                                            //if(age<15) {
                                              //  dip += (2.0-dip)*(15-age)/15.0;
                                            //}

                                            //if(theta<=lat_max+0.1 && theta>=lat_min-0.1 && phi-(fi+dip*(1.0-r))>0.0) {
					    if(phi-(fi+dip*(1.0-r))>0.0){
                                                if(temp1>=0.0 && temp1<dist1 && r>d_wdg && flag2==0 && r<1.0) { // && r<1.0
                                                  if(E->T[m][node]>=0.7*E->control.lith_age_mantle_temp) {
                                                    if(E->T[m][node]>0.9*E->control.lith_age_mantle_temp)
                                                        temp2 = 0.001*EEta[m][(i-1)*vpts+jj];
                                                    else if(E->T[m][node]>0.8*E->control.lith_age_mantle_temp)
                                                        temp2 = 0.01*EEta[m][(i-1)*vpts+jj];
                                                    else if(E->T[m][node]>=0.7*E->control.lith_age_mantle_temp)
                                                        temp2 = 0.1*EEta[m][(i-1)*vpts+jj];
                                                    if(temp2 < 0.01)
                                                        temp2 = 0.01;
                                                    if(temp1>dist1*3/4.0)
                                                        temp2 += (EEta[m][(i-1)*vpts+jj]-temp2)*cos(pow((dist1-temp1)/(dist1/4.0),0.3)*1.57);
                                                    EEta[m][(i-1)*vpts+jj] = temp2;
                                                  }
                                                  else
                                                   flag2 = 1;
                                                }
                                              } //end of mantle wedge

					    /* add a localized craton (CP) as a high-viscosity body */
                                            //dist1 = 0.13;
					    dist1 = 0.05;
                                            temp1 = phi-(fi+dist1+2*(1.0-r));
                                            /*if(r>0.97 && temp1>=0.0 && phi<=E->control.fi_max-0.3) {
                                                  temp2 = E->viscosity.N0[1] + 50*(cos((1.0-r)/0.03*3.1416)+1)/2.0;
                                                  if(r<0.98) temp2 *= 1.0-cos(pow((r-0.975)/0.005,3)*1.57);
                                                      EEta[m][(i)*vpts+jj] = temp2;
                                                  if(r<0.98)
                                                      temp2 = (cos(pow((r-0.97)/0.01,3)*3.1416)+1)/2.0;

                                                  if(EEta[m][(i)*vpts+jj]<0.01) EEta[m][(i)*vpts+jj]=0.01;

                                                  if(E->T[m][node]>=0.5*E->control.lith_age_mantle_temp) {
                                                    if(E->T[m][node]>0.8*E->control.lith_age_mantle_temp)
                                                        temp2 = 0.001*EEta[m][(i-1)*vpts+jj];
                                                    else if(E->T[m][node]>0.7*E->control.lith_age_mantle_temp)
                                                        temp2 = 0.01*EEta[m][(i-1)*vpts+jj];
                                                    else if(E->T[m][node]>=0.5*E->control.lith_age_mantle_temp)
                                                        temp2 = 0.1*EEta[m][(i-1)*vpts+jj];
                                                    if(temp2 < 0.01)
                                                        temp2 = 0.01;
                                                    EEta[m][(i-1)*vpts+jj] = temp2;
                                                  }
                                            } */ // end of CP


                                    } // end of caps == 1

                              } // end of rheo_trick >= 3

			      /* add a weak boundary surrounding the box */
                       /* temp1 = E->viscosity.N0[l-1];
                        if(l<3) 
                            temp = 0.01;
                        else 
                            temp = 0.01*temp1;
                        if(phi > E->control.fi_max-0.1)
                            if(phi > E->control.fi_max-0.05) {
                                temp2 = temp1+(temp-temp1)*(cos(pow((phi-(E->control.fi_max-0.05))/0.05,10)*3.14159)+1)/2;
                                EEta[m][(i-1)*vpts + jj] = temp2; 
                            }
			    else {
                                temp2 = temp1+(temp-temp1)*(cos(pow((phi-(E->control.fi_max-0.05))/0.05,3)*3.14159)+1)/2;
                                if(EEta[m][(i-1)*vpts + jj] > temp2)
                                  EEta[m][(i-1)*vpts + jj] = temp2;
                            }
                        if(phi < E->control.fi_min+0.1)
                           if(phi < E->control.fi_min+0.05) {
                                temp2 = temp1+(temp-temp1)*(cos(pow((E->control.fi_min+0.05-phi)/0.05,10)*3.14159)+1)/2;
                                EEta[m][(i-1)*vpts + jj] = temp2;
                           }
			   else {
                                temp2 = temp1+(temp-temp1)*(cos(pow((E->control.fi_min+0.5-phi)/0.05,3)*3.14159)+1)/2;
                                if(EEta[m][(i-1)*vpts + jj] > temp2)
                                  EEta[m][(i-1)*vpts + jj] = temp2;
                            } */


                      } // end of jj

             } // end of iii

      } // end of jjj

      /* define weak plate boundaries from imposed surface velocities */
      //if(E->mesh.topvbc == 1) {
        // weak_plate_boundary_from_velocities(E,EEta);
      // }

      break;


    case 4:

        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(i=1;i<=nel;i++)   {
                l = E->mat[m][i];
		if(E->control.mat_control) /* moved this up here TWB */
		  tempa = E->viscosity.N0[l-1] * E->VIP[m][i];
		else
		  tempa = E->viscosity.N0[l-1];

                j = 0;

                for(kk=1;kk<=ends;kk++) {
                    TT[kk] = E->T[m][E->ien[m][i].node[kk]];
                    zz[kk] = (1.-E->sx[m][3][E->ien[m][i].node[kk]]);
                }

                for(jj=1;jj<=vpts;jj++) {
                    temp=0.0;
                    zzz=0.0;
                    for(kk=1;kk<=ends;kk++)   {
                        TT[kk]=max(TT[kk],zero);
                        temp += min(TT[kk],one) * E->N.vpt[GNVINDEX(kk,jj)];
                        zzz += zz[kk] * E->N.vpt[GNVINDEX(kk,jj)];
                    }

                    /* The viscosity formulation (dimensional) is: visc=visc0*exp[(Ea+p*Va)/R*T]
                       Typical values for dry upper mantle are: Ea = 300 KJ/mol ; Va = 1.e-5 m^3/mol
                       T=T0+DT*T'; where DT - temperature contrast (from Rayleigh number)
                       T' - nondimensional temperature; T0 - surface tempereture (273 K)
                       T=DT*[(T0/DT) + T'] => visc=visc0*exp{(Ea+p*Va)/R*DT*[(T0/DT) + T']}
                       visc=visc0*exp{[(Ea/R*DT) + (p*Va/R*DT)]/[(T0/DT) + T']}
                       so: E->viscosity.E = Ea/R*DT ; E->viscosity.Z = Va/R*DT
                       p = zzz and E->viscosity.T = T0/DT */

		    EEta[m][ (i-1)*vpts + jj ] = tempa*
		      exp( (E->viscosity.E[l-1] +  E->viscosity.Z[l-1]*zzz )
			   / (E->viscosity.T[l-1]+temp) );

                }
            }
        break;


    case 5:			/* this still needs to be documented, who wrote this? */

        /* same as rheol 3, except alternative margin, VIP, formulation */
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(i=1;i<=nel;i++)   {
                l = E->mat[m][i];

                tempa = E->viscosity.N0[l-1];
                j = 0;

                for(kk=1;kk<=ends;kk++) {
                    TT[kk] = E->T[m][E->ien[m][i].node[kk]];
                    zz[kk] = (1.-E->sx[m][3][E->ien[m][i].node[kk]]);
                }

                for(jj=1;jj<=vpts;jj++) {
                    temp=0.0;
                    zzz=0.0;
                    for(kk=1;kk<=ends;kk++)   {
                        TT[kk]=max(TT[kk],zero);
                        temp += min(TT[kk],one) * E->N.vpt[GNVINDEX(kk,jj)];
                        zzz += zz[kk] * E->N.vpt[GNVINDEX(kk,jj)];
                    }

                    if(E->control.mat_control==0)
                        EEta[m][ (i-1)*vpts + jj ] = tempa*
                            exp( E->viscosity.E[l-1]/(temp+E->viscosity.T[l-1])
                                 - E->viscosity.E[l-1]/(one +E->viscosity.T[l-1]) );

                    if(E->control.mat_control==1) {
                    /*    visc1 = E->VIP[m][i];
                        visc2 = 2.0/(1./visc1 + 1.);
                        tempa_exp = tempa*
	                exp( E->viscosity.E[l-1]/(temp+E->viscosity.T[l-1])
		           - E->viscosity.E[l-1]/(one +E->viscosity.T[l-1]) );
                        visc1 = tempa*E->viscosity.max_value;
                        if(tempa_exp > visc1) tempa_exp=visc1;
                        EEta[m][ (i-1)*vpts + jj ] = visc2*tempa_exp;
                     */
		     
                       visc1 = E->VIP[m][i];
                       visc2 = tempa*
	               exp( E->viscosity.E[l-1]/(temp+E->viscosity.T[l-1])
		          - E->viscosity.E[l-1]/(one +E->viscosity.T[l-1]) );
                       if(visc1 <= 0.95 && E->sx[m][3][E->ien[m][i].node[jj]]>0.97) visc2=visc1;
                       EEta[m][ (i-1)*vpts + jj ] = visc2;
		      

                      }
		     
 		     // Lijun's mod
		     fi_trench = 0.75;
		     dip = 15.0;
		     dip *=3.14159/180.0;
		     fi = E->sx[m][2][E->ien[m][i].node[jj]];
		     rad = E->sx[m][3][E->ien[m][i].node[jj]];
                     if(fabs((1.0-rad)-tan(dip)*(fi-fi_trench-0.015))<0.005 && rad>0.98 || fi < 0.02 && rad>0.99 || fi > 0.77 && rad>0.985 && rad<0.995 || fi > 0.95 && rad>0.99){
		    //if(fi > fi_trench && rad>0.99 || fi < 0.02 && rad>0.99){
			
                          EEta[m][ (i-1)*vpts + jj ] = 0.001;
			  //E->T[m][E->ien[m][i].node[jj]] = 1.0-0.5*(1.0-E->T[m][E->ien[m][i].node[jj]]);
			  //E->T[m][E->ien[m][i].node[jj]] = 1.0-0.5*(1.0-temp);
                          //fprintf(stderr,"visc=%f\n",EEta[m][ (i-1)*vpts + jj ]);
                     }
                     // End
		     

                }
            }
        break;


    case 6:			/* eta = N_0 exp(E(T_0-T) + (1-z) Z_0 ) */

        for(m=1;m <= E->sphere.caps_per_proc;m++)
	  for(i=1;i <= nel;i++)   {

	    l = E->mat[m][i] - 1;

	    if(E->control.mat_control)
	      tempa = E->viscosity.N0[l] * E->VIP[m][i];
	    else
	      tempa = E->viscosity.N0[l];
	    j = 0;

	    for(kk=1;kk<=ends;kk++) {
	      TT[kk] = E->T[m][E->ien[m][i].node[kk]];
	      zz[kk] = (1.0 - E->sx[m][3][E->ien[m][i].node[kk]]);
	    }

	    for(jj=1;jj <= vpts;jj++) {
	      temp=0.0;zzz=0.0;
	      for(kk=1;kk <= ends;kk++)   {
		TT[kk]=max(TT[kk],zero);
		temp += min(TT[kk],one) * E->N.vpt[GNVINDEX(kk,jj)];
		zzz += zz[kk] * E->N.vpt[GNVINDEX(kk,jj)];
	      }
	      EEta[m][ (i-1)*vpts + jj ] = tempa*
		exp( E->viscosity.E[l]*(E->viscosity.T[l] - temp) +
		     zzz *  E->viscosity.Z[l]);
	      //fprintf(stderr,"N0 %11g T %11g T0 %11g E %11g z %11g km Z %11g mat: %i log10(eta): %11g\n",
	      //      tempa,temp,E->viscosity.T[l],E->viscosity.E[l], zzz *6371 ,E->viscosity.Z[l],l+1,log10(EEta[m][ (i-1)*vpts + jj ]));
	    }
	  }
        break;




    }

    return;
}


/*
void visc_from_S(E,EEta,propogate)
     struct All_variables *E;
     float **EEta;
     int propogate;
{
    float one,two,scale,stress_magnitude,depth,exponent1;
    float *eedot;

    void strain_rate_2_inv();
    int m,e,l,z,jj,kk;

    const int vpts = vpoints[E->mesh.nsd];
    const int nel = E->lmesh.nel;

    eedot = (float *) malloc((2+nel)*sizeof(float));
    one = 1.0;
    two = 2.0;

    for(m=1;m<=E->sphere.caps_per_proc;m++)  {
      if(E->viscosity.sdepv_visited){
	
        // get second invariant for all elements 
        strain_rate_2_inv(E,m,eedot,1);
      }else{
	for(e=1;e<=nel;e++)	// initialize with unity if no velocities around 
	  eedot[e] = 1.0;
	E->viscosity.sdepv_visited = 1;

      }
        // eedot cannot be too small, or the viscosity will go to inf 
	for(e=1;e<=nel;e++){
	  eedot[e] = max(eedot[e], 1.0e-16);
	}

	// (Lijun put) Formula: eta = pow(eta,1/n)*pow(epcilon_dot,1/n-1)
        for(e=1;e<=nel;e++)   {
            exponent1= one/E->viscosity.sdepv_expt[E->mat[m][e]-1];
            scale=pow(eedot[e],exponent1-one);
            for(jj=1;jj<=vpts;jj++)
                EEta[m][(e-1)*vpts + jj] = scale*pow(EEta[m][(e-1)*vpts+jj],exponent1);
        }
    }

    free ((void *)eedot);
    return;
}
*/

void visc_from_S(E,EEta,propogate)
     struct All_variables *E;
     float **EEta;
     int propogate;
{
    float one,two,scale,stress_magnitude,depth,exponent1;
    float *eedot;
    double t1, t2;

    void strain_rate_2_inv();
    int m,e,l,z,jj,kk;

    const int vpts = vpoints[E->mesh.nsd];
    const int nel = E->lmesh.nel;

    eedot = (float *) malloc((2+nel)*sizeof(float));
    one = 1.0;
    two = 2.0;

    for(m=1;m<=E->sphere.caps_per_proc;m++)  {
      if(E->viscosity.sdepv_visited){

        /* get second invariant for all elements */
        strain_rate_2_inv(E,m,eedot,1);
      }else{
        for(e=1;e<=nel;e++)     /* initialize with unity if no velocities around */
          eedot[e] = 1.0;
        E->viscosity.sdepv_visited = 1;

      }
        /* eedot cannot be too small, or the viscosity will go to inf */
        for(e=1;e<=nel;e++){
          eedot[e] = max(eedot[e], 1.0e-16);
          // Lijun add t1
          t1=max(eedot[e],1e-16);
        }

        /* Lijun add normalization */
        MPI_Allreduce(&t1, &t2,1,MPI_DOUBLE,MPI_MAX,E->parallel.world);
        /*if(E->parallel.me == 0) {
           fprintf(E->fp,"peak strain rate @ proc %d = %f\n", E->parallel.me, t1);
           fprintf(E->fp,"global max strain rate = %f\n", t2);
        }*/
        if(t2 > 1.0) {
           for(e=1;e<=nel;e++)
             eedot[e] = eedot[e]/t2*10;
        }


        // (Lijun add Formula) eta = pow(eta,1/n)*pow(epcilon_dot,(1-n)/n)
        for(e=1;e<=nel;e++)   {
            exponent1= one/E->viscosity.sdepv_expt[E->mat[m][e]-1];
            //fprintf(E->fp,"exponet1=%f\n",exponent1);
            scale=pow(eedot[e],exponent1-one);
            for(jj=1;jj<=vpts;jj++)
                if(scale<1.0)
                    EEta[m][(e-1)*vpts + jj] = scale*pow(EEta[m][(e-1)*vpts+jj],exponent1);
        }
    }

    free ((void *)eedot);
    return;
}



void visc_from_P(E,EEta) /* "plasticity" implementation

			 viscosity will be limited by a yield stress

			 \sigma_y  = min(a + b * (1-r), y)

			 where a,b,y are parameters input via pdepv_a,b,y

			 and

			 \eta_y = \sigma_y / (2 \eps_II)

			 where \eps_II is the second invariant. Then

			 \eta_eff = (\eta_0 \eta_y)/(\eta_0 + \eta_y)

			 for pdepv_eff = 1

			 or

			 \eta_eff = min(\eta_0,\eta_y)

			 for pdepv_eff = 0

			 where \eta_0 is the regular viscosity


			 TWB

			 */
     struct All_variables *E;
     float **EEta;
{
    float *eedot,zz[9],zzz,tau,eta_p,eta_new;
    int m,e,l,z,jj,kk;

    const int vpts = vpoints[E->mesh.nsd];
    const int nel = E->lmesh.nel;
    const int ends = enodes[E->mesh.nsd];

    void strain_rate_2_inv();

    eedot = (float *) malloc((2+nel)*sizeof(float));

    for(m=1;m<=E->sphere.caps_per_proc;m++)  {

      if(E->viscosity.pdepv_visited){

        strain_rate_2_inv(E,m,eedot,1);	/* get second invariant for all elements */

      }else{
	for(e=1;e<=nel;e++)	/* initialize with unity if no velocities around */
	  eedot[e] = 1.0;
	if(m == E->sphere.caps_per_proc)
	  E->viscosity.pdepv_visited = 1;
	if((E->parallel.me == 0)&&(E->control.verbose)){
	  for(e=0;e < E->viscosity.num_mat;e++)
	    fprintf(stderr,"num mat: %i a: %g b: %g y: %g\n",
		    e,E->viscosity.pdepv_a[e],E->viscosity.pdepv_b[e],E->viscosity.pdepv_y[e]);
	}
      }

      for(e=1;e <= nel;e++)   {	/* loop through all elements */

	l = E->mat[m][e] -1 ;	/* material of this element */

	for(kk=1;kk <= ends;kk++) /* nodal depths */
	  zz[kk] = (1.0 - E->sx[m][3][E->ien[m][e].node[kk]]); /* for depth, zz = 1 - r */

	for(jj=1;jj <= vpts;jj++){ /* loop through integration points */

	  zzz = 0.0;		/* get mean depth of integration point */
	  for(kk=1;kk<=ends;kk++)
	    zzz += zz[kk] * E->N.vpt[GNVINDEX(kk,jj)];

	  /* depth dependent yield stress */
	  tau = E->viscosity.pdepv_a[l] + zzz * E->viscosity.pdepv_b[l];

	  /* min of depth dep. and constant yield stress */
	  tau = min(tau,  E->viscosity.pdepv_y[l]);

	  /* yield viscosity */
	  eta_p = tau/(2.0 * eedot[e] + 1e-7) + E->viscosity.pdepv_offset;


	  if(E->viscosity.pdepv_eff){
	    /* two dashpots in series */
	    eta_new  = 1.0/(1.0/EEta[m][ (e-1)*vpts + jj ] + 1.0/eta_p);
	  }else{
	    /* min viscosities*/
	    eta_new  = min(EEta[m][ (e-1)*vpts + jj ], eta_p);
	  }
	  //fprintf(stderr,"z: %11g mat: %i a: %11g b: %11g y: %11g ee: %11g tau: %11g eta_p: %11g eta_new: %11g eta_old: %11g\n",
	  //zzz,l,E->viscosity.pdepv_a[l], E->viscosity.pdepv_b[l],E->viscosity.pdepv_y[l],
	  //eedot[e],tau,eta_p,eta_new,EEta[m][(e-1)*vpts + jj]);
	  EEta[m][(e-1)*vpts + jj] = eta_new;
        } /* end integration point loop */
      }	/* end element loop */

    } /* end caps loop */
    free ((void *)eedot);
    return;
}


void visc_from_C_Liu( E, EEta)
     struct All_variables *E;
     float **EEta;
{
  double vmean,cc_loc[10],CC[10][9],cbackground;
  int m,i,jj,flavor,numtracers;
  int iempty = 0;

  const int vpts = vpoints[E->mesh.nsd];
  const int nel = E->lmesh.nel;
  const int ends = enodes[E->mesh.nsd];

  for(m=1;m<=E->sphere.caps_per_proc;m++)
     for(i=1;i<=nel;i++)   {
        numtracers = 0;

        for (flavor=0; flavor<E->trace.nflavors; flavor++)
            numtracers += E->trace.ntracer_flavor[m][flavor][i];

        if (numtracers == 0) {
            iempty++;
            /* fprintf(E->trace.fpt, "No tracer in element %d!\n", e); */
            continue;
        }

        for(jj=1;jj<=vpts;jj++) {
            if (E->control.tracer) {
               if (E->trace.ntracer_flavor[m][1][i] > 0.0) { // number of flavor 1 tracers per elem.
                   vmean = E->trace.ntracer_flavor[m][1][i] / (double)numtracers;
                   //EEta[m][ (i-1)*vpts + jj ] *= vmean*0.001;
                   EEta[m][ (i-1)*vpts + jj ] *= 0.001;
               }
            }
        } // end of jj
     } // end of i loop

     return;
}


/*

   multiple with c(mpositional factor which is determined by a geometric
mean average from the tracer composition, assuming two flavors and
compositions between zero and unity
*/
void visc_from_C( E, EEta)
     struct All_variables *E;
     float **EEta;
{
  double vmean,cc_loc[10],CC[10][9],cbackground,r;
  int m,l,z,jj,kk,i,p,q,nd;

  //void compute_elemental_composition_ratio();

  const int vpts = vpoints[E->mesh.nsd];
  const int nel = E->lmesh.nel;
  const int ends = enodes[E->mesh.nsd];


  //E->composition.ncomp = E->trace.nflavors;
  //fprintf(stderr,"before calling comp.\n");
  //fprintf(stderr,"E->composition.ncomp=%d\n",E->composition.ncomp);

  //compute_elemental_composition_ratio(E);

  /*if(E->parallel.me==0) {
    fprintf(E->fp,"icomp_rheo=%d\n",E->composition.icompositional_rheology);
    fprintf(E->fp,"Number of composition is %d\n",E->composition.ncomp);
    fprintf(E->fp,"Number of flavors is %d\n", E->trace.nflavors);
    for(m=1;m <= E->sphere.caps_per_proc;m++) 
      for(p=0; p<E->composition.ncomp; p++) {
	fprintf(E->fp,"E->viscosity.cdepv_ff[%d] is %f\n", p,E->viscosity.cdepv_ff[p]);
      }
  }
  fflush(E->fp);
  */  

  for(m=1;m <= E->sphere.caps_per_proc;m++)  {
    for(i = 1; i <= nel; i++){
      /* determine composition of each of the nodes of the
         element */
        for(p=0; p<E->composition.ncomp; p++) {
            for(kk = 1; kk <= ends; kk++){
                CC[p][kk] = E->composition.comp_node[m][p][E->ien[m][i].node[kk]];
                if(CC[p][kk] < 0) CC[p][kk]=0.0;
                if(CC[p][kk] > 1) CC[p][kk]=1.0;
            }
        }
    
        for(jj = 1; jj <= vpts; jj++) {
            /* concentration of background material */
            cbackground = 1;
            for(p=0; p<E->composition.ncomp; p++) {
                /* compute mean composition  */
                cc_loc[p] = 0.0;
                for(kk = 1; kk <= ends; kk++) {
                    cc_loc[p] += CC[p][kk] * E->N.vpt[GNVINDEX(kk, jj)];
                }
                cbackground -= cc_loc[p];
            }

            /* geometric mean of viscosity */
            vmean = cbackground * E->viscosity.cdepv_ff[0];
            for(p=0; p<E->composition.ncomp; p++) {
                vmean += cc_loc[p] * E->viscosity.cdepv_ff[p+1];
            }
            vmean = exp(vmean);

            /* multiply the viscosity with this prefactor */
            r = E->sx[m][3][E->ien[m][i].node[1]];
            if(r>0.98) EEta[m][ (i-1)*vpts + jj ] *= vmean;

        } /* end jj loop */
    } /* end el loop */
  } /* end cap */
  return;
}



void strain_rate_2_inv(E,m,EEDOT,SQRT)
     struct All_variables *E;
     float *EEDOT;
     int m,SQRT;
{
    void get_global_shape_fn();
    void velo_from_element();
    void construct_c3x3matrix_el();
    void get_ba_p();

    struct Shape_function GN;
    struct Shape_function_dA dOmega;
    struct Shape_function_dx GNx;

    double edot[4][4], rtf[4][9];
    double theta;
    double ba[9][9][4][7];
    float VV[4][9], Vxyz[7][9], dilation[9];

    int e, i, j, p, q, n;

    const int nel = E->lmesh.nel;
    const int dims = E->mesh.nsd;
    const int ends = enodes[dims];
    const int lev = E->mesh.levmax;
    const int ppts = ppoints[dims];
    const int sphere_key = 1;


    for(e=1; e<=nel; e++) {

        /* get shape function on presure gauss points */
        get_global_shape_fn(E, e, &GN, &GNx, &dOmega, 2,
                                sphere_key, rtf, lev, m);

        velo_from_element(E, VV, m, e, sphere_key);

        theta = rtf[1][1];


        /* Vxyz is the strain rate vector, whose relationship with
         * the strain rate tensor (e) is that:
         *    Vxyz[1] = e11
         *    Vxyz[2] = e22
         *    Vxyz[3] = e33
         *    Vxyz[4] = 2*e12
         *    Vxyz[5] = 2*e13
         *    Vxyz[6] = 2*e23
         * where 1 is theta, 2 is phi, and 3 is r
         */
        for(j=1; j<=ppts; j++) {
            Vxyz[1][j] = 0.0;
            Vxyz[2][j] = 0.0;
            Vxyz[3][j] = 0.0;
            Vxyz[4][j] = 0.0;
            Vxyz[5][j] = 0.0;
            Vxyz[6][j] = 0.0;
            dilation[j] = 0.0;
        }

        if ((theta < 0.09) || (theta > 3.05)) {
            /* When the element is close to the poles, use a more
             * precise method to compute the strain rate. */

            if ((e-1)%E->lmesh.elz==0) {
                construct_c3x3matrix_el(E,e,&E->element_Cc,&E->element_Ccx,lev,m,1);
            }

            get_ba_p(&(E->N), &GNx, &E->element_Cc, &E->element_Ccx,
                     rtf, E->mesh.nsd, ba);

            for(j=1;j<=ppts;j++)
                for(p=1;p<=6;p++)
                    for(i=1;i<=ends;i++)
                        for(q=1;q<=dims;q++) {
                            Vxyz[p][j] += ba[i][j][q][p] * VV[q][i];
                        }

        }
        else {
            for(j=1; j<=ppts; j++) {
                for(i=1; i<=ends; i++) {
                    Vxyz[1][j] += (VV[1][i] * GNx.ppt[GNPXINDEX(0, i, j)]
                                   + VV[3][i] * E->N.ppt[GNPINDEX(i, j)])
                        * rtf[3][j];
                    Vxyz[2][j] += ((VV[2][i] * GNx.ppt[GNPXINDEX(1, i, j)]
                                    + VV[1][i] * E->N.ppt[GNPINDEX(i, j)]
                                    * cos(rtf[1][j])) / sin(rtf[1][j])
                                   + VV[3][i] * E->N.ppt[GNPINDEX(i, j)])
                        * rtf[3][j];
                    Vxyz[3][j] += VV[3][i] * GNx.ppt[GNPXINDEX(2, i, j)];

                    Vxyz[4][j] += ((VV[1][i] * GNx.ppt[GNPXINDEX(1, i, j)]
                                    - VV[2][i] * E->N.ppt[GNPINDEX(i, j)]
                                    * cos(rtf[1][j])) / sin(rtf[1][j])
                                   + VV[2][i] * GNx.ppt[GNPXINDEX(0, i, j)])
                        * rtf[3][j];
                    Vxyz[5][j] += VV[1][i] * GNx.ppt[GNPXINDEX(2, i, j)]
                        + rtf[3][j] * (VV[3][i] * GNx.ppt[GNPXINDEX(0, i, j)]
                                       - VV[1][i] * E->N.ppt[GNPINDEX(i, j)]);
                    Vxyz[6][j] += VV[2][i] * GNx.ppt[GNPXINDEX(2, i, j)]
                        + rtf[3][j] * (VV[3][i]
                                       * GNx.ppt[GNPXINDEX(1, i, j)]
                                       / sin(rtf[1][j])
                                       - VV[2][i] * E->N.ppt[GNPINDEX(i, j)]);
                }
            }
        } /* end of else */

        if(E->control.inv_gruneisen != 0) {
            for(j=1; j<=ppts; j++)
                dilation[j] = (Vxyz[1][j] + Vxyz[2][j] + Vxyz[3][j]) / 3.0;
        }

        edot[1][1] = edot[2][2] = edot[3][3] = 0;
        edot[1][2] = edot[1][3] = edot[2][3] = 0;

        /* edot is 2 * (the deviatoric strain rate tensor) */
        for(j=1; j<=ppts; j++) {
            edot[1][1] += 2.0 * (Vxyz[1][j] - dilation[j]);
            edot[2][2] += 2.0 * (Vxyz[2][j] - dilation[j]);
            edot[3][3] += 2.0 * (Vxyz[3][j] - dilation[j]);
            edot[1][2] += Vxyz[4][j];
            edot[1][3] += Vxyz[5][j];
            edot[2][3] += Vxyz[6][j];
        }

        EEDOT[e] = edot[1][1] * edot[1][1]
            + edot[1][2] * edot[1][2] * 2.0
            + edot[2][2] * edot[2][2]
            + edot[2][3] * edot[2][3] * 2.0
            + edot[3][3] * edot[3][3]
            + edot[1][3] * edot[1][3] * 2.0;
    }

    if(SQRT)
	for(e=1;e<=nel;e++)
	    EEDOT[e] =  sqrt(0.5 *EEDOT[e]);
    else
	for(e=1;e<=nel;e++)
	    EEDOT[e] *=  0.5;

    return;
}


static void apply_low_visc_wedge_channel(struct All_variables *E, float **evisc)
{
    void parallel_process_termination();

    int i,j,m;
    const int vpts = vpoints[E->mesh.nsd];
    float *F;

    /* low viscosity channel/wedge require tracers to work */
    if(E->control.tracer == 0) {
        if(E->parallel.me == 0) {
            fprintf(stderr, "Error: low viscosity channel/wedge is turned on, "
                   "but tracer is off!\n");
            fprintf(E->fp, "Error: low viscosity channel/wedge is turned on, "
                   "but tracer is off!\n");
            fflush(E->fp);
        }
        parallel_process_termination();
    }


    F = (float *)malloc((E->lmesh.nel+1)*sizeof(float));
    for(i=1 ; i<=E->lmesh.nel ; i++)
        F[i] = 0.0;

    /* if low viscosity channel ... */
    if(E->viscosity.channel)
        low_viscosity_channel_factor(E, F);


    /* if low viscosity wedge ... */
    if(E->viscosity.wedge)
        low_viscosity_wedge_factor(E, F);


    for(i=1 ; i<=E->lmesh.nel ; i++) {
        if (F[i] != 0.0)
            for(m = 1 ; m <= E->sphere.caps_per_proc ; m++) {
                for(j=1;j<=vpts;j++) {
                    evisc[m][(i-1)*vpts + j] = F[i];
            }
        }
    }


    free(F);

    return;
}




static void low_viscosity_channel_factor(struct All_variables *E, float *F)
{
    int i, ii, k, m, e, ee;
    int nz_min[NCS], nz_max[NCS];
    const int flavor = 0;
    double rad_mean, rr;

    for(m=1; m<=E->sphere.caps_per_proc; m++) {
        /* find index of radius corresponding to lv_min_radius */
        for(e=1; e<=E->lmesh.elz; e++) {
            rad_mean = 0.5 * (E->sx[m][3][E->ien[m][e].node[1]] +
                              E->sx[m][3][E->ien[m][e].node[8]]);
            if(rad_mean >= E->viscosity.lv_min_radius) break;
        }
        nz_min[m] = e;

        /* find index of radius corresponding to lv_max_radius */
        for(e=E->lmesh.elz; e>=1; e--) {
            rad_mean = 0.5 * (E->sx[m][3][E->ien[m][e].node[1]] +
                              E->sx[m][3][E->ien[m][e].node[8]]);
            if(rad_mean <= E->viscosity.lv_max_radius) break;
        }
        nz_max[m] = e;
    }



    for(m=1; m<=E->sphere.caps_per_proc; m++) {
        for(k=1; k<=E->lmesh.elx*E->lmesh.ely; k++) {
            for(i=nz_min[m]; i<=nz_max[m]; i++) {
                e = (k-1)*E->lmesh.elz + i;

                rad_mean = 0.5 * (E->sx[m][3][E->ien[m][e].node[1]] +
                                  E->sx[m][3][E->ien[m][e].node[8]]);

                /* loop over elements below e */
                for(ii=i; ii>=nz_min[m]; ii--) {
                    ee = (k-1)*E->lmesh.elz + ii;

                    rr = 0.5 * (E->sx[m][3][E->ien[m][ee].node[1]] +
                                E->sx[m][3][E->ien[m][ee].node[8]]);

                    /* if ee has tracers in it and is within the channel */
                    if((E->trace.ntracer_flavor[m][flavor][ee] > 0) &&
                       (rad_mean <= rr + E->viscosity.lv_channel_thickness)) {
                           F[e] = E->viscosity.lv_reduction;
                           break;
                       }
                }
            }
        }
    }


    /** debug **
    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(e=1; e<=E->lmesh.nel; e++)
            fprintf(stderr, "lv_reduction: %d %e\n", e, F[e]);
    /**/

    return;
}


static void low_viscosity_wedge_factor(struct All_variables *E, float *F)
{
    int i, ii, k, m, e, ee;
    int nz_min[NCS], nz_max[NCS];
    const int flavor = 0;
    double rad_mean, rr;

    for(m=1; m<=E->sphere.caps_per_proc; m++) {
        /* find index of radius corresponding to lv_min_radius */
        for(e=1; e<=E->lmesh.elz; e++) {
            rad_mean = 0.5 * (E->sx[m][3][E->ien[m][e].node[1]] +
                              E->sx[m][3][E->ien[m][e].node[8]]);
            if(rad_mean >= E->viscosity.lv_min_radius) break;
        }
        nz_min[m] = e;

        /* find index of radius corresponding to lv_max_radius */
        for(e=E->lmesh.elz; e>=1; e--) {
            rad_mean = 0.5 * (E->sx[m][3][E->ien[m][e].node[1]] +
                              E->sx[m][3][E->ien[m][e].node[8]]);
            if(rad_mean <= E->viscosity.lv_max_radius) break;
        }
        nz_max[m] = e;
    }



    for(m=1; m<=E->sphere.caps_per_proc; m++) {
        for(k=1; k<=E->lmesh.elx*E->lmesh.ely; k++) {
            for(i=nz_min[m]; i<=nz_max[m]; i++) {
                e = (k-1)*E->lmesh.elz + i;

                rad_mean = 0.5 * (E->sx[m][3][E->ien[m][e].node[1]] +
                                  E->sx[m][3][E->ien[m][e].node[8]]);

                /* loop over elements below e */
                for(ii=i; ii>=nz_min[m]; ii--) {
                    ee = (k-1)*E->lmesh.elz + ii;

                    /* if ee has tracers in it */
                    if(E->trace.ntracer_flavor[m][flavor][ee] > 0) {
                        F[e] = E->viscosity.lv_reduction;
                        break;
                    }
                }
            }
        }
    }


    /** debug **
    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(e=1; e<=E->lmesh.nel; e++)
            fprintf(stderr, "lv_reduction: %d %e\n", e, F[e]);
    /**/

    return;
}


static void weak_plate_boundary_from_velocities(E,EEta)
    struct All_variables *E;
    float **EEta; 
{                 
    int m,i,j,k,l,p,q,t,jj,node,nodeg,el,bd_i,bd_j,bd_k,n11,n12,n21,n22;
    int nox,noy,nox1,noy1,noz1,lev,flag;
    float v_grd11,v_grd12,v_grd21,v_grd22,age;
    float theta,phi,r,tempa,temp,temp1,temp2;
    float find_age_in_MY(); 
                  
    const int vpts = vpoints[E->mesh.nsd];
    const int dims=E->mesh.nsd;
              
    nox=E->mesh.nox;
    noy=E->mesh.noy;
    nox1=E->lmesh.nox;
    noy1=E->lmesh.noy;
    noz1=E->lmesh.noz; 
    lev=E->mesh.levmax;
    age=find_age_in_MY(E);
    temp = 0.2;  
               
    for(m=1;m<=E->sphere.caps_per_proc;m++)
      for(k=1;k<=noy1;k++)
        for(j=1;j<=nox1;j++)    {
           nodeg = E->lmesh.nxs+j-1 + (E->lmesh.nys+k-2)*nox;
           temp1 = -1.0;
               
           for(p=1;p<=5;p++) {
              n11=nodeg-p; n12=nodeg+p; n21=nodeg-p*nox; n22=nodeg+p*nox;
              v_grd11 = fabs(E->velo_1[n11]-E->velo_1[n11+1]);
              v_grd12 = fabs(E->velo_1[n12]-E->velo_1[n12-1]);
              v_grd21 = fabs(E->velo_1[n21]-E->velo_1[n21+nox]);
              v_grd22 = fabs(E->velo_1[n22]-E->velo_1[n22-nox]);
              if(v_grd11>temp || v_grd12>temp || v_grd21>temp || v_grd22>temp) {
                temp1 = 0.01 + (10.0-0.01)*(1.0-cos(pow((p-1)/4.0,5)*3.14159/2));
              }
              if(temp1<0.0) {
                for(q=1;q<=p;q++) {
                  n11=nodeg-p-q*nox; n12=nodeg-p+q*nox;
                  n21=nodeg+p-q*nox; n22=nodeg+p+q*nox;
                  v_grd11 = fabs(E->velo_1[n11]-E->velo_1[n11+1+nox]);
                  v_grd12 = fabs(E->velo_1[n12]-E->velo_1[n12+1-nox]);
                  v_grd21 = fabs(E->velo_1[n21]-E->velo_1[n21-1+nox]);
                  v_grd22 = fabs(E->velo_1[n22]-E->velo_1[n22-1-nox]);
                  if(v_grd11>temp || v_grd12>temp || v_grd21>temp || v_grd22>temp)
                    temp1 = 0.01 + (10.0-0.01)*(1.0-cos(pow((q+p-1)/(4.0+q),5)*3.14159/2));
                  }
              }
              if(temp1<0.0) {
                for(t=1;t<=p;t++) {
                  n11=nodeg-t-p*nox; n12=nodeg-t+p*nox;
                  n21=nodeg+t-p*nox; n22=nodeg+t+p*nox;
                  v_grd11 = fabs(E->velo_1[n11]-E->velo_1[n11+1+nox]);
                  v_grd12 = fabs(E->velo_1[n12]-E->velo_1[n12+1-nox]);
                  v_grd21 = fabs(E->velo_1[n21]-E->velo_1[n21-1+nox]);
                  v_grd22 = fabs(E->velo_1[n22]-E->velo_1[n22-1-nox]);
                  if(v_grd11>temp || v_grd12>temp || v_grd21>temp || v_grd22>temp)
                    temp1 = 0.01 + (10.0-0.01)*(1.0-cos(pow((t+p-1)/(4.0+t),5)*3.14159/2));
                  }
              }
              if(temp1>0.0) p += 10;
           } //end of p

           if(temp1 > 0.0) {
             flag = 0; //0 mean not right on the subduction zone
             for(i=noz1;i>=1;i--)  {
               if(k == 1)
                 bd_k = 1;
               else
                 bd_k = 2;
               if(j == 1)
                 bd_j = 1;
               else
                 bd_j = 2;
               if(i == 1)
                 bd_i = 0;
               else
                 bd_i = 1;

               el = (i-bd_i) + (j-bd_j)*(noz1-1) + (k-bd_k)*(noz1-1)*(nox1-1);
               node = i + (j-1)*noz1 + (k-1)*noz1*nox1;

               theta = E->sx[m][1][node];
               phi = E->sx[m][2][node];
               r = E->sx[m][3][node];

               l = E->mat[m][el];
               tempa = E->viscosity.N0[l-1];

               if(r > 0.99) {
                    temp2 = temp1 + (E->viscosity.N0[1]-temp1)*(1.0-cos((1.0-r)/0.01*1.57));
                    for(jj=1;jj<=vpts;jj++) {
                        //temp2 = temp1 + (EEta[m][(el-1)*vpts+jj]-temp1)*(1.0-cos((1.0-r)/0.01*1.57));
                        if(i==noz1 && EEta[m][(el-1)*vpts+jj]<0.1)
                            flag = 1;
                        if(EEta[m][(el-1)*vpts+jj]>temp2 && flag==0) // && theta<=lat_max && theta>=lat_min)
                           //if(theta<lat_max+0.01 || r>0.99)
                                EEta[m][(el-1)*vpts+jj] = temp2;
                    }
               }
             } //end of i
           }
        } //end of j

  return;
}


