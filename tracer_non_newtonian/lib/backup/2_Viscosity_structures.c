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

    if(E->viscosity.CDEPV)	/* compositional prefactor */
      visc_from_C(E,evisc);

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
    float zzz,zz[9],lon_temp,r_temp,r_min,lon_max,guide_dep;
    float visc1, visc2, tempa_exp,*PB1[4],*PB2[4];
    float fi, theta, rad, r, fi_trench, the_trench, dip;
    float find_age_in_MY();
    float T_max, T_min,age,age1,age2;
    char input_s[1000],pb_1[255],pb_2[255],bounds[255];
    double Tmaxd();
    const int vpts = vpoints[E->mesh.nsd];
    const int ends = enodes[E->mesh.nsd];
    const int nel = E->lmesh.nel;
    const int dims=E->mesh.nsd;
    FILE *fp,*fp1,*fp2,*fp3,*fp4;

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

        age=find_age_in_MY(E);
        intage = age;
        age1=1.0*age;
        age2=age1+1.0;

        fp2=fopen("rheo.dat","r");
        if(fp2==NULL) rheo_trick = 0;
        fscanf(fp2,"%d",&(rheo_trick));
        fclose(fp2);

        if(rheo_trick >1 ) {
          if(rheo_trick>2){
            if(E->sphere.caps > 1 && rheo_trick>2) {
                cap = E->sphere.capid[1] - 1;
                sprintf(pb_1,"/home/lijun/jobs/BC_files/pb_file/pb%0.0f.%d.dat",age,cap);
                if(rheo_trick==3)
                        sprintf(pb_2,"../clr_pos/clr_%0.0f.xy",age);
            }
            else if(E->sphere.caps == 1)
                sprintf(pb_1,"/home/lijun/jobs/BC_files/pb_file/pb%0.0f.dat",age);

            if(E->parallel.me==0) {
                fprintf(E->fp,"Plate boundary: Current boundary = %s\n",pb_1);
                if(rheo_trick==3)
                        fprintf(E->fp,"C.P. boundary = %s\n",pb_2);
            }
          }
          else if(rheo_trick==2) {
            sprintf(pb_1,"../clr_pos/clr_%0.0f.xy",age);
          }

          if(rheo_trick>2){
            nnn=159;
            ttt=129;
          }
          else if(rheo_trick==2){
            nnn=33;
            ttt=nnn;
          }

          for(i=1;i<=dims;i++) {
              PB1[i]=(float*) malloc ((nnn)*sizeof(float));
              PB2[i]=(float*) malloc ((nnn)*sizeof(float));
          }

          fp1=fopen(pb_1,"r");
          if(rheo_trick==3)
                fp2=fopen(pb_2,"r");

          for(j=1;j<=nnn;j++)   {
                if(j<=ttt) {
                    if(rheo_trick>2) {
                        fscanf(fp1,"%f %f",&(PB1[1][j]),&(PB1[2][j]));
                        if(rheo_trick==3 && j<=33)
                               fscanf(fp2,"%f %f %f",&(PB2[1][j]),&(PB2[2][j]),&(PB2[3][j]));
                    }
                    else if(rheo_trick==2)
                        fscanf(fp1,"%f %f %f",&(PB1[1][j]),&(PB1[2][j]),&(PB1[3][j]));
                }
                else {
                    PB1[1][j] = PB1[1][1]-(j-ttt)*(PB1[1][2]-PB1[1][1]);
                    PB1[2][j] = PB1[2][1];
                }
                if(E->sphere.caps == 1)  // Regional model
                        if(fabs(PB1[1][j]-E->sx[1][1][1])<0.005) {
                                theta_start = j;
                        }
          }
          fclose(fp1);
          if (rheo_trick==3)
                fclose(fp2);
        }


	for(m=1;m<=E->sphere.caps_per_proc;m++)
          for(kkk=1;kkk<=E->lmesh.noy;kkk++)
            for(jjj=1;jjj<=E->lmesh.nox;jjj++){
                if(E->sphere.caps > 1 && rheo_trick != 0) {
                    flag = 0;
                    flag1 = 0;
                    temp1 = 0.1;
                    for(i=1;i<nnn;i++) {
                        node = 1 + (jjj-1)*E->lmesh.noz + (kkk-1)*E->lmesh.noz*E->lmesh.nox;
                        if(fabs(PB1[1][i]-E->sx[m][1][node])<0.01)
                            flag = i;
                        if(rheo_trick==3 && i<=33){
                            if(fabs(PB2[1][i]-E->sx[m][1][node])<=temp1){
                                flag1 = i;
                                temp1 = fabs(PB2[1][i]-E->sx[m][1][node]);
                            }
                        }
                    }
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

                //    theta = E->sx[m][1][node];
                //    fi = E->sx[m][2][node];
                    r = E->sx[m][3][node];

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
                            //temp += min(TT[kk],one) * E->N.vpt[GNVINDEX(kk,jj)];
                            temp += min(TT[kk],T_max) * E->N.vpt[GNVINDEX(kk,jj)];
                            zzz += zz[kk] * E->N.vpt[GNVINDEX(kk,jj)];
                        }

                        if(E->control.mat_control==0) {
                                //Lijun addition: maintain the same slab stiffness for all visc_um values
                                if(r > 0.91 && r < 1.0)
                                  if((E->control.TBCbotval-E->T[m][node]) > 0.8*(E->control.TBCbotval-T_min))
                                        tempa = 1.0*T_max/E->control.TBCbotval;
                                //End of Lijun's addition
                                EEta[m][ (i-1)*vpts + jj ] = tempa*exp( E->viscosity.E[l-1]/((temp-T_min)/(one-T_min)*one+E->viscosity.T[l-1]) - E->viscosity.E[l-1]/(one +E->viscosity.T[l-1]) );
                                if(E->T[m][node] < E->control.TBCbotval && r > 0.91)
                                        EEta[m][ (i-1)*vpts + jj ] *=1;
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
                                        if(rheo_trick == 4 || rheo_trick == 10) {
                                                lon_dcr=0.2;
                                                lon_inc1 = 1.5;
                                        }
                                        else if(rheo_trick >= 5) {
                                                lon_dcr=0.3;
                                                lon_inc = 0.5;
                                                if(rheo_trick == 6){
                                                    lon_inc=0.15;
                                                    temp1 = 100;
                                                }
                                                if(rheo_trick == 5){
                                                    lon_inc = 1.5;
                                                    temp1 = 100;
                                                    guide_dep = 0.97;
                                                }
                                                if(rheo_trick == 7){
                                                    lon_inc = 1.5;
                                                    lon_dcr = 0.3;
                                                }
                                                if(rheo_trick == 8){
                                                     lon_inc=0.05;
						     lon_dcr = 0.05;
                                                     guide_dep = 0.94;
                                                }
						if(rheo_trick == 9){
                                                     temp1 = 50;
                                                     lon_inc=0.15;
                                                     lon_dcr=0.1;
                                                     guide_dep = 0.89;
                                                }
                                                if(rheo_trick == 11){
                                                    lon_inc=1.2;
                                                    lon_inc1 = lon_inc;
                                                }
                                        }

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
                                    if(rheo_trick!=4){
                                        if(E->sx[m][2][node] > 0.0)
                                                fi = PB1[2][theta_start+jjj-1] + 2*3.1415926;
                                        else
                                                fi = E->sx[m][2][node];

                                        if((E->sx[m][2][node] > (fi-0.1)) && (r > 0.975))
                                                EEta[m][ (i-1)*vpts + jj ] = 100.0;
                                        if((E->sx[m][2][node] > (fi-0.05)) && (r > 0.98))
                                                EEta[m][ (i-1)*vpts + jj ] = 0.1;
                                    }
                                    else {
                                        temp1 = 500;
                                        if(E->sx[m][2][node]>fi-0.1 && E->sx[m][2][node
]<fi+lon_inc1 && r>guide_dep) {
                                            if(temp<0.5*E->control.TBCbotval && flag2<1*vpts) {
                                                EEta[m][ (i-1)*vpts + jj ] = temp1;
                                                flag2 +=1;
                                            }

					    // free subduction param.
                                       /*     if(E->sx[m][2][node]<fi+0.03 && r>0.97)
                                                EEta[m][ (i-1)*vpts + jj ] = temp1;
                                            if(fabs(E->sx[m][2][node]-fi-0.02)<0.02 && fabs(E->sx[m][3][node]-(0.99-0.92)/lon_inc1*(fi+lon_inc1-E->sx[m][2][node])-0.92)<=0.01) { // || E->sx[m][2][node]<fi && E->sx[m][3][node]>0.98) {
                                                EEta[m][ (i)*vpts + jj ] = temp2;
                                                EEta[m][ (i+1)*vpts + jj ] = temp2;
                                                EEta[m][ (i-1)*vpts + jj ] = temp2;
                                            }   
                                         */
                                            if (EEta[m][ (i+1)*vpts + jj ] == temp2 && EEta[m][ (i)*vpts + jj ] != temp2) {
                                                EEta[m][ (i)*vpts + jj ] = temp1;
                                                EEta[m][ (i-1)*vpts + jj ] = temp1;
                                                EEta[m][ (i+1)*vpts + jj ] = temp1;
                                            }                                         }
                                    }
                                }

                        }

                    }
                }
            }

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
	
        /* get second invariant for all elements */
        strain_rate_2_inv(E,m,eedot,1);
      }else{
	for(e=1;e<=nel;e++)	/* initialize with unity if no velocities around */
	  eedot[e] = 1.0;
	E->viscosity.sdepv_visited = 1;

      }
        /* eedot cannot be too small, or the viscosity will go to inf */
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

/*

multiply with compositional factor which is determined by a geometric
mean average from the tracer composition, assuming two flavors and
compositions between zero and unity

*/
void visc_from_C( E, EEta)
     struct All_variables *E;
     float **EEta;
{
  float comp,comp_fac,CC[9],tcomp;
  double vmean,cc_loc;
  int m,l,z,jj,kk,i;


  const int vpts = vpoints[E->mesh.nsd];
  const int nel = E->lmesh.nel;
  const int ends = enodes[E->mesh.nsd];
  if(E->trace.nflavors != 2)
    myerror(E,"sorry, CDEPV only supports two flavors");
  if(E->composition.ncomp != 1)
    myerror(E,"CDEPV only supports one composition yet");

  for(m=1;m <= E->sphere.caps_per_proc;m++)  {
    for(i = 1; i <= nel; i++){
      /* determine composition of each of the nodes of the
	 element */
      for(kk = 1; kk <= ends; kk++){
	CC[kk] = E->composition.comp_node[m][0][E->ien[m][i].node[kk]];
	if(CC[kk] < 0)CC[kk]=0.0;
	if(CC[kk] > 1)CC[kk]=1.0;
      }
      for(jj = 1; jj <= vpts; jj++){
	/* compute mean composition  */
	cc_loc = 0.0;
	for(kk = 1; kk <= ends; kk++)
	  cc_loc += CC[kk] * E->N.vpt[GNVINDEX(kk, jj)];

	/* geometric mean of viscosity */
	vmean = exp(cc_loc  * E->viscosity.cdepv_ff[1] +
		    (1.0-cc_loc) * E->viscosity.cdepv_ff[0]);
	/* multiply the viscosity with this prefactor */
	EEta[m][ (i-1)*vpts + jj ] *= vmean;
      } /* end jj loop */
    } /* end el loop */
  } /* end cap */
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
