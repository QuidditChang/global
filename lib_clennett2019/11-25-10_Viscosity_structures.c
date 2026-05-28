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

#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "parsing.h"


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
  }

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
    input_float("sdepv_misfit",&(E->viscosity.sdepv_misfit),"0.001",m);
    input_float_vector("sdepv_expt",E->viscosity.num_mat,(E->viscosity.sdepv_expt),m);
  }

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

  input_string("Viscosity",E->viscosity.STRUCTURE,NULL,m);
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
    void apply_viscosity_smoother();
    void visc_to_node_interpolate();
    void visc_from_nodes_to_gint();
    void visc_from_gint_to_nodes();


    int i,j,m;
    float temp1,temp2,*vvvis;
    double *TG;

    const int vpts = vpoints[E->mesh.nsd];

    if(E->viscosity.TDEPV)
       visc_from_T(E,evisc,propogate);
    else
       visc_from_mat(E,evisc);

    if(E->viscosity.SDEPV)
       visc_from_S(E,evisc,propogate);

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

    /*
 if (E->control.verbose)  {
    fprintf(E->fp_out,"output_evisc \n");
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
      fprintf(E->fp_out,"output_evisc for cap %d\n",E->sphere.capid[m]);
      for(i=1;i<=E->lmesh.nel;i++)
        fprintf(E->fp_out,"%d %d %f %f\n",i,E->mat[m][i],evisc[m][(i-1)*vpts+1],evisc[m][(i-1)*vpts+7]);
      }
    fflush(E->fp_out);
    }
    */

    visc_from_gint_to_nodes(E,evisc,visc,E->mesh.levmax);

    visc_from_nodes_to_gint(E,visc,evisc,E->mesh.levmax);

/*    visc_to_node_interpolate(E,evisc,visc);
*/

/*    for(m=1;m<=E->sphere.caps_per_proc;m++) {
      for(i=1;i<=E->lmesh.nel;i++)
	if (i%E->lmesh.elz==0) {
          fprintf(E->fp_out,"%.4e %.4e %.4e %5d %2d\n",E->eco[m][i].centre[1],E->eco[m][i].centre[2],log10(evisc[m][(i-1)*vpts+1]),i,E->mat[m][i]);

	  }
        }  */
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
    int m,i,j,k,l,jj,kk,iii,jjj,kkk,imark,intage,node,nnn,ttt,count,bd_i,bd_j,bd_k,nd_core,el_core,lv,ndlv,el,ie;
    int theta_start,fi_start,rheo_trick,cap,flag,flag1,flag2,*PB3[4],exponent1;
    int nodeg,node_coarse,node_fine,start_lev,start_node;
    float lat_max,lat_min,lat_max1,lat_min1,lat_max2,lat_min2,lon_inc,lon_inc1,lon_dcr;
    float find_age_in_MY();
    float zero,e_6,one,eta0,Tave,depth,temp,tempa,temp1,temp2,temp_core,TT[9],viscE,viscT;
    float zzz,zz[9],theta,phi,r,fi_1,fi_2,fi,tt1,tt2,ttt1,ttt2,*PB1[4],*PB2[4];
    float visc1, visc2, tempa_exp, *eedot, scale,x,y,z,x0,y0,z0,x1,y1,z1,distance,dist1;
    double t1, t2;
    const int vpts = vpoints[E->mesh.nsd];
    const int ends = enodes[E->mesh.nsd];
    const int nel = E->lmesh.nel;
    float T_max, T_min,age,age1,age2;
    char input_s[1000],pb_1[255],pb_2[255],bounds[255];
    double Tmaxd();
    FILE *fp,*fp1,*fp2,*fp3,*fp4;
    const int dims=E->mesh.nsd;

    void strain_rate_2_inv();
    float theta_locater();

    e_6 = 1.e-6;
    one = E->control.lith_age_mantle_temp;
    //kkk = (E->lmesh.noy-1)*(E->lmesh.noz-1)*(E->lmesh.nox-1);
    //fprintf(stderr,"nel=%d, nno=%d\n",nel,kkk);
    zero = 0.0;
    imark = 0;

    switch (E->viscosity.RHEOL)   {
    case 1:
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

    case 2:
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

    case 3:

      fp=fopen("filter_limit","r");
      if(!fp) {
	      T_max=Tmaxd(E,E->T);
	      T_min=Tmind(E,E->T);
	      T_max *= 1.1; T_min *= 0.9;
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

      age=find_age_in_MY(E);
      intage = age;
      age1=1.0*intage;
      age2=age1+1.0;

      fp2=fopen("rheo.dat","r");
      if(fp2==NULL) 
	rheo_trick = 0;
      else {
	fscanf(fp2,"%d",&(rheo_trick));
      	fclose(fp2);
      }
      fprintf(E->fp,"Rheo_trick=%d\n",rheo_trick);

      if(rheo_trick >1 ) {
	if(rheo_trick>2){
	    if(E->sphere.caps > 1) {
		//cap = E->sphere.capid[1] - 1; 
		//sprintf(pb_1,"/home/lijun/jobs/BC_files/pb_file/pb%0.0f.%d.dat",age,cap);
		sprintf(pb_1,"/share/home/01523/lijun/work/jobs/BC_files/Nam_subduct/subduction_sL.%0.0f.xy",age1);
		sprintf(pb_2,"/share/home/01523/lijun/work/jobs/BC_files/Nam_subduct/subduction_sL.%0.0f.xy",age2);
	    }
	    else if(E->sphere.caps == 1) {
		//sprintf(pb_1,"/home/lijun/jobs/BC_files/pb_file/pb%0.0f.dat",age);
		sprintf(pb_1,"/home/01523/lijun/work/BC_files/Nam_subduct/subduction_sL.%0.0f.xy",age1);
		sprintf(pb_2,"/home/01523/lijun/work/BC_files/Nam_subduct/subduction_sL.%0.0f.xy",age2);
	    }
            if(E->parallel.me==0)
		fprintf(E->fp,"Plate boundary: %s,\n %s\n",pb_1,pb_2);
	}

	/* define the maximum array length */
	if(rheo_trick>2){
	    nnn=150;
	    ttt=150;
	}
        for(i=1;i<=dims;i++) {
              PB1[i]=(float*) malloc ((nnn+1)*sizeof(float));
	      PB2[i]=(float*) malloc ((nnn+1)*sizeof(float));
	}
	/* read in the plate boundary locations */
	fp1=fopen(pb_1,"r");	
	fp2=fopen(pb_2,"r");
        for(j=1;j<=nnn;j++)   {
		if(j<=ttt) {
		    if(rheo_trick>2) {
			fscanf(fp1,"%f %f",&(tt1),&(tt2));
			fscanf(fp2,"%f %f",&(ttt1),&(ttt2));
			if(tt1) {
                    	    PB1[1][j] = tt1;
                    	    PB1[2][j] = tt2;
                    	    tt1 = -10000.0;
                	}   
		    }
                    if(ttt1) {
                    	PB2[1][j] = ttt1;
                    	PB2[2][j] = ttt2;
                    	ttt1 = -10000.0;
                    }
		} //end of if(j)
	} //end of for (j)
        fclose(fp1); 
	fclose(fp2);

        // find the max/min latitude of plate boundary
	lat_max1 = 0.0;
        lat_min1 = 3.14159;
        for(i=1;i<nnn;i++) {
             if(fabs(PB1[1][i])<4.0 && PB1[1][i]>=lat_max1)
                   lat_max1 = PB1[1][i];
             if(fabs(PB1[1][i])<4.0 && PB1[1][i]<=lat_min1)
                   lat_min1 = PB1[1][i];
        }
        fprintf(E->fp,"lat_max1=%f, lat_min1=%f\n",lat_max1,lat_min1);
        lat_max2 = 0.0;
        lat_min2 = 3.14159;
        for(i=1;i<nnn;i++) {
             if(fabs(PB2[1][i])<4.0 && PB2[1][i]>=lat_max2)
                   lat_max2 = PB2[1][i];
             if(fabs(PB2[1][i])<4.0 && PB2[1][i]<=lat_min2)
                   lat_min2 = PB2[1][i];
        }
        fprintf(E->fp,"lat_max2=%f, lat_min2=%f\n",lat_max2,lat_min2);
        lat_max = lat_max1 + (lat_max2-lat_max1)*(age-age1)/(age2-age1);
        lat_min = lat_min1 + (lat_min2-lat_min1)*(age-age1)/(age2-age1);
        fprintf(E->fp,"lat_max=%f, lat_min=%f\n",lat_max,lat_min);
      }

      /* loop through all nodes & apply rheo_trick */ 
      for(m=1;m<=E->sphere.caps_per_proc;m++)
	for(kkk=1;kkk<=E->lmesh.noy;kkk++) 
	    for(jjj=1;jjj<=E->lmesh.nox;jjj++){
		if(E->sphere.caps >= 1 && rheo_trick != 0) {
		    //find two nearest plate boundary  lat's to a given mesh point 
		    node = 1 + (jjj-1)*E->lmesh.noz + (kkk-1)*E->lmesh.noz*E->lmesh.nox;
		    fi_1 = fi_2 = 0.0;
		    fi_1 = theta_locater(PB1,nnn,E->sx[m][1][node]);
                    fi_2 = theta_locater(PB2,nnn,E->sx[m][1][node]);
                    fi = fi_1 + (fi_2-fi_1)*(age-age1)/(age2-age1);
		    if(fi != 0.0) flag = 1;
		}

		nodeg=E->lmesh.nxs-1+jjj+(E->lmesh.nys+kkk-2)*E->mesh.nox;
		nd_core = -1;
		el_core = -1;
		temp_core = 0.8*E->control.lith_age_mantle_temp; // temperature at the slab center
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
                    
		    theta = E->sx[m][1][node];
		    phi = E->sx[m][2][node];
		    r = E->sx[m][3][node];

	  	    l = E->mat[m][i];
	  	    tempa = E->viscosity.N0[l-1];
		    /* creat a smooth radial visc. profile */
		    if(r>0.993){
                        tempa = E->viscosity.N0[0] + (E->viscosity.N0[1]-E->viscosity.N0[0])*(cos((r-0.993)/0.007*3.14159)+1)/2;
			viscE = E->viscosity.E[0] + (E->viscosity.E[1]-E->viscosity.E[0])*(cos((r-0.993)/0.007*3.14159)+1)/2;
			viscT = E->viscosity.T[0] + (E->viscosity.T[1]-E->viscosity.T[0])*(cos((r-0.993)/0.007*3.14159)+1)/2;
		    }
                    else if(r<0.98 && r>0.92) {
                        tempa = E->viscosity.N0[1] + (E->viscosity.N0[2]-E->viscosity.N0[1])*(cos((r-0.92)/0.06*3.14159)+1)/2;
			viscE = E->viscosity.E[1] + (E->viscosity.E[2]-E->viscosity.E[1])*(cos((r-0.92)/0.06*3.14159)+1)/2;
                        viscT = E->viscosity.T[1] + (E->viscosity.T[2]-E->viscosity.T[1])*(cos((r-0.92)/0.06*3.14159)+1)/2;
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


		    if(E->parallel.me_loc[3]==E->parallel.nprocz-1 && rheo_trick>=6) {
			if(E->age_t[nodeg]<5.0/E->data.scalet && E->sphere.cap[m].VB[2][node]*E->data.timedir>0.0 && theta<=lat_max+0.015 && theta>=lat_min) 
                            E->sphere.cap[m].VB[2][node] *= (cos(fabs(E->age_t[nodeg]*E->data.scalet)/5.0*3.14159)+1)/2.0+1;
			//if(phi<fi & phi>fi-0.05) {
			if(fabs(phi-fi)<0.03) {
			   if(phi<fi)
			       E->sphere.cap[m].VB[2][node] *= cos(fabs(phi-fi)/0.03*1.57)+1;
			   //else
			   //    E->sphere.cap[m].VB[3][node] = E->data.timedir*(-1000)*(cos(fabs(phi-fi)/0.03*1.57)+0.0);
			}
		    } 

		    /*if(E->parallel.me_loc[3]==E->parallel.nprocz-1 && rheo_trick>=6) {
			if(E->sphere.cap[m].VB[2][node]*E->data.timedir>0.0 && theta<=lat_max+0.015 && theta>=lat_min) {
			    flag1 = 0;
			    for(lv=E->mesh.gridmax;lv>=E->mesh.gridmin;lv--) {
				if(lv==E->mesh.gridmax) {
				  E->NODE[lv][m][node] = (INTX | INTY | INTZ);

				  E->NODE[lv][m][node] = E->NODE[lv][m][node] & (~VBX);
				  E->NODE[lv][m][node] = E->NODE[lv][m][node] | VBZ;
				  E->sphere.cap[m].VB[3][node] = 0.0;
				  E->NODE[lv][m][node] = E->NODE[lv][m][node] & (~VBY);
				  E->NODE[lv][m][node] = E->NODE[lv][m][node] | SBX;
				  E->sphere.cap[m].VB[1][node] = 0.0;
				  E->NODE[lv][m][node] = E->NODE[lv][m][node] & (~SBZ);
				  E->NODE[lv][m][node] = E->NODE[lv][m][node] | SBY;
				  E->sphere.cap[m].VB[2][node] = 0.0;
				  start_lev = lv;
				  start_node = node;
				  flag1 = 1;
				  fprintf(stderr,"level=%d, node=%d\n",lv,node);
				}
				else if (flag1==1){
				  flag1 = 0;
			    	  for(el=1;el<=E->lmesh.NEL[lv];el++)
			    	    for(ie=1;ie<=enodes[dims];ie++)       {
			    	      node_coarse = E->IEN[lv][m][el].node[ie];
			    	      node_fine=E->IEN[start_lev][m][E->EL[lv][m][el].sub[ie]].node[ie];
				      if(node_fine==start_node && flag1==0) {
					E->NODE[lv][m][node_coarse] = (INTX | INTY | INTZ);
					fprintf(stderr,"level=%d,start_node=%d,end_node=%d\n",lv,start_node,node_coarse);
					E->NODE[lv][m][node_coarse] = E->NODE[lv][m][node_coarse] & (~VBX);
					E->NODE[lv][m][node_coarse] = E->NODE[lv][m][node_coarse] | VBZ;
					E->NODE[lv][m][node_coarse] = E->NODE[lv][m][node_coarse] & (~VBY);
					E->NODE[lv][m][node_coarse] = E->NODE[lv][m][node_coarse] | SBX;
					E->NODE[lv][m][node_coarse] = E->NODE[lv][m][node_coarse] & (~SBZ);
					E->NODE[lv][m][node_coarse] = E->NODE[lv][m][node_coarse] | SBY;
					start_lev = lv;
					start_node = node_coarse;
					flag1 = 1;
					//break;
			    	      }
				    } //end of for(ie)
				} //end of else
			    } //end of for (lv)
			}
                    }  */

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
				//EEta[m][ (i-1)*vpts + jj ] = tempa*exp( E->viscosity.E[l-1]/(temp+E->viscosity.T[l-1]) - E->viscosity.E[l-1]/(one +E->viscosity.T[l-1]) );
				EEta[m][ (i-1)*vpts + jj ] = tempa*exp(viscE/(temp+viscT) - viscE/(one+viscT) );
				if(EEta[m][(i-1)*vpts+jj]/tempa>50 && fabs(phi-fi)<0.2) 
				    EEta[m][(i-1)*vpts+jj]=50*tempa;
			}	
	    	        if(E->control.mat_control==1) {
				//Lijun addition: maintain the same slab stiffness for all visc_um values
				if(r > 0.91 && r < 1.0)
				  if((E->control.lith_age_mantle_temp-E->T[m][node]) > 0.8*(E->control.lith_age_mantle_temp-T_min))
					tempa = 1.0*T_max/E->control.lith_age_mantle_temp;
                                //End of Lijun's addition
		      	    	EEta[m][ (i-1)*vpts + jj ] = tempa*E->VIP[m][i]*exp( E->viscosity.E[l-1]/((temp-T_min)/(one-T_min)*one+E->viscosity.T[l-1]) - E->viscosity.E[l-1]/(one +E->viscosity.T[l-1]) );
			}
			

			else if(rheo_trick == 2) { //lateral viscosity variation underneath continent
				temp1 = 200.0;
                                temp2 = 0.1;
				fi = PB1[2][flag]+2*3.1415926;
				if(r>0.97 && E->sx[m][2][node]>=fi && E->sx[m][2][node]<=fi+PB1[3][flag] && age >50) {
					EEta[m][ (i)*vpts + jj ] = temp1;
				}
			}
			else if(rheo_trick >= 3) { // stress guide + lateral viscosity variation
				temp1 = 100.0;
				temp2 = 0.1;
				if(flag > 0 && E->sphere.caps >= 1) {
					lon_inc = 1.5;
					lon_inc1 = lon_inc;
					lon_dcr = 0.15;
					if(rheo_trick==6) lon_dcr = 0.01;
                                        if(rheo_trick==7) lon_dcr = 0.03;
                                        if(rheo_trick>=8) lon_dcr = 0.1;


					if(rheo_trick==5 && r>0.975 && theta<=lat_max+0.04 && theta>=lat_min){
                                            if (phi>fi-0.02 && phi<fi+lon_dcr) {
                                                x0 = 0.97*sin(theta)*cos(fi+0.0315-(lat_max-theta)*0.04/(lat_max-lat_min));
                                                y0 = 0.97*sin(theta)*sin(fi+0.0315-(lat_max-theta)*0.04/(lat_max-lat_min));
                                                z0 = 0.97*cos(theta);
                                                x1 = 0.97*sin(theta)*cos(fi-0.0475-(lat_max-theta)*0.04/(lat_max-lat_min));
                                                y1 = 0.97*sin(theta)*sin(fi-0.0475-(lat_max-theta)*0.04/(lat_max-lat_min));
                                                z1 = 0.97*cos(theta);
                                                x = r*sin(theta)*cos(phi);
                                                y = r*sin(theta)*sin(phi);
                                                z = r*cos(theta);
                                                distance = sqrt(pow(x-x0,2)+pow(y-y0,2)+pow(z-z0,2));
                                                dist1 = sqrt(pow(x-x1,2)+pow(y-y1,2)+pow(z-z1,2));
                                                //temp1 = fabs(distance-0.08);
                                                temp1 = fabs((distance+dist1)-0.04055*2);
                                                temp2 = 0.003 + (1.0-r)/0.025*0.015;
                                                lat_max1 = lat_max+0.04*(1.0-r)/(1.0-0.975);
                                                if(temp1<temp2) {
                                                  temp2 = (1.0 + cos(temp1/temp2*3.14))/2; // *cos((r-0.97)/0.03*1.57);
                                                  if(theta>=lat_max1-0.005)
                                                    temp2 *= cos((theta-lat_max1+0.005)*1.5708/0.005);
                                                  if(theta<=lat_min+0.005)
                                                    temp2 *= cos((theta-lat_min-0.005)*1.5708/0.005);
                                                  temp1 = cos((r-0.975)/0.025*1.57); //weighting function (radial)
                                                  temp2 = E->control.lith_age_mantle_temp*(1.0-temp2/2); //defined Temp.
                                                  if(theta<=lat_max1)
                                                    E->T[m][node] = temp2*(1.0-temp1)+E->T[m][node]*temp1;
                                                  else if(E->T[m][node]<E->control.lith_age_mantle_temp)
                                                    E->T[m][node] = E->control.lith_age_mantle_temp;
                                                }
                                                else if(temp1<0.04 && r<0.995 && E->T[m][node]<E->control.lith_age_mantle_temp)
                                                  E->T[m][node] = E->control.lith_age_mantle_temp+(E->T[m][node]-E->control.lith_age_mantle_temp)*(temp1-temp2)/(0.04-temp2);
						if(r>0.993 && fabs(phi-fi)<0.06) {
						  temp1 = E->control.lith_age_mantle_temp + (0.0-E->control.lith_age_mantle_temp)*(E->sx[m][3][node]-0.993)/(E->sphere.ro-0.993);
						  //if(fabs(E->T[m][node]-E->control.lith_age_mantle_temp)<0.05)
						  if(fabs(E->T[m][node]>temp1))
							E->T[m][node] = temp1;
						}
                                            }
					    /*if(r>0.99 && fabs(phi-fi)<0.06) {
					    	temp1 = E->control.lith_age_mantle_temp + (0.0-E->control.lith_age_mantle_temp)*(E->sx[m][3][node]-0.99)/(E->sphere.ro-0.99);
                    			    	if(fabs(E->T[m][node]-E->control.lith_age_mantle_temp)<0.05)
                        			    e->T[m][node] = temp1;
					    }*/
                                            /*if (r>0.995 && fabs(phi-fi)<0.05) {
                                                temp1 = E->control.lith_age_mantle_temp*sin((r-0.995)/0.005*1.57);
                                                if(E->T[m][node]<temp1) E->T[m][node]=temp1;
                                            }*/
                                        }
                                        //end of rheo_trick = 5

					if(rheo_trick>=6 && r>0.97 && theta<=lat_max+0.015 && theta>=lat_min+0.02) {
                                            /* add a weak slab top */
					    temp1 = phi-(fi+3*(1.0-r));
                                            if(temp1>-0.08 && temp1<0.0 && r>0.98) { // && r<0.999){
                                                temp2 = 1.0 - cos((r-0.98)/0.02*1.57); //smooth in depth
                                                temp2 = 0.01+(1.0-(temp2*cos((pow(temp1/0.08,3)*1.57)))); //in fi
						if(EEta[m][(i-1)*vpts+jj]>temp2) EEta[m][(i-1)*vpts+jj]=temp2;
						   EEta[m][(i-1)*vpts+jj]=EEta[m][(i-1)*vpts+jj];
                                            }
					    /* add a weak subduction zone */
					    temp1 = phi-(fi+4*(1.0-r));
                                            if(temp1>-0.02 && temp1<0.04 && r>0.985){
						if(temp1>-0.02 && temp1<0.0) temp2 = (cos(pow(temp1/0.02,2)*3.1416)+1)/2;
						if(temp1>=0.0 && temp1<0.04)  temp2 = (cos(pow(temp1/0.04,1)*3.1416)+1)/2;
                                                temp2 = 0.01 + (1.0 - temp2);
                                                if(EEta[m][(i-1)*vpts+jj]>temp2) EEta[m][(i-1)*vpts+jj]=temp2;
						if(r<0.99 && EEta[m][(i-1)*vpts+jj]<E->viscosity.N0[1])
                                                        EEta[m][(i-1)*vpts+jj] = E->viscosity.N0[1];
                                            }
					    if(temp1<0.0 && E->T[m][node]<temp_core) {
						nd_core = node; 
						el_core = i;
						temp_core = E->T[m][node];
					    }
					    /* add a weak MOR */
					    if(E->age_t[nodeg]<5.0/E->data.scalet && r>0.98) { // && phi>fi-0.14) {
						temp1 = 1.0-(cos(pow(E->age_t[nodeg]/(5.0/E->data.scalet),2)*3.14159)+1)/2;
						EEta[m][(i-1)*vpts+jj] *= temp1;
						if(EEta[m][(i-1)*vpts+jj] < 0.02)
						   EEta[m][(i-1)*vpts+jj] = 0.02;
						if(r<0.99) {
						  temp1 = EEta[m][(i-1)*vpts+jj];
						  temp2 = temp1 + (E->viscosity.N0[1]-temp1)*(0.99-r)/(0.99-0.98);
                                                  EEta[m][(i-1)*vpts+jj] = temp2;
						}
					    }
					    /* add a weak transform fault */
					    if(fabs(theta-lat_max)<0.015 && r>0.98 && phi<fi && phi>fi-0.12) {
						temp1 = 1.0-(cos((theta-lat_max)/0.015*3.14159)+1)/2;
                                                temp2 = temp1*E->viscosity.N0[0];
						if(EEta[m][(i-1)*vpts+jj]>temp1) EEta[m][(i-1)*vpts+jj]=temp1;
                                                if(EEta[m][(i-1)*vpts+jj]<E->viscosity.N0[1])
                                                        EEta[m][(i-1)*vpts+jj] = E->viscosity.N0[1];
					    }
                                            /* add a weak slab hinge */
                                            temp1 = fabs(phi-(fi-(1.0-r)));
                                            if(temp1<0.03 && r>0.98) { //add a 10% weaker hinge
                                                temp2 = 0.1 + (1.0-cos(pow(temp1/0.03,0.5)*1.57));
                                                if(temp2<1.0) {
                                                    EEta[m][(i-1)*vpts+jj] *= temp2;
                                                    if(r<0.99 && EEta[m][(i-1)*vpts+jj]<E->viscosity.N0[1])
                                                        EEta[m][(i-1)*vpts+jj] = E->viscosity.N0[1];
                                                }
                                            }
					    /* add a weak layer atop slab */
					    temp1 = phi-fi;
					    if(rheo_trick==8 && flag2<=2*vpts && temp1>-0.1 && temp1<0.0)
                                                if(E->T[m][node]<=0.95*E->control.lith_age_mantle_temp && r>0.99) {
                                                    //temp2 = 1.0 - cos((flag2/vpts-2)/2.0*1.571);
                                                    EEta[m][(i-1)*vpts+jj] *= 0.001;
                                                    EEta[m][(i-2)*vpts+jj] *= 0.01;
                                                    EEta[m][(i-3)*vpts+jj] *= 0.1;
                                                    //fprintf(stderr,"iii=%d,flag2/vpts=%d\n",iii,flag2/vpts);
                                                    if(iii<=E->lmesh.noz) EEta[m][(i)*vpts+jj] *= 0.001;
                                                    if(iii<=E->lmesh.noz-1) EEta[m][(i+1)*vpts+jj] *= 0.01;
                                                    if(iii<=E->lmesh.noz-2) EEta[m][(i+2)*vpts+jj] *= 0.1;
                                                    flag2 += 1;
                                                }
					} // end of rheo_trick>=6

                                } //end of caps > 1

			} //end of rheo_trick>=3
			/* add a weak boundary surrounding the box */
			temp1 = E->viscosity.N0[l-1];
                        if(theta > E->control.theta_max-0.2) 
			    if(theta > E->control.theta_max-0.15) {
				temp2 = temp1 + (0.01-temp1)*(cos(pow((theta-(E->control.theta_max-0.15))/0.15,10)*3.14159)+1)/2;
                                EEta[m][(i-1)*vpts + jj] = temp2;
			    }
			    else {
			        temp2 = temp1 + (0.01-temp1)*(cos(pow((theta-(E->control.theta_max-0.15))/0.05,3)*3.14159)+1)/2;
			        if(EEta[m][(i-1)*vpts + jj] > temp2) 
				  EEta[m][(i-1)*vpts + jj] = temp2;
                            }
                        if(theta < E->control.theta_min+0.2) 
			    if(theta < E->control.theta_min+0.15) {
				temp2 = temp1+(0.01-temp1)*(cos(pow((E->control.theta_min+0.15-theta)/0.15,10)*3.14159)+1)/2;
                                EEta[m][(i-1)*vpts + jj] = temp2;
			    }
			    else {
                                temp2 = temp1+(0.01-temp1)*(cos(pow((E->control.theta_min+0.15-theta)/0.05,3)*3.14159)+1)/2;
			        if(EEta[m][(i-1)*vpts + jj] > temp2)
                                  EEta[m][(i-1)*vpts + jj] = temp2;
                            } 
			if(phi > E->control.fi_max-0.4)
			    if(phi > E->control.fi_max-0.3) {
				temp2 = temp1+(0.01-temp1)*(cos(pow((phi-(E->control.fi_max-0.3))/0.3,10)*3.14159)+1)/2;
                                if(theta <= E->control.theta_max-0.03 && theta >= E->control.theta_min+0.03)
                                  EEta[m][(i-1)*vpts + jj] = temp2; 
			    }
			    else if(theta <= E->control.theta_max-0.05 && theta >= E->control.theta_min+0.05){
                                temp2 = temp1+(0.01-temp1)*(cos(pow((phi-(E->control.fi_max-0.3))/0.1,3)*3.14159)+1)/2;
			        if(EEta[m][(i-1)*vpts + jj] > temp2)
                                  EEta[m][(i-1)*vpts + jj] = temp2;
                            }
                        if(phi < E->control.fi_min+0.2)
			   if(phi < E->control.fi_min+0.15) {
				temp2 = temp1+(0.01-temp1)*(cos(pow((E->control.fi_min+0.15-phi)/0.15,10)*3.14159)+1)/2;
				if(theta <= E->control.theta_max-0.03 && theta >= E->control.theta_min+0.03)
                                  EEta[m][(i-1)*vpts + jj] = temp2;
			   }
			   else if(theta <= E->control.theta_max-0.05 && theta >= E->control.theta_min+0.05){
                                temp2 = temp1+(0.01-temp1)*(cos(pow((E->control.fi_min+0.15-phi)/0.05,3)*3.14159)+1)/2;
			        if(EEta[m][(i-1)*vpts + jj] > temp2)
                                  EEta[m][(i-1)*vpts + jj] = temp2;
                            } 
                        /*if(r<0.88) {
			        temp2 = temp1+(0.01-temp1)*(cos((r-0.84)/0.04*3.14159)+1)/2;
			        if(EEta[m][(i-1)*vpts + jj] > temp2)
                                  EEta[m][(i-1)*vpts + jj] = temp2;
                         }*/
			// if(EEta[m][(i-1)*vpts + jj] < 0.01)
			//      EEta[m][(i-1)*vpts + jj] = 0.01;
			// end of wk bndry


		    } //end of jj (nodal loop)
		} // end of iii (radious)

		/* add a strong slab core */
		for(iii=E->lmesh.noz;iii>=1;iii--) {
                    node = iii + (jjj-1)*E->lmesh.noz + (kkk-1)*E->lmesh.noz*E->lmesh.nox;
                    theta = E->sx[m][1][node];
                    phi = E->sx[m][2][node];
                    r = E->sx[m][3][node];

		    lat_max1 = lat_max+0.04*(1.0-r)/(1.0-0.975);
		    temp1 = fabs(r-E->sx[m][3][nd_core]);
		    temp = 0.002; // + (1.0-r)/0.02*0.002;
		    i = node-nd_core;
		    if(rheo_trick>=7 && r>0.97 && r<0.997 && theta<=lat_max && theta>=lat_min && temp1<=temp && phi>fi-0.04) {
			for(jj=1;jj<=vpts;jj++) {
			   temp2 = (1.0 + cos(temp1/temp*3.14))/2; // *cos((r-0.97)/0.03*1.57);
                           if(theta>=lat_max1-0.005)
                              temp2 *= cos((theta-lat_max1+0.005)*1.5708/0.005);
                           if(theta<=lat_min+0.005)
                              temp2 *= cos((theta-lat_min-0.005)*1.5708/0.005);
                           temp2 = E->viscosity.N0[0]*(temp2); //define a strong slab core
                           if(theta<=lat_max1 && EEta[m][(el_core-2+i)*vpts+jj]<temp2)
                               EEta[m][(el_core-2+i)*vpts+jj] = temp2;
			} // end of jj
		    } // end of rheo_trick=6
		} // end of for (iii)

            } // end of jjj (latitude)

      if(rheo_trick != 0){
      	for(i=1;i<=dims;i++) {
          free ((void *) PB1[i]);
          free ((void *) PB2[i]);
      	}
      }

      break;

    case 4:

      for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(i=1;i<=nel;i++)   {
          l = E->mat[m][i];
          tempa = E->viscosity.N0[l-1];
          j = 0;

          for(kk=1;kk<=ends;kk++) {
            TT[kk] = E->T[m][E->ien[m][i].node[kk]];
            zz[kk] = (1.-pow(E->sx[m][3][E->ien[m][i].node[kk]],10));
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

/*
            if(E->control.mat_control==0)
              EEta[m][ (i-1)*vpts + jj ] = tempa*
                exp( (E->viscosity.E[l-1] +  E->viscosity.Z[l-1]*zzz )
                         / (E->viscosity.T[l-1]+temp) );



            if(E->control.mat_control==1)
              EEta[m][ (i-1)*vpts + jj ] = tempa*E->VIP[m][i]*
                exp( (E->viscosity.E[l-1] +  E->viscosity.Z[l-1]*zzz )
                         / (E->viscosity.T[l-1]+temp) );
*/

	    // lijun's definition
	    if(E->control.mat_control==0)
              EEta[m][ (i-1)*vpts + jj ] = tempa*
                exp( (E->viscosity.E[l-1] +  E->viscosity.Z[l-1]*zzz )
                         / (E->viscosity.T[l-1]+temp) );

          
          
            if(E->control.mat_control==1)
              EEta[m][ (i-1)*vpts + jj ] = tempa*E->VIP[m][i]*
                exp( (E->viscosity.E[l-1] +  E->viscosity.Z[l-1]*zzz )
                         / (E->viscosity.T[l-1]+temp) );
          
	    // end

	    }
        }
      break;


    case 5:

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
               visc1 = E->VIP[m][i];
               visc2 = 2.0/(1./visc1 + 1.);
               tempa_exp = tempa*
	          exp( E->viscosity.E[l-1]/(temp+E->viscosity.T[l-1])
		     - E->viscosity.E[l-1]/(one +E->viscosity.T[l-1]) );
               visc1 = tempa*E->viscosity.max_value;
               if(tempa_exp > visc1) tempa_exp=visc1;
               EEta[m][ (i-1)*vpts + jj ] = visc2*tempa_exp;
               /* if(E->parallel.me == 0 && visc1 < 1.0e-03)
                  fprintf(stderr,"%f  %f   %e  %e  %e\n",zzz,temp,visc1,visc2,
                          EEta[m][ (i-1)*vpts + jj ]); */
              }

	    }
        }
      break;

    case 6:

    	/* eta = N_0 exp(E/(T+T_0) - E/(1+T_0)) 

	plus weak zones in the lithspere reciprocal of strain rate

        */
	eedot = (float *) malloc((2+nel)*sizeof(float));

        for(m=1;m<=E->sphere.caps_per_proc;m++) {
            for(i=1;i<=nel;i++)   {
                l = E->mat[m][i] - 1;
                if(E->control.mat_control) /* switch moved up here TWB */
                  tempa = E->viscosity.N0[l] * E->VIP[m][i];
                else
                  tempa = E->viscosity.N0[l];

                for(kk=1;kk<=ends;kk++) {
                  TT[kk] = E->T[m][E->ien[m][i].node[kk]];
                }

                for(jj=1;jj<=vpts;jj++) {
                    temp=0.0;
                    for(kk=1;kk<=ends;kk++)   { /* took out
                                                   computation of
                                                   depth, not needed
                                                   TWB */
                      TT[kk]=max(TT[kk],zero);
                      temp += min(TT[kk],one) * E->N.vpt[GNVINDEX(kk,jj)];
                    }
                    EEta[m][ (i-1)*vpts + jj ] = tempa*
                      exp( E->viscosity.E[l]/(temp+E->viscosity.T[l])
                           - E->viscosity.E[l]/(one +E->viscosity.T[l]) );
                }
            }

            // add low viscosity zones reciprocal of strain rate
            /* get second invariant for all elements */
	    if(E->monitor.solution_cycles == 0)
                for(i=1;i<=nel;i++) /* initialize with unity if no velocities around */
                    eedot[i] = 1.0;
            else {
                strain_rate_2_inv(E,m,eedot,1);
                t1 = 0.0;
                for(i=1;i<=nel;i++) { /* eedot cannot be too small, or the viscosity will go to inf */
                   eedot[i] = max(eedot[i], 1.0e-16);
                   if(t1<eedot[i]) t1 = eedot[i]; // && E->mat[m][i]==1)  t1 = eedot[i];
                }
                MPI_Allreduce(&t1, &t2,1,MPI_DOUBLE,MPI_MAX,E->parallel.world);
                if(E->parallel.me == 0) {
                  if(E->monitor.solution_cycles == 1) {
                    fp=fopen("Max_strain_rate","w");
                    fprintf(fp,"%f\n",t2);
                    fclose(fp);
                  }
                  else {
                    fp=fopen("Max_strain_rate","r");
                    fscanf("%f",&t2);
                    fclose(fp);
                  }
                  fprintf(E->fp,"peak strain rate @ proc %d = %f\n", E->parallel.me, t1);
                  fprintf(E->fp,"global max strain rate = %f\n", t2);
                }
	        /* normalization */
                if(t2 > 1.0) {
           	  for(i=1;i<=nel;i++) {
             	     eedot[i] = pow(eedot[i]/t2,3)*1e5;
             	     if(eedot[i]>1.0) eedot[i] = pow(eedot[i],0.33);
           	  }
                }

                // eta = pow(eta,1/n)*pow(epcilon_dot,(1-n)/n)
	        fprintf(E->fp,"the_max=%f,the_min=%f\n",E->control.theta_max,E->control.theta_min);
	        fprintf(E->fp,"fi_max=%f,fi_min=%f\n",E->control.fi_max,E->control.fi_min);
                for(i=1;i<=nel;i++)   {
            	  exponent1= one/E->viscosity.sdepv_expt[E->mat[m][i]-1];
            	  //scale=pow(eedot[e],exponent1-one);
            	  scale=1.0/eedot[i];
		  node = E->ien[m][i].node[1];
            	  for(jj=1;jj<=vpts;jj++) {
                    if(scale < 1 && E->mat[m][i]==1) 
			if(E->sx[m][1][node]<E->control.theta_max-0.2 && E->sx[m][1][node]>E->control.theta_min+0.2
		  	   && E->sx[m][2][node]<E->control.fi_max-0.3 && E->sx[m][2][node]>E->control.fi_min+0.25)
                    		EEta[m][(i-1)*vpts + jj] = scale*pow(EEta[m][(i-1)*vpts+jj],exponent1);
                    		//EEta[m][(e-1)*vpts + jj] *= scale;
		    if(E->sx[m][1][node]>E->control.theta_max-0.1)
			EEta[m][(i-1)*vpts + jj] *= 0.001+(1.0-cos((E->sx[m][1][node]-E->control.theta_max)/0.1*1.57));
		    if(E->sx[m][1][node]<E->control.theta_min+0.1)
			EEta[m][(i-1)*vpts + jj] *= 0.001+(1.0-cos((E->sx[m][1][node]-E->control.theta_min)/0.1*1.57));
		    if(E->sx[m][2][node]>E->control.fi_max-0.1)
			EEta[m][(i-1)*vpts + jj] *= 0.001+(1.0-cos((E->sx[m][2][node]-E->control.fi_max)/0.1*1.57));
		    if(E->sx[m][2][node]<E->control.fi_min+0.1)
                        EEta[m][(i-1)*vpts + jj] *= 0.001+(1.0-cos((E->sx[m][2][node]-E->control.fi_min)/0.1*1.57));
		  }
		}
            } // end of else (i.e., E->monitor.solution_cycles > 0)
    	} // end of m

    	free ((void *)eedot);
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
      strain_rate_2_inv(E,m,eedot,1);

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



void strain_rate_2_inv(E,m,EEDOT,SQRT)
     struct All_variables *E;
     float *EEDOT;
     int m,SQRT;
{
    void get_global_shape_fn();
    void velo_from_element();

    struct Shape_function GN;
    struct Shape_function_dA dOmega;
    struct Shape_function_dx GNx;

    double edot[4][4],dudx[4][4],rtf[4][9];
    float VV[4][9];

    int e,i,p,q,n,nel,k;

    const int dims = E->mesh.nsd;
    const int ends = enodes[dims];
    const int lev = E->mesh.levmax;
    const int nno = E->lmesh.nno;
    const int vpts = vpoints[dims];
    const int sphere_key = 0;

    nel = E->lmesh.nel;

    for(e=1;e<=nel;e++) {

      get_global_shape_fn(E,e,&GN,&GNx,&dOmega,2,sphere_key,rtf,lev,m);

      velo_from_element(E,VV,m,e,sphere_key);

      for(p=1;p<=dims;p++)
        for(q=1;q<=dims;q++)
           dudx[p][q] = 0.0;

      for(i=1;i<=ends;i++)
        for(p=1;p<=dims;p++)
           for(q=1;q<=dims;q++)
              dudx[p][q] += VV[p][i] * GNx.ppt[GNPXINDEX(q-1,i,1)];

      for(p=1;p<=dims;p++)
        for(q=1;q<=dims;q++)
            edot[p][q] = dudx[p][q] + dudx[q][p];

      if (dims==2)
         EEDOT[e] = edot[1][1]*edot[1][1] + edot[2][2]*edot[2][2]
                  + edot[1][2]*edot[1][2]*2.0;

      else if (dims==3)
         EEDOT[e] = edot[1][1]*edot[1][1] + edot[1][2]*edot[1][2]*2.0
                  + edot[2][2]*edot[2][2] + edot[2][3]*edot[2][3]*2.0
                  + edot[3][3]*edot[3][3] + edot[1][3]*edot[1][3]*2.0;

      }

    if(SQRT)
	for(e=1;e<=nel;e++)
	    EEDOT[e] =  sqrt(0.5 *EEDOT[e]);
    else
	for(e=1;e<=nel;e++)
	    EEDOT[e] *=  0.5;

    return;
}



void visc_to_node_interpolate(E,evisc,visc)
 struct All_variables *E;
 float **evisc,**visc;
{

/*  void exchange_node_f(); */
/*  void get_global_shape_fn(); */
/*  void return_horiz_ave_f(); */
/*  void sphere_interpolate(); */
/*  void print_interpolated(); */
/*  void gather_TG_to_me0(); */
/*  void parallel_process_termination(); */
/*  int i,j,k,e,node,snode,m,nel2; */
/*    FILE *fp; */
/*    char output_file[255]; */

/*  float *TG,t,f,rad, Szz; */

/*  double time1,CPU_time0(),tww[9],rtf[4][9]; */

/*  struct Shape_function GN; */
/*  struct Shape_function_dA dOmega; */
/*  struct Shape_function_dx GNx; */

/*  const int dims=E->mesh.nsd,dofs=E->mesh.dof; */
/*  const int vpts=vpoints[dims]; */
/*  const int ppts=ppoints[dims]; */
/*  const int ends=enodes[dims]; */
/*  const int nno=E->lmesh.nno; */
/*  const int lev=E->mesh.levmax; */


/*     TG =(float *)malloc((E->sphere.nsf+1)*sizeof(float)); */
/*     for (i=E->sphere.nox;i>=1;i--) */
/*       for (j=1;j<=E->sphere.noy;j++)  { */
/*            node = i + (j-1)*E->sphere.nox; */
/* 	   TG[node] = 0.0; */
/*   	   m = E->sphere.int_cap[node]; */
/* 	   e = E->sphere.int_ele[node]; */

/* 	   if (m>0 && e>0) { */
/* 	      e=e+E->lmesh.elz-1; */
/* 	      TG[node] = log10(evisc[m][(e-1)*vpts+1]); */
/* 	      } */
/* 	   } */

/*     gather_TG_to_me0(E,TG); */

/*     if (E->parallel.me==E->parallel.nprocz-1)  { */
/*      sprintf(output_file,"%s.evisc_intp",E->control.data_file); */
/*      fp=fopen(output_file,"w"); */

/*     rad = 180/M_PI; */
/*     for (i=E->sphere.nox;i>=1;i--) */
/*       for (j=1;j<=E->sphere.noy;j++)  { */
/*            node = i + (j-1)*E->sphere.nox; */
/*            t = 90-E->sphere.sx[1][node]*rad; */
/* 	   f = E->sphere.sx[2][node]*rad; */
/* 	   fprintf (fp,"%.3e %.3e %.4e\n",f,t,TG[node]); */
/* 	   } */
/*       fclose(fp); */
/*      } */

/*  free((void *)TG); */

   return;
   }

float theta_locater(PB,n,theta)
    float **PB,theta;
    int n;
{
    int flag,flag1,i;
    float temp1,fi;

    flag = 0;
    flag1 = 0;
    temp1 = 0.3;
    for(i=1;i<n;i++)
        if(fabs(PB[1][i]-theta)<=temp1) {
            flag = i;
            temp1 = fabs(PB[1][i]-theta);
        }
    temp1 = 0.3;
    for(i=1;i<n;i++)
        if(fabs(PB[1][i]-theta)<=temp1 && i!=flag) {
           flag1 = i;
           temp1 = fabs(PB[1][i]-theta);
        }
     //assign the longitude value to the plate-boundary marker (fi)
     if(flag) {
        if(flag1 && (PB[1][flag1]-PB[1][flag]))
            fi = PB[2][flag]+(PB[2][flag1]-PB[2][flag])*(theta-PB[1][flag])/(PB[1][flag1]-PB[1][flag])+2*3.1415926;
        else
            fi = PB[2][flag]+2*3.1415926;
     }
     else
        if(flag1) fi = PB[2][flag1]+2*3.1415926;

    return(fi);
}

