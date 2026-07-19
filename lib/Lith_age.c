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

#include <math.h>

#include "global_defs.h"

/*#include "age_related.h"*/
#include "parallel_related.h"
#include "parsing.h"
#include "lith_age.h"

float find_age_in_MY();
void lith_age_update_tbc(struct All_variables *E);

static float effective_plate_age_nd(const struct All_variables *E,
                                    float age_nd)
{
  const float max_age_nd = E->control.max_plate_age_Ma / E->data.scalet;

  return age_nd > max_age_nd ? max_age_nd : age_nd;
}


void lith_age_input(struct All_variables *E)
{
  int m = E->parallel.me;

  E->control.lith_age = 0;
  E->control.lith_age_time = 0;
  E->control.temperature_bound_adj = 0;

  input_int("lith_age",&(E->control.lith_age),"0",m);
  input_float("mantle_temp",&(E->control.lith_age_mantle_temp),"1.0",m);
  input_float("bottom_tbl_thickness",&(E->control.bottom_tbl_thickness),"0.0",m);
  input_float("bottom_tbl_diffusivity_ratio",&(E->control.bottom_tbl_diffusivity_ratio),"1.0",m);

  if (E->control.lith_age) {
    input_int("lith_age_time",&(E->control.lith_age_time),"0",m);
    input_string("lith_age_file",E->control.lith_age_file,"",m);
    input_float("lith_age_depth",&(E->control.lith_age_depth),"0.0471",m);
    input_float("max_plate_age_Ma",&(E->control.max_plate_age_Ma),"70.0",m);
    if(E->control.max_plate_age_Ma <= 0.0) {
      fprintf(stderr,"max_plate_age_Ma must be positive; found %g\n",
              E->control.max_plate_age_Ma);
      parallel_process_termination();
    }

    input_int("temperature_bound_adj",&(E->control.temperature_bound_adj),"0",m);
    if (E->control.temperature_bound_adj) {
      input_float("depth_bound_adj",&(E->control.depth_bound_adj),"0.1570",m);
      input_float("width_bound_adj",&(E->control.width_bound_adj),"0.08727",m);
    }
  }
  return;
}


void lith_age_init(struct All_variables *E)
{
  char output_file[255];
  FILE *fp1;
  int node, i, j, output;

  int gnox, gnoy;
  gnox=E->mesh.nox;
  gnoy=E->mesh.noy;

  //if (E->parallel.me == 0 ) fprintf(stderr,"INSIDE lith_age_init\n");
  E->age_t=(float*) malloc((gnox*gnoy+1)*sizeof(float));

  if(E->control.lith_age_time==1)   {
    /* to open files every timestep */
    E->control.lith_age_old_cycles = E->monitor.solution_cycles;
    output = 1;
    (E->solver.lith_age_read_files)(E,output);
  }
  else {
    /* otherwise, just open for the first timestep */
    /* NOTE: This is only used if we are adjusting the boundaries */
    sprintf(output_file,"%s",E->control.lith_age_file);
    fp1=fopen(output_file,"r");
    if (fp1 == NULL) {
      fprintf(E->fp,"(Boundary_conditions #1) Can't open %s\n",output_file);
      parallel_process_termination();
    }
    for(i=1;i<=gnoy;i++)
      for(j=1;j<=gnox;j++) {
	node=j+(i-1)*gnox;
	fscanf(fp1,"%f",&(E->age_t[node]));
	/* Input is Ma.  scalet is Ma per nondimensional thermal time. */
	E->age_t[node]=E->age_t[node]/E->data.scalet;
      }
    fclose(fp1);
  } /* end E->control.lith_age_time == false */
}


void lith_age_construct_tic(struct All_variables *E)
{
  int i, j, k, m, node, nodeg, radial_index;
  int nox, noy, noz, gnox;
  double r1, temp, temp1, temp2, temp3, temp4, dist, dist2;
  double bottom_tbl_age_nd, bottom_tbl_age_ma, bottom_tbl_delta, eta;
  double background_temperature, background_cmb_k, background_surface_k;
  double surface_hsc_delta;
  const double bottom_tbl_erf_inv = 1.8214;
  float age,age1,theta,phi,assim_depth,asm_depth,flag_depth2;
  void temperatures_conform_bcs();

  noy=E->lmesh.noy;
  nox=E->lmesh.nox;
  noz=E->lmesh.noz;

  gnox=E->mesh.nox;
  if(!E->refstate.has_temperature) {
    fprintf(stderr,
            "lith_age_construct_tic requires an initial background geotherm "
            "(strict ALA column 3; legacy extended column 7)\n");
    parallel_process_termination();
  }
  if(E->control.ala_pressure_buoyancy &&
     (fabs(E->control.TBCtopval) > 1.0e-6 ||
      fabs(E->control.TBCbotval-1.0) > 1.0e-6)) {
    fprintf(stderr,
            "Strict ALA Ttop/Tbottom scaling requires toptbcval=0 and "
            "bottbcval=1; found %e and %e\n",
            E->control.TBCtopval, E->control.TBCbotval);
    parallel_process_termination();
  }

  background_cmb_k = E->data.Ttop
                     + E->refstate.temperature_cmb*E->data.ref_temperature;
  background_surface_k = E->data.Ttop
                         + E->refstate.temperature_surface
                           *E->data.ref_temperature;
  surface_hsc_delta = (background_surface_k-E->data.Ttop)
                      / E->data.ref_temperature;
  if(E->control.ala_pressure_buoyancy &&
     (!isfinite(background_cmb_k) || !isfinite(background_surface_k) ||
      background_cmb_k >= E->data.Tbottom ||
      background_surface_k <= E->data.Ttop)) {
    fprintf(stderr,
            "Invalid K-background endpoints for anomaly superposition: "
            "T_K_CMB=%f Tbottom=%f T_K_surface=%f Ttop=%f\n",
            background_cmb_k, E->data.Tbottom,
            background_surface_k, E->data.Ttop);
    parallel_process_termination();
  }

  /* Establish the complete initial Katsura background geotherm first. This is
     initialization data, not an active thermodynamic ALA reference profile.
     Interior samples are copied directly; only the endpoint-closed strict-ALA
     rows use the unclosed values recovered by read_refstate(). */
  for(m=1;m<=E->sphere.caps_per_proc;m++)
    for(i=1;i<=noy;i++)
      for(j=1;j<=nox;j++)
	for(k=1;k<=noz;k++) {
	  node=k+(j-1)*noz+(i-1)*nox*noz;
	  radial_index = k+E->lmesh.nzs-1;
	  background_temperature = E->refstate.temperature[k];
	  if(E->control.ala_pressure_buoyancy && radial_index == 1)
	    background_temperature = E->refstate.temperature_cmb;
	  else if(E->control.ala_pressure_buoyancy &&
	          radial_index == E->mesh.noz)
	    background_temperature = E->refstate.temperature_surface;
	  E->T[m][node] = background_temperature;
	}

  /* Add the MATLAB-style bottom TBL as a decaying CMB perturbation on the
     local Katsura background.  H_TBL is the height where about 1 percent
     of the CMB anomaly remains, so sqrt(kappa*t)=H_TBL/(2*1.8214).
     It is not a cutoff; erfc() decays naturally above that height. */
  if(E->control.bottom_tbl_thickness > 0.0) {
    if(E->control.bottom_tbl_thickness > E->sphere.ro-E->sphere.ri ||
       E->control.bottom_tbl_diffusivity_ratio <= 0.0) {
      fprintf(stderr,
              "Invalid bottom TBL controls: thickness=%e "
              "diffusivity_ratio=%e\n",
              E->control.bottom_tbl_thickness,
              E->control.bottom_tbl_diffusivity_ratio);
      parallel_process_termination();
    }

    bottom_tbl_age_nd = pow(E->control.bottom_tbl_thickness /
                            (2.0 * bottom_tbl_erf_inv), 2.0);
    bottom_tbl_age_ma = bottom_tbl_age_nd * E->data.scalet /
                        E->control.bottom_tbl_diffusivity_ratio;
    /* Ttop/Tbottom own the physical scale.  bottbcval is validated above but
       is not allowed to redefine the physical CMB temperature. */
    bottom_tbl_delta = (E->data.Tbottom-background_cmb_k)
                       / E->data.ref_temperature;

    if(E->parallel.me == 0) {
      fprintf(stderr,
              "Bottom TBL: H=%e (%f km) diffusivity_ratio=%f "
              "inferred_age=%f Ma T_bottom_K=%f T_K_CMB=%f "
              "delta_T_K=%f delta_T_nd=%e\n",
              E->control.bottom_tbl_thickness,
              E->control.bottom_tbl_thickness * E->data.radius_km,
              E->control.bottom_tbl_diffusivity_ratio,
              bottom_tbl_age_ma, E->data.Tbottom, background_cmb_k,
              E->data.Tbottom-background_cmb_k, bottom_tbl_delta);
    }

    for(m=1;m<=E->sphere.caps_per_proc;m++)
      for(i=1;i<=noy;i++)
        for(j=1;j<=nox;j++)
	  for(k=1;k<=noz;k++) {
	    node=k+(j-1)*noz+(i-1)*nox*noz;
	    r1=E->sx[m][3][node];
	    eta=(r1-E->sphere.ri)*0.5/sqrt(bottom_tbl_age_nd);
	    E->T[m][node] += bottom_tbl_delta * erfc(eta);
	  }
  }

  /* Preserve the existing shallow HSC, initial-slab, mask and BC overlay. */
  for(m=1;m<=E->sphere.caps_per_proc;m++)
    for(i=1;i<=noy;i++)
      for(j=1;j<=nox;j++)
	for(k=1;k<=noz;k++)  {
	  nodeg=E->lmesh.nxs-1+j+(E->lmesh.nys+i-2)*gnox;
	  node=k+(j-1)*noz+(i-1)*nox*noz;
	  r1=E->sx[m][3][node];

          age1=E->age_t[nodeg];

	  dist=fabs(E->flag_depth2[nodeg])*1000.0/6371.0;
	  if(dist>0.1 && age1>1000.0/E->data.scalet) /* Jiashun set craton to 2500 Ma */
	      asm_depth=1.0;
	  else if(E->flag_depth2[nodeg]>=0.0)
	      asm_depth=1.0;
          else if (dist<=0.1)
              asm_depth=0.5;
          else
              asm_depth=1.0;
	  assim_depth = asm_depth*E->control.lith_age_depth;

          age1=effective_plate_age_nd(E,age1);
	  if( r1 >= E->sphere.ro-assim_depth)
	    { /* if closer than (lith_age_depth) from top */
	      /* zero age surface has mantle temperature */
	      if(age1>0.0) //Lijun: define continent to be 0 age
                    temp = (E->sphere.ro-r1) *0.5 /sqrt(age1);
                else
                    temp = 10.0;
	      if(E->control.ala_pressure_buoyancy) {
	        /* Surface HSC is an anomaly superposed on the local K background.
	           The dimensional operation T_K-A_surface is algebraically
	           equivalent to subtracting its normalized amplitude here. */
	        temp1 = surface_hsc_delta * erfc(temp);
		E->T[m][node] -= temp1;
	      }
	      else {
	        temp1 = E->control.lith_age_mantle_temp * erf(temp);
		E->T[m][node] = temp1;
	      }
	    }
	/* add initial slabs */
	flag_depth2=E->flag_depth2[nodeg]*1000.0/6371.0;
	if(flag_depth2<0.0 && flag_depth2>-0.05 && r1>=0.95) {
                temp1=0.0025-flag_depth2*flag_depth2;
                temp2=0.046*0.046-flag_depth2*flag_depth2;
                temp3=0.0009-flag_depth2*flag_depth2;
                temp4=(r1-0.95)*(r1-0.95);
                if(temp4<=temp1 && temp4>=temp3) {
                    dist2=0.05-sqrt(temp4+flag_depth2*flag_depth2);
                    temp=dist2*0.5/sqrt(50.0/E->data.scalet);
                    E->T[m][node]=E->control.lith_age_mantle_temp * erf(temp);
                }
           }

	}

  /* modify temperature BC to be concorded with read in T */
  lith_age_update_tbc(E);

  temperatures_conform_bcs(E);

  return;
}


void lith_age_update_tbc(struct All_variables *E)
{
  int i, j, k, m, node;
  int nox, noy, noz;
  double r1, rout, rin;
  const float e_4=1.e-4;

  noy = E->lmesh.noy;
  nox = E->lmesh.nox;
  noz = E->lmesh.noz;
  rout = E->sphere.ro;
  rin = E->sphere.ri;

  for(m=1;m<=E->sphere.caps_per_proc;m++)
    for(i=1;i<=noy;i++)
      for(j=1;j<=nox;j++)
	for(k=1;k<=noz;k++)  {
	  node=k+(j-1)*noz+(i-1)*nox*noz;
	  r1=E->sx[m][3][node];

	  if(fabs(r1-rout)>=e_4 && fabs(r1-rin)>=e_4)  {
	    E->sphere.cap[m].TB[1][node]=E->T[m][node];
	    E->sphere.cap[m].TB[2][node]=E->T[m][node];
	    E->sphere.cap[m].TB[3][node]=E->T[m][node];
	  }
	}

  return;
}


void lith_age_temperature_bound_adj(struct All_variables *E, int lv)
{
  int i,j,node,nno,ii,jj,kk,nodeg,nnn,ttt,intage;
  float tt1,tt2,ttt1,ttt2,ttt3,fff2,fff3,*PB1[4],*PB2[4],fi,fi_1,fi_2,gap,assim_depth,asm_depth,dist;
  float theta,phi,lat_max,lat_max1,lat_max2,lat_min,lat_min1,lat_min2,age,age1,age2;
  char pb_1[255],pb_2[255];
  FILE *fp1,*fp2;
  float find_age_in_MY();
  const int dims=E->mesh.nsd;

  float theta_locater();

  nno=E->lmesh.nno;

/* NOTE: To start, the relevent bits of "node" are zero. Thus, they only
get set to TBX/TBY/TBZ if the node is in one of the bounding regions.
Also note that right now, no matter which bounding region you are in,
all three get set to true. CPC 6/20/00 */

  if (E->control.temperature_bound_adj) {
    ttt2=E->control.theta_min + E->control.width_bound_adj;
    ttt3=E->control.theta_max - E->control.width_bound_adj;
    fff2=E->control.fi_min + E->control.width_bound_adj;
    fff3=E->control.fi_max - E->control.width_bound_adj;

    if(lv==E->mesh.gridmax)
      for(j=1;j<=E->sphere.caps_per_proc;j++)
	for(node=1;node<=E->lmesh.nno;node++)  {
	  if( ((E->sx[j][1][node]<=ttt2) && (E->sx[j][3][node]>=E->sphere.ro-E->control.depth_bound_adj)) || ((E->sx[j][1][node]>=ttt3) && (E->sx[j][3][node]>=E->sphere.ro-E->control.depth_bound_adj)) )
	    /* if < (width) from x bounds AND (depth) from top */
	    {
	      E->node[j][node]=E->node[j][node] | TBX;
	      E->node[j][node]=E->node[j][node] & (~FBX);
	      E->node[j][node]=E->node[j][node] | TBY;
	      E->node[j][node]=E->node[j][node] & (~FBY);
	      E->node[j][node]=E->node[j][node] | TBZ;
	      E->node[j][node]=E->node[j][node] & (~FBZ);
	    }

	  if( ((E->sx[j][2][node]<=fff2) && (E->sx[j][3][node]>=E->sphere.ro-E->control.depth_bound_adj)) )
	    /* if fi is < (width) from side AND z is < (depth) from top */
	    {
	      E->node[j][node]=E->node[j][node] | TBX;
	      E->node[j][node]=E->node[j][node] & (~FBX);
	      E->node[j][node]=E->node[j][node] | TBY;
	      E->node[j][node]=E->node[j][node] & (~FBY);
	      E->node[j][node]=E->node[j][node] | TBZ;
	      E->node[j][node]=E->node[j][node] & (~FBZ);
	    }

	  if( ((E->sx[j][2][node]>=fff3) && (E->sx[j][3][node]>=E->sphere.ro-E->control.depth_bound_adj)) )
	    /* if fi is < (width) from side AND z is < (depth) from top */
	    {
	      E->node[j][node]=E->node[j][node] | TBX;
	      E->node[j][node]=E->node[j][node] & (~FBX);
	      E->node[j][node]=E->node[j][node] | TBY;
	      E->node[j][node]=E->node[j][node] & (~FBY);
	      E->node[j][node]=E->node[j][node] | TBZ;
	      E->node[j][node]=E->node[j][node] & (~FBZ);
	    }

	}
  } /* end E->control.temperature_bound_adj */

  //if (E->control.lith_age_time && E->data.timedir>0) {
  if (E->control.lith_age_time) {

    if(lv==E->mesh.gridmax)
      for(j=1;j<=E->sphere.caps_per_proc;j++)
	//for(node=1;node<=E->lmesh.nno;node++)  {
	for(ii=1;ii<=E->lmesh.noy;ii++)
          for(jj=1;jj<=E->lmesh.nox;jj++) {
	    node = 1 + (jj-1)*E->lmesh.noz + (ii-1)*E->lmesh.noz*E->lmesh.nox;
            for(kk=1;kk<=E->lmesh.noz;kk++)  {
                nodeg=E->lmesh.nxs-1+jj+(E->lmesh.nys+ii-2)*E->mesh.nox;
                node=kk+(jj-1)*E->lmesh.noz+(ii-1)*E->lmesh.nox*E->lmesh.noz;
		theta = E->sx[j][1][node];
        	phi = E->sx[j][2][node];
		dist=fabs(E->flag_depth2[nodeg]);
		if(E->flag_depth2[nodeg]*1000.0/6371.0>=-0.1 && E->flag_depth2[nodeg]<=0.0 && dist*1000.0/6371.0<=0.1)
			asm_depth=0.003;
	        else if (dist<0.06)
			asm_depth=0.003;
                else if (dist<=0.1)
                        asm_depth=0.003+(1.0-0.003)*(dist-0.06)/(0.1-0.06);
                else
                        asm_depth=1.0;
                assim_depth=asm_depth*E->control.lith_age_depth;

		/*if(theta<E->control.theta_max-0.10 && theta>E->control.theta_min+0.10
                      && phi<E->control.fi_max-0.10 && phi>E->control.fi_min+0.10) {*/
                    if(E->sx[j][3][node]>=E->sphere.ro-assim_depth) {
		    /*if(E->sx[j][3][node]>=E->sphere.ro-E->control.lith_age_depth) { */
                    // if closer than (lith_age_depth) from top
                        E->node[j][node]=E->node[j][node] | TBX;
                        E->node[j][node]=E->node[j][node] & (~FBX);
                        E->node[j][node]=E->node[j][node] | TBY;
                        E->node[j][node]=E->node[j][node] & (~FBY);
                        E->node[j][node]=E->node[j][node] | TBZ;
                        E->node[j][node]=E->node[j][node] & (~FBZ);
		    }
		//}
		else if(E->sx[j][3][node] > E->sphere.ri + 1.e-4) {
                        E->node[j][node]=E->node[j][node] | FBX;
                        E->node[j][node]=E->node[j][node] & (~TBX);
                        E->node[j][node]=E->node[j][node] | FBY;
                        E->node[j][node]=E->node[j][node] & (~TBY);
                        E->node[j][node]=E->node[j][node] | FBZ;
                        E->node[j][node]=E->node[j][node] & (~TBZ);
                    }

            }//end of for(kk)
	  }//end of for(jj)

  } // end E->control.lith_age_time
  

  return;
}


void lith_age_conform_tbc(struct All_variables *E)
{
  int m,j,node,nox,noz,noy,gnox,gnoy,gnoz,nodeg,i,k;
  float ttt2,ttt3,fff2,fff3;
  float r1,t1,f1,t0,temp,temp1,age1,norm;
  float depth,assim_depth,asm_depth,dist;
  float e_4;
  FILE *fp1;
  char output_file[255];
  int output;


  e_4=1.e-4;
  output = 0;

  gnox=E->mesh.nox;
  gnoy=E->mesh.noy;
  gnoz=E->mesh.noz;
  nox=E->lmesh.nox;
  noy=E->lmesh.noy;
  noz=E->lmesh.noz;

  if(E->control.lith_age_time==1)   {
    /* to open files every timestep */
    if (E->control.lith_age_old_cycles != E->monitor.solution_cycles) {
      /*update so that output only happens once*/
      output = 1;
      E->control.lith_age_old_cycles = E->monitor.solution_cycles;
    }
    if (E->parallel.me == 0) fprintf(stderr,"INSIDE lith_age_conform_tbc\n");
    (E->solver.lith_age_read_files)(E,output);
  }

  /* NOW SET THE TEMPERATURES IN THE BOUNDARY REGIONS */
  if(E->monitor.solution_cycles>1 && E->control.temperature_bound_adj) {
    ttt2=E->control.theta_min + E->control.width_bound_adj;
    ttt3=E->control.theta_max - E->control.width_bound_adj;
    fff2=E->control.fi_min + E->control.width_bound_adj;
    fff3=E->control.fi_max - E->control.width_bound_adj;

    for(m=1;m<=E->sphere.caps_per_proc;m++)
      for(i=1;i<=noy;i++)
	for(j=1;j<=nox;j++)
	  for(k=1;k<=noz;k++)  {
	    nodeg=E->lmesh.nxs-1+j+(E->lmesh.nys+i-2)*gnox;
	    node=k+(j-1)*noz+(i-1)*nox*noz;
	    t1=E->sx[m][1][node];
	    f1=E->sx[m][2][node];
	    r1=E->sx[m][3][node];

	    if(fabs(r1-E->sphere.ro)>=e_4 && fabs(r1-E->sphere.ri)>=e_4)  { // if NOT right on the boundary 
	      if( ((E->sx[m][1][node]<=ttt2) && (E->sx[m][3][node]>=E->sphere.ro-E->control.depth_bound_adj)) || ((E->sx[m][1][node]>=ttt3) && (E->sx[m][3][node]>=E->sphere.ro-E->control.depth_bound_adj))) {
		// if < (width) from x bounds AND (depth) from top 
		age1=effective_plate_age_nd(E,E->age_t[nodeg]);
		if(age1>0.0) //Lijun define continent to be 0 age
		    temp = (E->sphere.ro-r1) *0.5 /sqrt(age1)/3;
		else
		    temp = 10.0;
		//temp = 10.0;
		t0 = E->control.lith_age_mantle_temp * erf(temp);

		// keep the age the same!
		E->sphere.cap[m].TB[1][node] = t0;
		E->sphere.cap[m].TB[2][node] = t0;
		E->sphere.cap[m].TB[3][node] = t0;
	      }

	      if( ((E->sx[m][2][node]<=fff2) || (E->sx[m][2][node]>=fff3)) && (E->sx[m][3][node]>=E->sphere.ro-E->control.depth_bound_adj)) {
		// if < (width) from y bounds AND (depth) from top 


		// keep the age the same!
		age1=effective_plate_age_nd(E,E->age_t[nodeg]);
		if(age1>0.0)
		    temp = (E->sphere.ro-r1) *0.5 /sqrt(age1)/3;
		else
                    temp = 10.0;
		//temp = 10.0;
		t0 = E->control.lith_age_mantle_temp * erf(temp);

		E->sphere.cap[m].TB[1][node]=t0;
		E->sphere.cap[m].TB[2][node]=t0;
	 	E->sphere.cap[m].TB[3][node]=t0;
	      }

	    }

	  } // end k   

  }   //  end of solution cycles  && temperature_bound_adj
  


  /* NOW SET THE TEMPERATURES IN THE LITHOSPHERE IF CHANGING EVERY TIME STEP */
  if(E->monitor.solution_cycles>0 && E->control.lith_age_time)   {
    for(m=1;m<=E->sphere.caps_per_proc;m++)
      for(i=1;i<=noy;i++)
	for(j=1;j<=nox;j++)
	  for(k=1;k<=noz;k++)  {
	    nodeg=E->lmesh.nxs-1+j+(E->lmesh.nys+i-2)*gnox;
	    node=k+(j-1)*noz+(i-1)*nox*noz;
	    t1=E->sx[m][1][node];
	    f1=E->sx[m][2][node];
	    r1=E->sx[m][3][node];

            age1=E->age_t[nodeg];

	    if(fabs(r1-E->sphere.ro)>=e_4 && fabs(r1-E->sphere.ri)>=e_4)  { // if NOT right on the boundary 
	      dist=fabs(E->flag_depth2[nodeg]);
	      if(E->flag_depth2[nodeg]*1000.0/6371.0>=-0.1 && E->flag_depth2[nodeg]<=0.0 && dist*1000.0/6371.0<=0.1)
			asm_depth=0.003;
	      else if (dist<0.06)
                        asm_depth=0.003;
              else if (dist<=0.1)
                        asm_depth=0.003+(1.0-0.003)*(dist-0.06)/(0.1-0.06);
              else
                        asm_depth=1.0;
	      assim_depth=asm_depth*E->control.lith_age_depth;
	      age1=effective_plate_age_nd(E,age1);
	      if(  E->sx[m][3][node]>=E->sphere.ro-assim_depth ) {
		// if closer than (lith_age_depth) from top 

                depth=E->sphere.ro - E->sx[m][3][node];

		// set a new age from the file 
		if(age1>0.0) {
                    temp = (E->sphere.ro-r1) *0.5 /sqrt(age1);
		    t0 = E->control.lith_age_mantle_temp * erf(temp);
		}
                else {
                    temp = 10.0;
		    t0 = E->control.lith_age_mantle_temp * erf(temp);
		}
		//t0 = E->control.lith_age_mantle_temp;

		E->sphere.cap[m].TB[1][node]=t0;
		E->sphere.cap[m].TB[2][node]=t0;
		E->sphere.cap[m].TB[3][node]=t0; 
	      }
	    }
	  }     // end k  
  }   //  end of solution cycles  && lith_age_time
  

  return;
}


void assimilate_lith_conform_bcs(struct All_variables *E)
{
  float depth, dist, dist2, daf, assimilate_new_temp,temp,temp1,temp2,temp3,temp4,fi_1,fi_2,fi,*PB1[4],*PB2[4],asm_depth,assim_depth,flag_depth2;
  float theta,phi,lat_max,lat_max1,lat_max2,lat_min,lat_min1,lat_min2,age,age1,age2,rad;
  float tt1,tt2,ttt1,ttt2,gap,wid_assim,v_trench;
  int m,j,nno,node,nox,noz,noy,gnox,gnoy,gnoz,nodeg,ii,i,k,nnn,ttt,intage;
  char pb_1[255],pb_2[255];
  FILE *fp1,*fp2;
  float find_age_in_MY();
  const int dims=E->mesh.nsd;
  unsigned int type;

  nno=E->lmesh.nno;
  gnox=E->mesh.nox;
  gnoy=E->mesh.noy;
  gnoz=E->mesh.noz;
  nox=E->lmesh.nox;
  noy=E->lmesh.noy;
  noz=E->lmesh.noz;

  age=find_age_in_MY(E);
  intage=age;

  for(j=1;j<=E->sphere.caps_per_proc;j++)
      for(i=1;i<=noy;i++)
        for(m=1;m<=nox;m++) 
          for(k=1;k<=noz;k++)  {
	    nodeg=E->lmesh.nxs-1+m+(E->lmesh.nys+i-2)*gnox;
            node=k+(m-1)*noz+(i-1)*nox*noz;
	    theta = E->sx[j][1][node];
            phi = E->sx[j][2][node];
	
            type = (E->node[j][node] & (TBX | TBZ | TBY));

            switch (type) {
            case 0:  /* no match, next node */
                break;
            case TBX:
            	assimilate_new_temp = E->sphere.cap[j].TB[1][node];
            	break;
            case TBZ:
            	assimilate_new_temp = E->sphere.cap[j].TB[3][node];
            	break;
            case TBY:
            	assimilate_new_temp = E->sphere.cap[j].TB[2][node];
            	break;
            case (TBX | TBZ):     /* clashes ! */
            	assimilate_new_temp = 0.5 * (E->sphere.cap[j].TB[1][node] + E->sphere.cap[j].TB[3][node]);
            	break;
            case (TBX | TBY):     /* clashes ! */
            	assimilate_new_temp = 0.5 * (E->sphere.cap[j].TB[1][node] + E->sphere.cap[j].TB[2][node]);
            	break;
            case (TBZ | TBY):     /* clashes ! */
            	assimilate_new_temp = 0.5 * (E->sphere.cap[j].TB[3][node] + E->sphere.cap[j].TB[2][node]);
            	break;
            case (TBZ | TBY | TBX):     /* clashes ! */
            	assimilate_new_temp = 0.3333333 * (E->sphere.cap[j].TB[1][node] + E->sphere.cap[j].TB[2][node] + E->sphere.cap[j].TB[3][node]);
            	break;
            } /* end switch */

            depth = E->sphere.ro - E->sx[j][3][node];

            switch (type) {
            case 0:  /* no match, next node */
            	break;
            default:
		dist=fabs(E->flag_depth2[nodeg]);
		if(E->flag_depth2[nodeg]*1000.0/6371.0>=-0.1 && E->flag_depth2[nodeg]<=0.0 && dist*1000.0/6371.0<=0.1)
			asm_depth=0.003;
		else if (dist<0.06)
                        asm_depth=0.003;
                else if (dist<=0.1)
                        asm_depth=0.003+(1.0-0.003)*(dist-0.06)/(0.1-0.06);
                else
                        asm_depth=1.0;
                assim_depth=asm_depth*E->control.lith_age_depth;
                if(depth <= assim_depth) {
                    daf = 0.5*depth/assim_depth;
                    E->T[j][node] = daf*E->T[j][node] + (1.0-daf)*assimilate_new_temp;
                }
		/* add initial slabs */
		/*if(intage<E->trench_visit_age && E->monitor.solution_cycles>0) {
        	flag_depth2=E->flag_depth2[nodeg]*1000.0/6371.0;
		rad = E->sx[j][3][node];
        	if(flag_depth2<0.0 && flag_depth2>-0.03 && rad>=0.97) {
        	        temp1=0.0009-flag_depth2*flag_depth2;
        	        temp2=0.026*0.026-flag_depth2*flag_depth2;
        	        temp3=0.0001-flag_depth2*flag_depth2;
        	        temp4=(rad-0.97)*(rad-0.97);
        	        if(temp4<=temp1 && temp4>=temp3) {
        	            dist2=0.03-sqrt(temp4+flag_depth2*flag_depth2);
        	            temp=dist2*0.5/sqrt(50.0/E->data.scalet);
        	            E->T[m][node]=E->control.lith_age_mantle_temp * erf(temp);
        	        }
        	   }
		}*/
            } /* end switch */
	    /*rad = E->sx[j][3][node];
	    if(intage<E->trench_visit_age && intage==120) {
		dist2 = pow(theta-125.0/180.0*3.1415926,2.0)/pow(15.0/180.0*3.1415926,2.0) + \
			pow(phi-220.0/180.0*3.1415926,2.0)/pow(20.0/180.0*3.1415926,2.0) + \
			pow(rad-0.95,2.0)/pow(0.03,2.0);
		if(dist2<=1.0)
		    E->T[j][node] = E->T[j][node] + 0.2 * cos(pow(dist2, 2.0)*3.1415926/2.0);
	    }*/

          } /* next node */

return;
}
