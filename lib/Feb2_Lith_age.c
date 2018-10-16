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


void lith_age_input(struct All_variables *E)
{
  int m = E->parallel.me;

  E->control.lith_age = 0;
  E->control.lith_age_time = 0;
  E->control.temperature_bound_adj = 0;

  input_int("lith_age",&(E->control.lith_age),"0",m);
  input_float("mantle_temp",&(E->control.lith_age_mantle_temp),"1.0",m);

  if (E->control.lith_age) {
    input_int("lith_age_time",&(E->control.lith_age_time),"0",m);
    input_string("lith_age_file",E->control.lith_age_file,"",m);
    input_float("lith_age_depth",&(E->control.lith_age_depth),"0.0471",m);

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
	E->age_t[node]=E->age_t[node]*E->data.scalet;
      }
    fclose(fp1);
  } /* end E->control.lith_age_time == false */
}


void lith_age_construct_tic(struct All_variables *E)
{
  int i, j, k, m, node, nodeg;
  int nox, noy, noz, gnox, gnoy, gnoz;
  double r1, temp, temp1;
  float age,theta,phi,assim_depth,asm_depth;
  void temperatures_conform_bcs();

  noy=E->lmesh.noy;
  nox=E->lmesh.nox;
  noz=E->lmesh.noz;

  gnox=E->mesh.nox;
  gnoy=E->mesh.noy;
  gnoz=E->mesh.noz;

  for(m=1;m<=E->sphere.caps_per_proc;m++)
    for(i=1;i<=noy;i++)
      for(j=1;j<=nox;j++)
	for(k=1;k<=noz;k++)  {
	  nodeg=E->lmesh.nxs-1+j+(E->lmesh.nys+i-2)*gnox;
	  node=k+(j-1)*noz+(i-1)*nox*noz;
	  r1=E->sx[m][3][node];
	//  assim_depth = /*E->asm_depth[nodeg]**/E->control.lith_age_depth;
	  if (asm_depth<0.06)
              asm_depth=0.003;
          else if (asm_depth<=0.3)
              asm_depth=0.003+(asm_depth-0.06)/(0.3007-0.06);
          else
              asm_depth=1.0;
	  assim_depth = asm_depth*E->control.lith_age_depth;
	  E->T[m][node] = E->control.lith_age_mantle_temp;
	  if( r1 >= E->sphere.ro-assim_depth)
	    { /* if closer than (lith_age_depth) from top */
	      /* zero age surface has mantle temperature */
	      if(E->age_t[nodeg]>0.0 && E->age_t[nodeg]<180.0/E->data.scalet) //Lijun: define continent to be 0 age
                    temp = (E->sphere.ro-r1) *0.5 /sqrt(E->age_t[nodeg])/5;
                else
                    temp = 10.0;
	      //temp = 10.0;
	      temp1 = E->control.lith_age_mantle_temp * erf(temp);
	//      if(E->age_t[nodeg]>180.0/E->data.scalet && r1>=E->sphere.ro-1.0e-4) temp1 = 0.0;

	      if((E->node[m][node] & TBZ) == TBZ) 
		E->T[m][node] = temp1;

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
  float tt1,tt2,ttt1,ttt2,ttt3,fff2,fff3,*PB1[4],*PB2[4],fi,fi_1,fi_2,gap,assim_depth,asm_depth;
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
		asm_depth=fabs(E->asm_depth[nodeg]);
	        if (asm_depth<0.06)
			asm_depth=0.003;
                else if (asm_depth<=0.3)
                        asm_depth=0.003+(asm_depth-0.06)/(0.3007-0.06);
                else
                        asm_depth=1.0;
                assim_depth=asm_depth*E->control.lith_age_depth;

                if(E->sx[j][3][node]>=E->sphere.ro-assim_depth) {
                // if closer than (lith_age_depth) from top
                    if(E->age_t[nodeg]>0.0/E->data.scalet && E->age_t[nodeg]<180.0/E->data.scalet
                      && theta<E->control.theta_max-0.2 && theta>E->control.theta_min+0.2
      	              && phi<E->control.fi_max-0.2 && phi>E->control.fi_min+0.2) { //Lijun ed.
		    
                        E->node[j][node]=E->node[j][node] | TBX;
                        E->node[j][node]=E->node[j][node] & (~FBX);
                        E->node[j][node]=E->node[j][node] | TBY;
                        E->node[j][node]=E->node[j][node] & (~FBY);
                        E->node[j][node]=E->node[j][node] | TBZ;
                        E->node[j][node]=E->node[j][node] & (~FBZ);
                    }
		    /*else if(E->age_t[nodeg]>180.0/E->data.scalet && E->sx[j][3][node]>=E->sphere.ro-1.0e-4) {
			E->node[j][node]=E->node[j][node] | TBX;
                        E->node[j][node]=E->node[j][node] & (~FBX);
                        E->node[j][node]=E->node[j][node] | TBY;
                        E->node[j][node]=E->node[j][node] & (~FBY);
                        E->node[j][node]=E->node[j][node] | TBZ;
                        E->node[j][node]=E->node[j][node] & (~FBZ);
		    }*/
		    else {
			E->node[j][node]=E->node[j][node] | FBX;
                        E->node[j][node]=E->node[j][node] & (~TBX);
                        E->node[j][node]=E->node[j][node] | FBY;
                        E->node[j][node]=E->node[j][node] & (~TBY);
                        E->node[j][node]=E->node[j][node] | FBZ;
                        E->node[j][node]=E->node[j][node] & (~TBZ);
		    }
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
  float r1,t1,f1,t0,temp,temp1;
  float depth,assim_depth,asm_depth;
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
    //if (E->parallel.me == 0) fprintf(stderr,"INSIDE lith_age_conform_tbc\n");
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
		if(E->age_t[nodeg]>0.0 && E->age_t[nodeg]<180.0/E->data.scalet) //Lijun define continent to be 0 age with T_sfc=1.0
		    temp = (E->sphere.ro-r1) *0.5 /sqrt(E->age_t[nodeg])/3;
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
		if(E->age_t[nodeg]>0.0 && E->age_t[nodeg]<180.0/E->data.scalet)
		    temp = (E->sphere.ro-r1) *0.5 /sqrt(E->age_t[nodeg])/3;
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

	    if(fabs(r1-E->sphere.ro)>=e_4 && fabs(r1-E->sphere.ri)>=e_4)  { // if NOT right on the boundary 
	      asm_depth=fabs(E->asm_depth[nodeg]);
	      if (asm_depth<0.06)
                        asm_depth=0.003;
              else if (asm_depth<=0.3)
                        asm_depth=0.003+(asm_depth-0.06)/(0.3007-0.06);
              else
                        asm_depth=1.0;
	      assim_depth=asm_depth*E->control.lith_age_depth;
	      if(  E->sx[m][3][node]>=E->sphere.ro-assim_depth ) {
		// if closer than (lith_age_depth) from top 

                depth=E->sphere.ro - E->sx[m][3][node];

		// set a new age from the file 
//		temp = (E->sphere.ro-r1) *0.5 /sqrt(E->age_t[nodeg]);
		if(E->age_t[nodeg]>0.0 && E->age_t[nodeg]<180.0/E->data.scalet)
                    temp = (E->sphere.ro-r1) *0.5 /sqrt(E->age_t[nodeg])/5;
                else
                    temp = 10.0;
		t0 = E->control.lith_age_mantle_temp * erf(temp);
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
  float depth, daf, assimilate_new_temp,temp1,temp2,fi_1,fi_2,fi,*PB1[4],*PB2[4],asm_depth,assim_depth;
  float theta,phi,lat_max,lat_max1,lat_max2,lat_min,lat_min1,lat_min2,age,age1,age2;
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
	    age=find_age_in_MY(E);

            switch (type) {
            case 0:  /* no match, next node */
            	break;
            default:
		asm_depth=fabs(E->asm_depth[nodeg]);
		if (asm_depth<0.06)
                        asm_depth=0.003;
                else if (asm_depth<=0.3)
                        asm_depth=0.003+(asm_depth-0.06)/(0.3007-0.06);
                else
                        asm_depth=1.0;
                assim_depth=asm_depth*E->control.lith_age_depth;
                if(depth <= assim_depth) {
                    if(E->age_t[nodeg]>0.01/E->data.scalet && E->age_t[nodeg]<0.01/E->data.scalet/* && phi<fi-0.02*/) {
                        daf = 0.5*depth/assim_depth;
                        E->T[j][node] = daf*E->T[j][node] + (1.0-daf)*assimilate_new_temp;
                    }
                }
            } /* end switch */

          } /* next node */

return;
}

