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
#include "element_definitions.h"
#include "global_defs.h"
#include <math.h>

#include "lith_age.h"

/* ========================================== */

static void horizontal_bc();
static void velocity_apply_periodic_bcs();
static void temperature_apply_periodic_bcs();
static void velocity_refl_vert_bc();
static void temperature_refl_vert_bc();
void velocity_apply_internal_bcs(struct All_variables*);
void velocity_apply_internal_bcs_uniform(struct All_variables*);
void construct_id(struct All_variables*);
float find_age_in_MY();

/* ========================================== */

void regional_velocity_boundary_conditions(E)
     struct All_variables *E;
{
  void velocity_imp_vert_bc();
  void read_velocity_boundary_from_file();
  void velocity_apply_internal_bcs();
  void renew_top_velocity_boundary();
  void apply_side_sbc(); 

  int node,d,j,noz,lv;

  for(lv=E->mesh.gridmax;lv>=E->mesh.gridmin;lv--)
    for (j=1;j<=E->sphere.caps_per_proc;j++)     {
      noz = E->lmesh.NOZ[lv];

      if(E->mesh.topvbc == 0) {
	horizontal_bc(E,E->sphere.cap[j].VB,noz,1,0.0,VBX,0,lv,j);
	horizontal_bc(E,E->sphere.cap[j].VB,noz,3,0.0,VBZ,1,lv,j);
	horizontal_bc(E,E->sphere.cap[j].VB,noz,2,0.0,VBY,0,lv,j);
	horizontal_bc(E,E->sphere.cap[j].VB,noz,1,E->control.VBXtopval,SBX,1,lv,j);
	horizontal_bc(E,E->sphere.cap[j].VB,noz,3,0.0,SBZ,0,lv,j);
	horizontal_bc(E,E->sphere.cap[j].VB,noz,2,E->control.VBYtopval,SBY,1,lv,j);
	}
      else if(E->mesh.topvbc == 1) {
        horizontal_bc(E,E->sphere.cap[j].VB,noz,1,E->control.VBXtopval,VBX,1,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,noz,3,0.0,VBZ,1,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,noz,2,E->control.VBYtopval,VBY,1,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,noz,1,0.0,SBX,0,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,noz,3,0.0,SBZ,0,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,noz,2,0.0,SBY,0,lv,j);

	if(E->control.vbcs_file)   {
	  read_velocity_boundary_from_file(E);   /* read in the velocity boundary condition from file */
	}
      }
      else if(E->mesh.topvbc == 2) {
	/* This extra BC is for a open top */
        horizontal_bc(E,E->sphere.cap[j].VB,noz,1,0.0,VBX,0,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,noz,3,0.0,VBZ,0,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,noz,2,0.0,VBY,0,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,noz,1,E->control.VBXtopval,SBX,1,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,noz,3,0.0,SBZ,1,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,noz,2,E->control.VBYtopval,SBY,1,lv,j);
        }
      else if(E->mesh.topvbc == 3) {
        /* This applies a hybrid internal VBC */
        horizontal_bc(E,E->sphere.cap[j].VB,noz,3,0.0,VBZ,1,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,noz,3,0.0,SBZ,0,lv,j);

        if(E->control.vbcs_file)   {
          read_velocity_boundary_from_file(E);   /* read in the velocity boundary condition from file */
          velocity_apply_internal_bcs(E);
          construct_id(E);
        }
      }
      else if(E->mesh.topvbc == 4) {
        /* This applies a uniform internal VBC */
        horizontal_bc(E,E->sphere.cap[j].VB,noz,1,E->control.VBXtopval,VBX,1,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,noz,3,0.0,VBZ,1,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,noz,2,E->control.VBYtopval,VBY,1,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,noz,1,0.0,SBX,0,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,noz,3,0.0,SBZ,0,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,noz,2,0.0,SBY,0,lv,j);

        if(E->control.vbcs_file)   {
          read_velocity_boundary_from_file(E);   /* read in the velocity boundary condition from file */
          velocity_apply_internal_bcs_uniform(E);
          construct_id(E);
        }
      }


      if(E->mesh.botvbc == 0) {
        horizontal_bc(E,E->sphere.cap[j].VB,1,1,0.0,VBX,0,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,1,3,0.0,VBZ,1,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,1,2,0.0,VBY,0,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,1,1,E->control.VBXbotval,SBX,1,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,1,3,0.0,SBZ,0,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,1,2,E->control.VBYbotval,SBY,1,lv,j);
        }
      else if(E->mesh.botvbc == 1) {
        horizontal_bc(E,E->sphere.cap[j].VB,1,1,E->control.VBXbotval,VBX,1,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,1,3,0.0,VBZ,1,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,1,2,E->control.VBYbotval,VBY,1,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,1,1,0.0,SBX,0,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,1,3,0.0,SBZ,0,lv,j);
        horizontal_bc(E,E->sphere.cap[j].VB,1,2,0.0,SBY,0,lv,j);
        }
      }    /* end for j and lv */

      velocity_refl_vert_bc(E);

      if(E->control.side_sbcs)
	apply_side_sbc(E);

      if(E->control.verbose) {
	for (j=1;j<=E->sphere.caps_per_proc;j++)
	  for (node=1;node<=E->lmesh.nno;node++)
	    fprintf(E->fp_out,"m=%d VB== %d %g %g %g flag %u %u %u\n",j,node,E->sphere.cap[j].VB[1][node],E->sphere.cap[j].VB[2][node],E->sphere.cap[j].VB[3][node],E->node[j][node]&VBX,E->node[j][node]&VBY,E->node[j][node]&VBZ);
	fflush(E->fp_out);
      }
      /* If any imposed internal velocity structure it goes here */


   return;
}


/* ========================================== */

void regional_temperature_boundary_conditions(E)
     struct All_variables *E;
{
  void temperature_imposed_vert_bcs();
  void temperature_lith_adj();
  void temperatures_conform_bcs();
  int j,lev,noz;

  lev = E->mesh.levmax;


     temperature_refl_vert_bc(E);

  for (j=1;j<=E->sphere.caps_per_proc;j++)    {
    noz = E->lmesh.noz;
    if(E->mesh.toptbc == 1)    {
      horizontal_bc(E,E->sphere.cap[j].TB,noz,3,E->control.TBCtopval,TBZ,1,lev,j);
      horizontal_bc(E,E->sphere.cap[j].TB,noz,3,E->control.TBCtopval,FBZ,0,lev,j);
      }
    else   {
      horizontal_bc(E,E->sphere.cap[j].TB,noz,3,E->control.TBCtopval,TBZ,0,lev,j);
      horizontal_bc(E,E->sphere.cap[j].TB,noz,3,E->control.TBCtopval,FBZ,1,lev,j);
      }

    if(E->mesh.bottbc == 1)    {
      horizontal_bc(E,E->sphere.cap[j].TB,1,3,E->control.TBCbotval,TBZ,1,lev,j);
      horizontal_bc(E,E->sphere.cap[j].TB,1,3,E->control.TBCbotval,FBZ,0,lev,j);
      }
    else        {
      horizontal_bc(E,E->sphere.cap[j].TB,1,3,E->control.TBCbotval,TBZ,0,lev,j);
      horizontal_bc(E,E->sphere.cap[j].TB,1,3,E->control.TBCbotval,FBZ,1,lev,j);
      }

    if((E->control.temperature_bound_adj==1) || (E->control.lith_age_time==1))  {
/* set the regions in which to use lithosphere files to determine temperature
   note that this is called if the lithosphere age in inputted every time step
   OR it is only maintained in the boundary regions */
      lith_age_temperature_bound_adj(E,lev);
    }

    }     /* end for j */

   temperatures_conform_bcs(E);
   E->temperatures_conform_bcs = temperatures_conform_bcs;

   return; }

/* ========================================== */

static void velocity_refl_vert_bc(E)
     struct All_variables *E;
{
  int m,i,j,ii,jj;
  int node1,node2;
  int level,nox,noy,noz;
  const int dims=E->mesh.nsd;

 /*  for two YOZ planes   */


  if (E->parallel.me_loc[1]==0 || E->parallel.me_loc[1]==E->parallel.nprocx-1)
   for (m=1;m<=E->sphere.caps_per_proc;m++)
    for(j=1;j<=E->lmesh.noy;j++)
      for(i=1;i<=E->lmesh.noz;i++)  {
        node1 = i + (j-1)*E->lmesh.noz*E->lmesh.nox;
        node2 = node1 + (E->lmesh.nox-1)*E->lmesh.noz;

        ii = i + E->lmesh.nzs - 1;
        if (E->parallel.me_loc[1]==0 )  {
           E->sphere.cap[m].VB[1][node1] = 0.0;
           if((ii != 1) && (ii != E->mesh.noz))
              E->sphere.cap[m].VB[3][node1] = 0.0;
               }
        if (E->parallel.me_loc[1]==E->parallel.nprocx-1)  {
           E->sphere.cap[m].VB[1][node2] = 0.0;
           if((ii != 1) && (ii != E->mesh.noz))
              E->sphere.cap[m].VB[3][node2] = 0.0;
           }
        }      /* end loop for i and j */

/*  for two XOZ  planes  */


    if (E->parallel.me_loc[2]==0)
     for (m=1;m<=E->sphere.caps_per_proc;m++)
      for(j=1;j<=E->lmesh.nox;j++)
        for(i=1;i<=E->lmesh.noz;i++)       {
          node1 = i + (j-1)*E->lmesh.noz;
          ii = i + E->lmesh.nzs - 1;

          E->sphere.cap[m].VB[2][node1] = 0.0;
          if((ii != 1) && (ii != E->mesh.noz))
            E->sphere.cap[m].VB[3][node1] = 0.0;
          }    /* end of loop i & j */

    if (E->parallel.me_loc[2]==E->parallel.nprocy-1)
     for (m=1;m<=E->sphere.caps_per_proc;m++)
      for(j=1;j<=E->lmesh.nox;j++)
        for(i=1;i<=E->lmesh.noz;i++)       {
          node2 = (E->lmesh.noy-1)*E->lmesh.noz*E->lmesh.nox + i + (j-1)*E->lmesh.noz;
          ii = i + E->lmesh.nzs - 1;

          E->sphere.cap[m].VB[2][node2] = 0.0;
          if((ii != 1) && (ii != E->mesh.noz))
            E->sphere.cap[m].VB[3][node2] = 0.0;
          }    /* end of loop i & j */


  /* all vbc's apply at all levels  */
  for(level=E->mesh.levmax;level>=E->mesh.levmin;level--) {

    if ( (E->control.CONJ_GRAD && level==E->mesh.levmax) ||E->control.NMULTIGRID)  {
    noz = E->lmesh.NOZ[level] ;
    noy = E->lmesh.NOY[level] ;
    nox = E->lmesh.NOX[level] ;

     for (m=1;m<=E->sphere.caps_per_proc;m++)  {
       if (E->parallel.me_loc[1]==0 || E->parallel.me_loc[1]==E->parallel.nprocx-1) {
         for(j=1;j<=noy;j++)
          for(i=1;i<=noz;i++) {
          node1 = i + (j-1)*noz*nox;
          node2 = node1 + (nox-1)*noz;
          ii = i + E->lmesh.NZS[level] - 1;
          if (E->parallel.me_loc[1]==0 )  {
            E->NODE[level][m][node1] = E->NODE[level][m][node1] | VBX;
            E->NODE[level][m][node1] = E->NODE[level][m][node1] & (~SBX);
            if((ii!=1) && (ii!=E->mesh.NOZ[level])) {
               E->NODE[level][m][node1] = E->NODE[level][m][node1] & (~VBY);
               E->NODE[level][m][node1] = E->NODE[level][m][node1] | SBY;
               E->NODE[level][m][node1] = E->NODE[level][m][node1] & (~VBZ);
               E->NODE[level][m][node1] = E->NODE[level][m][node1] | SBZ;
               }
            }
          if (E->parallel.me_loc[1]==E->parallel.nprocx-1)  {
            E->NODE[level][m][node2] = E->NODE[level][m][node2] | VBX;
            E->NODE[level][m][node2] = E->NODE[level][m][node2] & (~SBX);
            if((ii!=1) && (ii!=E->mesh.NOZ[level])) {
              E->NODE[level][m][node2] = E->NODE[level][m][node2] & (~VBY);
              E->NODE[level][m][node2] = E->NODE[level][m][node2] | SBY;
              E->NODE[level][m][node2] = E->NODE[level][m][node2] & (~VBZ);
              E->NODE[level][m][node2] = E->NODE[level][m][node2] | SBZ;
                  }
            }
          }   /* end for loop i & j */

         }


      if (E->parallel.me_loc[2]==0)
        for(j=1;j<=nox;j++)
          for(i=1;i<=noz;i++) {
            node1 = i + (j-1)*noz;
            ii = i + E->lmesh.NZS[level] - 1;
            jj = j + E->lmesh.NXS[level] - 1;

            E->NODE[level][m][node1] = E->NODE[level][m][node1] | VBY;
            E->NODE[level][m][node1] = E->NODE[level][m][node1] & (~SBY);
            if((ii!= 1) && (ii != E->mesh.NOZ[level]))  {
                E->NODE[level][m][node1] = E->NODE[level][m][node1] & (~VBZ);
                E->NODE[level][m][node1] = E->NODE[level][m][node1] | SBZ;
                }
            if((jj!=1) && (jj!=E->mesh.NOX[level]) && (ii!=1) && (ii!=E->mesh.NOZ[level])){
                E->NODE[level][m][node1] = E->NODE[level][m][node1] & (~VBX);
                E->NODE[level][m][node1] = E->NODE[level][m][node1] | SBX;
                }
                }    /* end for loop i & j  */

      if (E->parallel.me_loc[2]==E->parallel.nprocy-1)
        for(j=1;j<=nox;j++)
          for(i=1;i<=noz;i++)       {
            node2 = (noy-1)*noz*nox + i + (j-1)*noz;
            ii = i + E->lmesh.NZS[level] - 1;
            jj = j + E->lmesh.NXS[level] - 1;
            E->NODE[level][m][node2] = E->NODE[level][m][node2] | VBY;
            E->NODE[level][m][node2] = E->NODE[level][m][node2] & (~SBY);
            if((ii!= 1) && (ii != E->mesh.NOZ[level]))  {
                E->NODE[level][m][node2] = E->NODE[level][m][node2] & (~VBZ);
                E->NODE[level][m][node2] = E->NODE[level][m][node2] | SBZ;
                }
            if((jj!=1) && (jj!=E->mesh.NOX[level]) && (ii!=1) && (ii!=E->mesh.NOZ[level])){
                E->NODE[level][m][node2] = E->NODE[level][m][node2] & (~VBX);
                E->NODE[level][m][node2] = E->NODE[level][m][node2] | SBX;
                }
            }

       }       /* end for m  */
       }
       }       /*  end for loop level  */

  return;
}

static void temperature_refl_vert_bc(E)
     struct All_variables *E;
{
  int i,j,m;
  int node1,node2;
  const int dims=E->mesh.nsd;

 /* Temps and bc-values  at top level only */

   if (E->parallel.me_loc[1]==0 || E->parallel.me_loc[1]==E->parallel.nprocx-1)
    for(m=1;m<=E->sphere.caps_per_proc;m++)
    for(j=1;j<=E->lmesh.noy;j++)
      for(i=1;i<=E->lmesh.noz;i++) {
        node1 = i + (j-1)*E->lmesh.noz*E->lmesh.nox;
        node2 = node1 + (E->lmesh.nox-1)*E->lmesh.noz;
        if (E->parallel.me_loc[1]==0 )                   {
          E->node[m][node1] = E->node[m][node1] & (~TBX);
          E->node[m][node1] = E->node[m][node1] | FBX;
          E->sphere.cap[m].TB[1][node1] = 0.0;
              }
        if (E->parallel.me_loc[1]==E->parallel.nprocx-1)   {
          E->node[m][node2] = E->node[m][node2] & (~TBX);
          E->node[m][node2] = E->node[m][node2] | FBX;
          E->sphere.cap[m].TB[1][node2] = 0.0;
              }
        }       /* end for loop i & j */

    if (E->parallel.me_loc[2]==0)
     for(m=1;m<=E->sphere.caps_per_proc;m++)
      for(j=1;j<=E->lmesh.nox;j++)
        for(i=1;i<=E->lmesh.noz;i++) {
          node1 = i + (j-1)*E->lmesh.noz;
          E->node[m][node1] = E->node[m][node1] & (~TBY);
              E->node[m][node1] = E->node[m][node1] | FBY;
              E->sphere.cap[m].TB[2][node1] = 0.0;
              }

    if (E->parallel.me_loc[2]==E->parallel.nprocy-1)
     for(m=1;m<=E->sphere.caps_per_proc;m++)
      for(j=1;j<=E->lmesh.nox;j++)
        for(i=1;i<=E->lmesh.noz;i++) {
          node2 = i +(j-1)*E->lmesh.noz + (E->lmesh.noy-1)*E->lmesh.noz*E->lmesh.nox;
          E->node[m][node2] = E->node[m][node2] & (~TBY);
          E->node[m][node2] = E->node[m][node2] | FBY;
          E->sphere.cap[m].TB[3][node2] = 0.0;
          }    /* end loop for i and j */

  return;
}


/*  =========================================================  */


static void horizontal_bc(E,BC,ROW,dirn,value,mask,onoff,level,m)
     struct All_variables *E;
     float *BC[];
     int ROW;
     int dirn;
     float value;
     unsigned int mask;
     char onoff;
     int level,m;

{
  int i,j,node,rowl;

    /* safety feature */
  if(dirn > E->mesh.nsd)
     return;

  if (ROW==1)
      rowl = 1;
  else
      rowl = E->lmesh.NOZ[level];

  if ( (ROW==1 && E->parallel.me_loc[3]==0) ||
       (ROW==E->lmesh.NOZ[level] && E->parallel.me_loc[3]==E->parallel.nprocz-1) ) {

    /* turn bc marker to zero */
    if (onoff == 0)          {
      for(j=1;j<=E->lmesh.NOY[level];j++)
    	for(i=1;i<=E->lmesh.NOX[level];i++)     {
    	  node = rowl+(i-1)*E->lmesh.NOZ[level]+(j-1)*E->lmesh.NOX[level]*E->lmesh.NOZ[level];
    	  E->NODE[level][m][node] = E->NODE[level][m][node] & (~ mask);
    	  }        /* end for loop i & j */
      }

    /* turn bc marker to one */
    else        {
      for(j=1;j<=E->lmesh.NOY[level];j++)
        for(i=1;i<=E->lmesh.NOX[level];i++)       {
    	  node = rowl+(i-1)*E->lmesh.NOZ[level]+(j-1)*E->lmesh.NOX[level]*E->lmesh.NOZ[level];
    	  E->NODE[level][m][node] = E->NODE[level][m][node] | (mask);

    	  if(level==E->mesh.levmax)   /* NB */
    	    BC[dirn][node] = value;
    	  }     /* end for loop i & j */
      }

    }             /* end for if ROW */

  return;
}


static void velocity_apply_periodic_bcs(E)
    struct All_variables *E;
{
  int n1,n2,level;
  int i,j,ii,jj;
  const int dims=E->mesh.nsd;

  fprintf(E->fp,"Periodic boundary conditions\n");

  return;
  }

static void temperature_apply_periodic_bcs(E)
    struct All_variables *E;
{
 const int dims=E->mesh.nsd;

 fprintf(E->fp,"pERIodic temperature boundary conditions\n");

  return;
  }


/* Lijun Liu defined internal vbc */
/* still not working perfect */
void velocity_apply_internal_bcs11(E)
    struct All_variables *E;
{
 int j,ii,jj,kk,el,ie,lv,flag1;
 int start_lev,start_node,node,nodeg,node_coarse,node_fine;
 float theta,phi,r;
 const int dims=E->mesh.nsd;
 
 
 for (j=1;j<=E->sphere.caps_per_proc;j++)
    for(kk=1;kk<=E->lmesh.noy;kk++)
      for(jj=1;jj<=E->lmesh.nox;jj++)
        for(ii=E->lmesh.noz;ii>=1;ii--){
          if(E->parallel.me_loc[3]==E->parallel.nprocz-1) {
            node = ii + (jj-1)*E->lmesh.noz + (kk-1)*E->lmesh.noz*E->lmesh.nox;
            nodeg = E->lmesh.nxs-1+jj+(E->lmesh.nys+kk-2)*E->mesh.nox;
            theta = E->sx[j][1][node];
            phi = E->sx[j][2][node];
            r = E->sx[j][3][node];
            
            /* if(theta<E->control.theta_max-0.4 && theta>E->control.theta_min+0.3 && phi<E->control.fi_max-0.5 && phi>E->control.fi_min+0.3) { */
               /* if(r>0.997) { */
            flag1 = 0;
            for(lv=E->mesh.gridmax;lv>=E->mesh.gridmin;lv--) {
               if(lv==E->mesh.gridmax) {
                   /*if(ii==E->lmesh.noz) {
                        E->NODE[lv][j][node] = E->NODE[lv][j][node] & (~VBX);
                        E->NODE[lv][j][node] = E->NODE[lv][j][node] & (~VBY);
                        E->NODE[lv][j][node] = E->NODE[lv][j][node] | SBX;
                        E->NODE[lv][j][node] = E->NODE[lv][j][node] | SBY;
                        start_lev = lv;
                        start_node = node;
                        flag1 = 1;
                    }*/
                    if(ii==E->lmesh.noz-5) { /* z-increment by a power of 2 plus 1 */
                        E->NODE[lv][j][node] = E->NODE[lv][j][node] | VBX;
                        E->NODE[lv][j][node] = E->NODE[lv][j][node] | VBY;
                        E->NODE[lv][j][node] = E->NODE[lv][j][node] & (~SBX);
                        E->sphere.cap[j].VB[1][node] = E->sphere.cap[j].VB[1][node-ii+E->lmesh.noz];
                        E->NODE[lv][j][node] = E->NODE[lv][j][node] & (~SBY);
                        E->sphere.cap[j].VB[2][node] = E->sphere.cap[j].VB[2][node-ii+E->lmesh.noz];
                        start_lev = lv;
                        start_node = node;
                        flag1 = 1;
                    }
               }
               else if (flag1==1){
                  flag1 = 0;
                  for(el=1;el<=E->lmesh.NEL[lv];el++)
                      for(ie=1;ie<=enodes[dims];ie++)       {
                         node_coarse = E->IEN[lv][j][el].node[ie];
                         node_fine=E->IEN[start_lev][j][E->EL[lv][j][el].sub[ie]].node[ie];
                         if(node_fine==start_node && flag1==0) {
                            E->NODE[lv][j][node_coarse] = E->NODE[start_lev][j][start_node];
                            start_lev = lv;
                            start_node = node_coarse;
                            flag1 = 1;
                          }
                      } /* end of for(ie) */
               } /* end of else if */
            } /* end of for (lv) */
          } /* end of if */
        } /* end of ii */

  return;
  }


void velocity_apply_internal_bcs_uniform(E)
    struct All_variables *E;
{
 int i,j,k,m,lv;
 int noxl,noyl,nozl,nodel,nodeg;
 float theta,phi,r;
 float fi,v_trench,age;
 FILE *fp1;

 for (m=1;m<=E->sphere.caps_per_proc;m++)
    for(lv=E->mesh.gridmax;lv>=E->mesh.gridmin;lv--) {
       noxl = E->lmesh.NOX[lv];
       noyl = E->lmesh.NOY[lv];
       nozl = E->lmesh.NOZ[lv];

         for(i=1;i<=noyl;i++)
            for(j=1;j<=noxl;j++)
               for(k=1;k<=nozl;k++) {
                   if(lv==E->mesh.gridmax) {
                      nodeg=E->lmesh.nxs-1+j+(E->lmesh.nys+i-2)*E->mesh.nox;
                   }
                   nodel =  k + (j-1) * nozl + (i-1)*noxl*nozl;
                   theta = E->SX[lv][m][1][nodel];
                   phi = E->SX[lv][m][2][nodel];
                   r = E->SX[lv][m][3][nodel];

                   if(r>0.993 && r<1.0) {
                      E->NODE[lv][m][nodel] = E->NODE[lv][m][nodel] | VBX;;
                      E->NODE[lv][m][nodel] = E->NODE[lv][m][nodel] | VBY;
		      E->NODE[lv][m][nodel] = E->NODE[lv][m][nodel] & (~VBZ);
                      E->NODE[lv][m][nodel] = E->NODE[lv][m][nodel] & (~SBX);
                      E->NODE[lv][m][nodel] = E->NODE[lv][m][nodel] & (~SBY);
		      E->NODE[lv][m][nodel] = E->NODE[lv][m][nodel] | SBZ;
                      if(lv==E->mesh.gridmax) {
                         E->sphere.cap[m].VB[1][nodel] = E->velo_1[nodeg];
                         E->sphere.cap[m].VB[2][nodel] = E->velo_2[nodeg];
			 E->sphere.cap[m].VB[2][nodel] = 0.0;
                      } /* gridmax loop */
                   } /* if(r) loop */
                   /*else {
                      E->NODE[lv][m][nodel] = E->NODE[lv][m][nodel] & (~VBX);
                      E->NODE[lv][m][nodel] = E->NODE[lv][m][nodel] & (~VBY);
		      E->NODE[lv][m][nodel] = E->NODE[lv][m][nodel] & (~VBZ);
                      E->NODE[lv][m][nodel] = E->NODE[lv][m][nodel] | SBX;
                      E->NODE[lv][m][nodel] = E->NODE[lv][m][nodel] | SBY;
		      E->NODE[lv][m][nodel] = E->NODE[lv][m][nodel] | SBZ;
                   }*/

               } /* k loop */
  } /* lv loop */

  return;
}


void velocity_apply_internal_bcs(E)
    struct All_variables *E;
{
 int i,j,k,m,lv;
 int noxl,noyl,nozl,nodel,nodeg;
 float theta,phi,r;
 float fi,v_trench,age;
 FILE *fp1;

  if (E->control.lith_age_time) {
    age=find_age_in_MY(E);

    fp1=fopen("velo_trench","r");
    if(fp1==NULL)
        v_trench = -2.0;
    else {
        fscanf(fp1,"%f",&(v_trench));
        fclose(fp1);
    }

    if (E->parallel.me == 0) {
        fprintf(stderr,"INSIDE lith_age_temperature_bound_adj\n");
    }

    //fprintf(stderr,"start_age is %f\n",E->control.start_age);
    fi = 0.5 + v_trench*(3.14159/(11.1*180.0))*(E->control.start_age-age);
  } /* end of if loop */



 for (m=1;m<=E->sphere.caps_per_proc;m++) 
    for(lv=E->mesh.gridmax;lv>=E->mesh.gridmin;lv--) {
       noxl = E->lmesh.NOX[lv];
       noyl = E->lmesh.NOY[lv];
       nozl = E->lmesh.NOZ[lv];

       //if(E->parallel.me_loc[3]==E->parallel.nprocz-1) { 
         for(i=1;i<=noyl;i++)
            for(j=1;j<=noxl;j++) 
               for(k=1;k<=nozl;k++) {
                   if(lv==E->mesh.gridmax) {
                      nodeg=E->lmesh.nxs-1+j+(E->lmesh.nys+i-2)*E->mesh.nox;
                   } 
                   nodel =  k + (j-1) * nozl + (i-1)*noxl*nozl;
                   theta = E->SX[lv][m][1][nodel];
                   phi = E->SX[lv][m][2][nodel];
                   r = E->SX[lv][m][3][nodel];
        
                   if(r>0.993 && r<1.1 && (phi<fi-0.05 || phi>fi+0.05)) { // && phi>0.4 && phi<0.5) {
                   //if(r>0.997 && phi>fi+0.0) {
                      E->NODE[lv][m][nodel] = E->NODE[lv][m][nodel] | VBX;;
                      E->NODE[lv][m][nodel] = E->NODE[lv][m][nodel] | VBY;
                      E->NODE[lv][m][nodel] = E->NODE[lv][m][nodel] & (~SBX);
                      E->NODE[lv][m][nodel] = E->NODE[lv][m][nodel] & (~SBY);
                      if(lv==E->mesh.gridmax) {
                         E->sphere.cap[m].VB[1][nodel] = E->velo_1[nodeg];
                         E->sphere.cap[m].VB[2][nodel] = E->velo_2[nodeg];
			 //E->sphere.cap[m].VB[3][nodel] = 0.0;
                      } /* gridmax loop */
                   } /* if(r) loop */ 
		   else {
		      E->NODE[lv][m][nodel] = E->NODE[lv][m][nodel] & (~VBX);
                      E->NODE[lv][m][nodel] = E->NODE[lv][m][nodel] & (~VBY);
                      E->NODE[lv][m][nodel] = E->NODE[lv][m][nodel] | SBX;
                      E->NODE[lv][m][nodel] = E->NODE[lv][m][nodel] | SBY;
		   }

               } /* k loop */ 
       //} /* E->parallel.me loop */     
  } /* lv loop */

  return;
}

/* version */
/* $Id: Regional_boundary_conditions.c 4437 2006-08-26 00:49:06Z tan2 $ */

/* End of file  */
