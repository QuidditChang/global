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
/* Functions to assemble the element k matrices and the element f vector.
   Note that for the regular grid case the calculation of k becomes repetitive
   to the point of redundancy. */

#include <math.h>
#include "element_definitions.h"
#include "global_defs.h"
#include "material_properties.h"
#include "ala_block_vector.h"
#include "ala_coupled_operator.h"

void myerror(struct All_variables *,char *);



void add_force(struct All_variables *E, int e, double elt_f[24], int m)
{
  const int dims=E->mesh.nsd;
  const int ends=enodes[E->mesh.nsd];
  int a, a1, a2, a3, p, node;

  for(a=1;a<=ends;a++)          {
    node = E->ien[m][e].node[a];
    p=(a-1)*dims;
    a1=E->id[m][node].doff[1];
    E->F[m][a1] += elt_f[p];
    a2=E->id[m][node].doff[2];
    E->F[m][a2] += elt_f[p+1];
    a3=E->id[m][node].doff[3];
    E->F[m][a3] += elt_f[p+2];
  }
}



/* ================================================================
   Function to assemble the global  F vector.
                     +
   Function to get the global H vector (mixed method driving terms)
   ================================================================ */

void assemble_forces(E,penalty)
     struct All_variables *E;
     int penalty;
{
  double elt_f[24];
  int m,a,e,i;

  void get_buoyancy();
  void get_elt_f();
  void get_elt_tr();
  void strip_bcs_from_residual();

  const int neq=E->lmesh.neq;
  const int nel=E->lmesh.nel;
  const int lev=E->mesh.levmax;

  get_buoyancy(E,E->buoyancy);

  for(m=1;m<=E->sphere.caps_per_proc;m++)    {

    for(a=0;a<neq;a++)
      E->F[m][a] = 0.0;

    for (e=1;e<=nel;e++)  {
      get_elt_f(E,e,elt_f,1,m);
      add_force(E, e, elt_f, m);
    }

    /* for traction bc */
    for(i=1; i<=E->boundary.nel; i++) {
      e = E->boundary.element[m][i];

      for(a=0;a<24;a++) elt_f[a] = 0.0;
      for(a=SIDE_BEGIN; a<=SIDE_END; a++)
	get_elt_tr(E, i, a, elt_f, m);

      add_force(E, e, elt_f, m);
    }
  }       /* end for m */

  (E->solver.exchange_id_d)(E, E->F, lev);
  strip_bcs_from_residual(E,E->F,lev);
  return;
}


void assemble_forces_pseudo_surf(E,penalty)
     struct All_variables *E;
     int penalty;
{
  double elt_f[24];
  int m,a,e,i;

  void get_buoyancy();
  void get_elt_f();
  void get_elt_tr_pseudo_surf();
  void strip_bcs_from_residual();

  const int neq=E->lmesh.neq;
  const int nel=E->lmesh.nel;
  const int lev=E->mesh.levmax;

  get_buoyancy(E,E->buoyancy);

  for(m=1;m<=E->sphere.caps_per_proc;m++)    {

    for(a=0;a<neq;a++)
      E->F[m][a] = 0.0;

    for (e=1;e<=nel;e++)  {
      get_elt_f(E,e,elt_f,1,m);
      add_force(E, e, elt_f, m);
    }

    /* for traction bc */
    for(i=1; i<=E->boundary.nel; i++) {
      e = E->boundary.element[m][i];

      for(a=0;a<24;a++) elt_f[a] = 0.0;
      for(a=SIDE_BEGIN; a<=SIDE_END; a++)
	get_elt_tr_pseudo_surf(E, i, a, elt_f, m);

      add_force(E, e, elt_f, m);
    }
  }       /* end for m */

  (E->solver.exchange_id_d)(E, E->F, lev);
  strip_bcs_from_residual(E,E->F,lev);
  return;
}


/*==============================================================
  Function to supply the element strain-displacement matrix Ba at velocity
  quadrature points, which is used to compute element stiffness matrix
  ==============================================================  */

void get_ba(struct Shape_function *N, struct Shape_function_dx *GNx,
       struct CC *cc, struct CCX *ccx, double rtf[4][9],
       int dims, double ba[9][9][4][7])
{
    int k, a, n;
    const int vpts = VPOINTS3D;
    const int ends = ENODES3D;

    double ra[9], isi[9], ct[9];
    double gnx0, gnx1, gnx2, shp, cc1, cc2, cc3;

    for(k=1;k<=vpts;k++) {
    ra[k] = rtf[3][k];
    isi[k] = 1.0 / sin(rtf[1][k]);
    ct[k] = cos(rtf[1][k]) * isi[k];
    }

    for(a=1;a<=ends;a++)
        for(k=1;k<=vpts;k++) {
            gnx0 = GNx->vpt[GNVXINDEX(0,a,k)];
            gnx1 = GNx->vpt[GNVXINDEX(1,a,k)];
            gnx2 = GNx->vpt[GNVXINDEX(2,a,k)];
            shp  = N->vpt[GNVINDEX(a,k)];
            for(n=1;n<=dims;n++) {
                cc1 = cc->vpt[BVINDEX(1,n,a,k)];
                cc2 = cc->vpt[BVINDEX(2,n,a,k)];
                cc3 = cc->vpt[BVINDEX(3,n,a,k)];

        ba[a][k][n][1] = ( gnx0 * cc1
                                   + shp * ccx->vpt[BVXINDEX(1,n,1,a,k)]
                                   + shp * cc3 ) * ra[k];

        ba[a][k][n][2] = ( shp * cc1 * ct[k]
                                   + shp * cc3
                                   + ( gnx1 * cc2
                                       + shp * ccx->vpt[BVXINDEX(2,n,2,a,k)] )
                                   * isi[k] ) * ra[k];

        ba[a][k][n][3] = gnx2 * cc3;

        ba[a][k][n][4] = ( gnx0 * cc2
                                   + shp * ccx->vpt[BVXINDEX(2,n,1,a,k)]
                                   - shp * cc2 * ct[k]
                                   + ( gnx1 * cc1
                                       + shp * ccx->vpt[BVXINDEX(1,n,2,a,k)] )
                                   * isi[k] ) * ra[k];

        ba[a][k][n][5] = gnx2 * cc1
                               + ( gnx0 * cc3
                                   + shp * ( ccx->vpt[BVXINDEX(3,n,1,a,k)]
                                             - cc1 ) ) * ra[k];

        ba[a][k][n][6] = gnx2 * cc2
                               - ra[k] * shp * cc2
                               + ( gnx1 * cc3
                                   + shp * ccx->vpt[BVXINDEX(3,n,2,a,k)] )
                               * isi[k] * ra[k];
        }
        }

    return;
}


/*==============================================================
  Function to supply the element strain-displacement matrix Ba at pressure
  quadrature points, which is used to compute strain rate
  ==============================================================  */

void get_ba_p(struct Shape_function *N, struct Shape_function_dx *GNx,
              struct CC *cc, struct CCX *ccx, double rtf[4][9],
              int dims, double ba[9][9][4][7])
{
    int k, a, n;
    const int ppts = PPOINTS3D;
    const int ends = ENODES3D;

    double ra[9], isi[9], ct[9];
    double gnx0, gnx1, gnx2, shp, cc1, cc2, cc3;

    for(k=1;k<=ppts;k++) {
        ra[k] = rtf[3][k];
        isi[k] = 1.0 / sin(rtf[1][k]);
        ct[k] = cos(rtf[1][k]) * isi[k];
    }

    for(k=1;k<=ppts;k++)
        for(a=1;a<=ends;a++) {
            gnx0 = GNx->ppt[GNPXINDEX(0,a,k)];
            gnx1 = GNx->ppt[GNPXINDEX(1,a,k)];
            gnx2 = GNx->ppt[GNPXINDEX(2,a,k)];
            shp  = N->ppt[GNPINDEX(a,k)];
            for(n=1;n<=dims;n++) {
                cc1 = cc->ppt[BPINDEX(1,n,a,k)];
                cc2 = cc->ppt[BPINDEX(2,n,a,k)];
                cc3 = cc->ppt[BPINDEX(3,n,a,k)];

        ba[a][k][n][1] = ( gnx0 * cc1
                           + shp * ccx->ppt[BPXINDEX(1,n,1,a,k)]
                           + shp * cc3 ) * ra[k];

        ba[a][k][n][2] = ( shp * cc1 * ct[k]
                           + shp * cc3
                           + ( gnx1 * cc2
                               + shp * ccx->ppt[BPXINDEX(2,n,2,a,k)] )
                           * isi[k] ) * ra[k];

        ba[a][k][n][3] = gnx2 * cc3;

        ba[a][k][n][4] = ( gnx0 * cc2
                           + shp * ccx->ppt[BPXINDEX(2,n,1,a,k)]
                           - shp * cc2 * ct[k]
                           + ( gnx1 * cc1
                               + shp * ccx->ppt[BPXINDEX(1,n,2,a,k)] )
                           * isi[k] ) * ra[k];

        ba[a][k][n][5] = gnx2 * cc1
                       + ( gnx0 * cc3
                           + shp * ( ccx->ppt[BPXINDEX(3,n,1,a,k)]
                                     - cc1 ) ) * ra[k];

        ba[a][k][n][6] = gnx2 * cc2
                       - ra[k] * shp * cc2
                       + ( gnx1 * cc3
                           + shp * ccx->ppt[BPXINDEX(3,n,2,a,k)] )
                       * isi[k] * ra[k];
            }
        }
    return;
}



/*==============================================================
  Function to supply the element k matrix for a given element e.
  ==============================================================  */

void get_elt_k(E,el,elt_k,lev,m,iconv)
     struct All_variables *E;
     int el,m;
     double elt_k[24*24];
     int lev, iconv;
{
    double bdbmu[4][4];
    int pn,qn,ad,bd;

    int a,b,i,j,i1,j1,k;
    double rtf[4][9],W[9];

    const double two = 2.0;
    const double two_thirds = 2.0/3.0;

    void get_global_shape_fn();
    void construct_c3x3matrix_el();
    struct Shape_function GN;
    struct Shape_function_dx GNx;
    struct Shape_function_dA dOmega;

    double ba[9][9][4][7]; /* integration points,node,3x6 matrix */

    const int nn=loc_mat_size[E->mesh.nsd];
    const int vpts = VPOINTS3D;
    const int ends = ENODES3D;
    const int dims=E->mesh.nsd;
    const int sphere_key = 1;

    get_global_shape_fn(E,el,&GN,&GNx,&dOmega,0,sphere_key,rtf,lev,m);

    if (iconv || (el-1)%E->lmesh.ELZ[lev]==0)
      construct_c3x3matrix_el(E,el,&E->element_Cc,&E->element_Ccx,lev,m,0);

    /* Note N[a].gauss_pt[n] is the value of shape fn a at the nth gaussian
       quadrature point. Nx[d] is the derivative wrt x[d]. */

    for(k=1;k<=vpts;k++) {
      W[k]=g_point[k].weight[dims-1]*dOmega.vpt[k]*E->EVI[lev][m][(el-1)*vpts+k];
    }

    get_ba(&(E->N), &GNx, &E->element_Cc, &E->element_Ccx,
           rtf, E->mesh.nsd, ba);

  for(a=1;a<=ends;a++)
    for(b=a;b<=ends;b++)   {
      bdbmu[1][1]=bdbmu[1][2]=bdbmu[1][3]=
      bdbmu[2][1]=bdbmu[2][2]=bdbmu[2][3]=
      bdbmu[3][1]=bdbmu[3][2]=bdbmu[3][3]=0.0;

      for(i=1;i<=dims;i++)
        for(j=1;j<=dims;j++)
          for(k=1;k<=VPOINTS3D;k++)
              bdbmu[i][j] += W[k] * ( two * ( ba[a][k][i][1]*ba[b][k][j][1] +
                                              ba[a][k][i][2]*ba[b][k][j][2] +
                                              ba[a][k][i][3]*ba[b][k][j][3] ) +
                                      ba[a][k][i][4]*ba[b][k][j][4] +
                                      ba[a][k][i][5]*ba[b][k][j][5] +
                                      ba[a][k][i][6]*ba[b][k][j][6] );

      if(E->control.inv_gruneisen != 0)
        for(i=1;i<=dims;i++)
          for(j=1;j<=dims;j++)
            for(k=1;k<=VPOINTS3D;k++)
                bdbmu[i][j] -= W[k] * two_thirds *
                    ( ba[a][k][i][1] + ba[a][k][i][2] + ba[a][k][i][3] ) *
                    ( ba[b][k][j][1] + ba[b][k][j][2] + ba[b][k][j][3] );


                /**/
      ad=dims*(a-1);
      bd=dims*(b-1);

      pn=ad*nn+bd;
      qn=bd*nn+ad;

      elt_k[pn       ] = bdbmu[1][1] ; /* above */
      elt_k[pn+1     ] = bdbmu[1][2] ;
      elt_k[pn+2     ] = bdbmu[1][3] ;
      elt_k[pn+nn    ] = bdbmu[2][1] ;
      elt_k[pn+nn+1  ] = bdbmu[2][2] ;
      elt_k[pn+nn+2  ] = bdbmu[2][3] ;
      elt_k[pn+2*nn  ] = bdbmu[3][1] ;
      elt_k[pn+2*nn+1] = bdbmu[3][2] ;
      elt_k[pn+2*nn+2] = bdbmu[3][3] ;

      elt_k[qn       ] = bdbmu[1][1] ; /* below diag */
      elt_k[qn+1     ] = bdbmu[2][1] ;
      elt_k[qn+2     ] = bdbmu[3][1] ;
      elt_k[qn+nn    ] = bdbmu[1][2] ;
      elt_k[qn+nn+1  ] = bdbmu[2][2] ;
      elt_k[qn+nn+2  ] = bdbmu[3][2] ;
      elt_k[qn+2*nn  ] = bdbmu[1][3] ;
      elt_k[qn+2*nn+1] = bdbmu[2][3] ;
      elt_k[qn+2*nn+2] = bdbmu[3][3] ;
                /**/

      } /*  Sum over all the a,b's to obtain full  elt_k matrix */

    return;
}


/* =============================================
   General calling function for del_squared:
   according to whether it should be element by
   element or node by node.
   ============================================= */

void assemble_del2_u(E,u,Au,level,strip_bcs)
     struct All_variables *E;
     double **u,**Au;
     int level;
     int strip_bcs;
{
  void e_assemble_del2_u();
  void n_assemble_del2_u();

  if(E->control.NMULTIGRID||E->control.NASSEMBLE)
    n_assemble_del2_u(E,u,Au,level,strip_bcs);
  else
    e_assemble_del2_u(E,u,Au,level,strip_bcs);

  return;
}

/* ======================================
   Assemble del_squared_u vector el by el
   ======================================   */

void e_assemble_del2_u(E,u,Au,level,strip_bcs)
  struct All_variables *E;
  double **u,**Au;
  int level;
  int strip_bcs;

{
  int  e,i,a,b,a1,a2,a3,ii,m,nodeb;
  void strip_bcs_from_residual();

  const int n=loc_mat_size[E->mesh.nsd];
  const int ends=enodes[E->mesh.nsd];
  const int dims=E->mesh.nsd;
  const int nel=E->lmesh.NEL[level];
  const int neq=E->lmesh.NEQ[level];

  for (m=1;m<=E->sphere.caps_per_proc;m++)   {
    for(i=0;i<neq;i++)
      Au[m][i] = 0.0;

    for(e=1;e<=nel;e++)   {
      for(a=1;a<=ends;a++) {
	ii = E->IEN[level][m][e].node[a];
	a1 = E->ID[level][m][ii].doff[1];
	a2 = E->ID[level][m][ii].doff[2];
	a3 = E->ID[level][m][ii].doff[3];
	for(b=1;b<=ends;b++) {
	        nodeb = E->IEN[level][m][e].node[b];
	        ii = (a*n+b)*dims-(dims*n+dims);
			/* i=1, j=1,2 */
	          /* i=1, j=1,2,3 */
		Au[m][a1] +=
		        E->elt_k[level][m][e].k[ii] *
			u[m][E->ID[level][m][nodeb].doff[1]]
		      + E->elt_k[level][m][e].k[ii+1] *
			u[m][E->ID[level][m][nodeb].doff[2]]
		      + E->elt_k[level][m][e].k[ii+2] *
			u[m][E->ID[level][m][nodeb].doff[3]];
		/* i=2, j=1,2,3 */
		Au[m][a2] +=
		        E->elt_k[level][m][e].k[ii+n] *
			u[m][E->ID[level][m][nodeb].doff[1]]
		      + E->elt_k[level][m][e].k[ii+n+1] *
			u[m][E->ID[level][m][nodeb].doff[2]]
		      + E->elt_k[level][m][e].k[ii+n+2] *
			u[m][E->ID[level][m][nodeb].doff[3]];
		/* i=3, j=1,2,3 */
		Au[m][a3] +=
		        E->elt_k[level][m][e].k[ii+n+n] *
			u[m][E->ID[level][m][nodeb].doff[1]]
		      + E->elt_k[level][m][e].k[ii+n+n+1] *
			u[m][E->ID[level][m][nodeb].doff[2]]
		      + E->elt_k[level][m][e].k[ii+n+n+2] *
			u[m][E->ID[level][m][nodeb].doff[3]];

 	    }         /* end for loop b */
        }             /* end for loop a */

       }          /* end for e */
     }         /* end for m  */

    (E->solver.exchange_id_d)(E, Au, level);

  if(strip_bcs)
     strip_bcs_from_residual(E,Au,level);

  return; }


/* ======================================================
   Assemble Au using stored, nodal coefficients.
   ====================================================== */

void n_assemble_del2_u(E,u,Au,level,strip_bcs)
     struct All_variables *E;
     double **u,**Au;
     int level;
     int strip_bcs;
{
    int m, e,i;
    int eqn1,eqn2,eqn3;

    double UU,U1,U2,U3;
    void strip_bcs_from_residual();

    int *C;
    higher_precision *B1,*B2,*B3;

    const int neq=E->lmesh.NEQ[level];
    const int nno=E->lmesh.NNO[level];
    const int dims=E->mesh.nsd;
    const int max_eqn = dims*14;


  for (m=1;m<=E->sphere.caps_per_proc;m++)  {

     for(e=0;e<=neq;e++)
	Au[m][e]=0.0;

     u[m][neq] = 0.0;

     for(e=1;e<=nno;e++)     {

       eqn1=E->ID[level][m][e].doff[1];
       eqn2=E->ID[level][m][e].doff[2];
       eqn3=E->ID[level][m][e].doff[3];

       U1 = u[m][eqn1];
       U2 = u[m][eqn2];
       U3 = u[m][eqn3];

       C=E->Node_map[level][m] + (e-1)*max_eqn;
       B1=E->Eqn_k1[level][m]+(e-1)*max_eqn;
       B2=E->Eqn_k2[level][m]+(e-1)*max_eqn;
       B3=E->Eqn_k3[level][m]+(e-1)*max_eqn;

       for(i=3;i<max_eqn;i++)  {
	  UU = u[m][C[i]];
  	  Au[m][eqn1] += B1[i]*UU;
  	  Au[m][eqn2] += B2[i]*UU;
  	  Au[m][eqn3] += B3[i]*UU;
       }
       for(i=0;i<max_eqn;i++)
          Au[m][C[i]] += B1[i]*U1+B2[i]*U2+B3[i]*U3;

       }     /* end for e */
     }     /* end for m */

     (E->solver.exchange_id_d)(E, Au, level);

    if (strip_bcs)
	strip_bcs_from_residual(E,Au,level);

    return;
}


void build_diagonal_of_K(E,el,elt_k,level,m)
     struct All_variables *E;
     int level,el,m;
     double elt_k[24*24];

{
    int a,a1,a2,p,node;

    const int n=loc_mat_size[E->mesh.nsd];
    const int dims=E->mesh.nsd;
    const int ends=enodes[E->mesh.nsd];

    for(a=1;a<=ends;a++) {
	    node=E->IEN[level][m][el].node[a];
	    /* dirn 1 */
	    a1 = E->ID[level][m][node].doff[1];
	    p=(a-1)*dims;
	    E->BI[level][m][a1] += elt_k[p*n+p];

	    /* dirn 2 */
	    a2 = E->ID[level][m][node].doff[2];
	    p=(a-1)*dims+1;
	    E->BI[level][m][a2] += elt_k[p*n+p];

	    /* dirn 3 */
	    a1 = E->ID[level][m][node].doff[3];
	    p=(a-1)*dims+2;
	    E->BI[level][m][a1] += elt_k[p*n+p];
            }

  return;
}

void build_diagonal_of_Ahat(E)
    struct All_variables *E;
{
    double assemble_dAhatp_entry();

    double BU;
    int m,e,npno,level;

  /* Initialize every level and cap before honoring the off switch.  The old
     early return initialized only the first level/cap and left the remaining
     entries stale even though all strict-ALA Krylov paths read BPI. */
  for(level=E->mesh.gridmin;level<=E->mesh.gridmax;level++)
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
      npno=E->lmesh.NPNO[level];
      for(e=1;e<=npno;e++)
        E->BPI[level][m][e]=1.0;
    }

  if(!E->control.precondition)
    return;

  for(level=E->mesh.gridmin;level<=E->mesh.gridmax;level++)
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
      npno=E->lmesh.NPNO[level];
      for(e=1;e<=npno;e++) {
        BU=assemble_dAhatp_entry(E,e,level,m);
        if(BU != 0.0)
          E->BPI[level][m][e]=1.0/BU;
        else
          E->BPI[level][m][e]=1.0;
      }
    }

  return;
}


void build_diagonal_of_leng_zhong_Ahat(E)
    struct All_variables *E;
{
    double assemble_leng_zhong_dAhatp_entry();
    double diagonal;
    int e, level, m, npno;

    for(level=E->mesh.gridmin;level<=E->mesh.gridmax;level++)
      for(m=1;m<=E->sphere.caps_per_proc;m++) {
        npno=E->lmesh.NPNO[level];
        for(e=1;e<=npno;e++)
          E->LZ_BPI[level][m][e]=1.0;
      }

    if(!E->control.ala_leng_zhong_2008 || !E->control.precondition)
      return;

    for(level=E->mesh.gridmin;level<=E->mesh.gridmax;level++)
      for(m=1;m<=E->sphere.caps_per_proc;m++) {
        npno=E->lmesh.NPNO[level];
        for(e=1;e<=npno;e++) {
          diagonal=assemble_leng_zhong_dAhatp_entry(E,e,level,m);
          if(!isfinite(diagonal) || diagonal<=0.0)
            myerror(E, "Leng-Zhong Stage 2 D-only BPI is not positive finite");
          E->LZ_BPI[level][m][e]=1.0/diagonal;
        }
      }

    return;
}


double assemble_Ahatp_jacobi_entry(E,e1,e2,level,m)
     struct All_variables *E;
     int e1,e2,level,m;
{
    int a,b,d,p1,p2,node1,node2,eqn;
    double g1,g2,value;
    const int ends=enodes[E->mesh.nsd];
    const int dims=E->mesh.nsd;

    value=0.0;
    for(a=1;a<=ends;a++) {
      node1=E->IEN[level][m][e1].node[a];
      for(b=1;b<=ends;b++) {
        node2=E->IEN[level][m][e2].node[b];
        if(node1 != node2)
          continue;
        for(d=1;d<=dims;d++) {
          p1=(a-1)*dims+d-1;
          p2=(b-1)*dims+d-1;
          g1=E->elt_del[level][m][e1].g[p1][0];
          g2=E->elt_del[level][m][e2].g[p2][0];
          if(E->control.ala_pressure_buoyancy) {
            g1 += E->elt_c[level][m][e1].c[p1][0];
            g2 += E->elt_c[level][m][e2].c[p2][0];
          }
          eqn=E->ID[level][m][node1].doff[d];
          value += g1*E->ALA_velocity_BI[level][m][eqn]*g2;
        }
      }
    }
    return(value);
}


void build_radial_line_Ahat_preconditioner(E)
     struct All_variables *E;
{
    int m,col,k,e,previous,level,elz,ncolumns;
    int local_fallback,global_fallback,local_columns,global_columns;
    double diagonal,offdiag,pivot,previous_pivot,line_max;
    double local_min,local_max,global_min,global_max;
    const double pivot_fraction=1.0e-12;

    level=E->mesh.levmax;
    elz=E->lmesh.ELZ[level];
    ncolumns=E->lmesh.ELX[level]*E->lmesh.ELY[level];
    local_columns=ncolumns*E->sphere.caps_per_proc;
    local_fallback=0;
    local_min=1.0e300;
    local_max=0.0;

    for(m=1;m<=E->sphere.caps_per_proc;m++)
      for(col=0;col<ncolumns;col++) {
        E->ALA_BPI_line_valid[level][m][col+1]=1;
        line_max=0.0;
        for(k=0;k<elz;k++) {
          e=col*elz+k+1;
          diagonal=1.0/E->BPI[level][m][e];
          E->ALA_BPI_line_diag[level][m][e]=diagonal;
          E->ALA_BPI_line_lower[level][m][e]=0.0;
          line_max=max(line_max,diagonal);
        }

        previous_pivot=0.0;
        for(k=0;k<elz;k++) {
          e=col*elz+k+1;
          diagonal=E->ALA_BPI_line_diag[level][m][e];
          if(k==0)
            pivot=diagonal;
          else {
            previous=e-1;
            offdiag=assemble_Ahatp_jacobi_entry(
                E,previous,e,level,m);
            E->ALA_BPI_line_lower[level][m][e]=
                offdiag/previous_pivot;
            pivot=diagonal-
                E->ALA_BPI_line_lower[level][m][e]*offdiag;
          }
          if(!isfinite(pivot) || pivot<=pivot_fraction*line_max) {
            E->ALA_BPI_line_valid[level][m][col+1]=0;
            local_fallback++;
            break;
          }
          E->ALA_BPI_line_diag[level][m][e]=pivot;
          previous_pivot=pivot;
          local_min=min(local_min,pivot);
          local_max=max(local_max,pivot);
        }
      }

    MPI_Allreduce(&local_fallback,&global_fallback,1,MPI_INT,MPI_SUM,
                  E->parallel.world);
    MPI_Allreduce(&local_columns,&global_columns,1,MPI_INT,MPI_SUM,
                  E->parallel.world);
    MPI_Allreduce(&local_min,&global_min,1,MPI_DOUBLE,MPI_MIN,
                  E->parallel.world);
    MPI_Allreduce(&local_max,&global_max,1,MPI_DOUBLE,MPI_MAX,
                  E->parallel.world);
    if(E->parallel.me==0) {
      fprintf(E->fp,
              "ALA radial-line preconditioner: pivot_range=(%e,%e) "
              "fallback_columns=%d/%d\n",
              global_min,global_max,global_fallback,global_columns);
      fprintf(stderr,
              "ALA radial-line preconditioner: pivot_range=(%e,%e) "
              "fallback_columns=%d/%d\n",
              global_min,global_max,global_fallback,global_columns);
    }
}


/* =====================================================
   Assemble grad(rho_ref*ez)*V element by element.
   Note that the storage is not zero'd before assembling.
   =====================================================  */

void assemble_c_u(struct All_variables *E,
                  double **U, double **result, int level)
{
    int e,j1,j2,j3,p,a,b,m;

    const int nel = E->lmesh.NEL[level];
    const int ends = enodes[E->mesh.nsd];
    const int dims = E->mesh.nsd;
    const int npno = E->lmesh.NPNO[level];

    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(a=1;a<=ends;a++) {
            p = (a-1)*dims;
            for(e=1;e<=nel;e++) {
                b = E->IEN[level][m][e].node[a];
                j1= E->ID[level][m][b].doff[1];
                j2= E->ID[level][m][b].doff[2];
                j3= E->ID[level][m][b].doff[3];
                result[m][e] += E->elt_c[level][m][e].c[p  ][0] * U[m][j1]
                              + E->elt_c[level][m][e].c[p+1][0] * U[m][j2]
                              + E->elt_c[level][m][e].c[p+2][0] * U[m][j3];
            }
        }

    return;
}



/* =====================================================
   Assemble div(rho_ref*V) = div(V) + grad(rho_ref*ez)*V
   element by element
   =====================================================  */

void assemble_div_rho_u(struct All_variables *E,
                        double **U, double **result, int level)
{
    void assemble_div_u();
    assemble_div_u(E, U, result, level);
    assemble_c_u(E, U, result, level);

    return;
}


/* Assemble the strict-ALA continuity operator with an explicitly selected
   finest-grid element beta field. This is diagnostic-only unless the selected
   field is also E->refstate.ala_beta, in which case it is algebraically
   identical to assemble_div_rho_u(). */
void assemble_div_rho_u_with_beta(struct All_variables *E,
                                  double **U, double **result, int level,
                                  double *fine_beta)
{
    int m,e,a,b,p,j1,j2,j3,nz,fine_first,fine_last,fine_nz;
    double beta,active_beta,dr,dr_total,scale;
    const int nel=E->lmesh.NEL[level];
    const int ends=enodes[E->mesh.nsd];
    const int dims=E->mesh.nsd;
    void assemble_div_u();

    assemble_div_u(E,U,result,level);
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(e=1;e<=nel;e++) {
            nz=((e-1)%E->lmesh.ELZ[level])+1;
            fine_first=((nz-1)*E->lmesh.elz)/E->lmesh.ELZ[level]+1;
            fine_last=(nz*E->lmesh.elz)/E->lmesh.ELZ[level];
            beta=0.0;
            active_beta=0.0;
            dr_total=0.0;
            for(fine_nz=fine_first;fine_nz<=fine_last;fine_nz++) {
                dr=E->sx[1][3][fine_nz+1]-E->sx[1][3][fine_nz];
                beta += fine_beta[fine_nz]*dr;
                active_beta += E->refstate.ala_beta[fine_nz]*dr;
                dr_total += dr;
            }
            beta /= dr_total;
            active_beta /= dr_total;
            if(!isfinite(beta) || beta<=0.0 ||
               !isfinite(active_beta) || active_beta<=0.0)
                myerror(E,"Invalid beta in ALA causal residual assembly");
            scale=beta/active_beta;
            for(a=1;a<=ends;a++) {
                p=(a-1)*dims;
                b=E->IEN[level][m][e].node[a];
                j1=E->ID[level][m][b].doff[1];
                j2=E->ID[level][m][b].doff[2];
                j3=E->ID[level][m][b].doff[3];
                result[m][e] += scale*(
                    E->elt_c[level][m][e].c[p][0]*U[m][j1]
                   +E->elt_c[level][m][e].c[p+1][0]*U[m][j2]
                   +E->elt_c[level][m][e].c[p+2][0]*U[m][j3]);
            }
        }
}


/* ==========================================
   Assemble a div_u vector element by element
   ==========================================  */

void assemble_div_u(struct All_variables *E,
                    double **U, double **divU, int level)
{
    int e,j1,j2,j3,p,a,b,m;

    const int nel=E->lmesh.NEL[level];
    const int ends=enodes[E->mesh.nsd];
    const int dims=E->mesh.nsd;
    const int npno=E->lmesh.NPNO[level];

  for(m=1;m<=E->sphere.caps_per_proc;m++)
    for(e=1;e<=npno;e++)
	divU[m][e] = 0.0;

  for(m=1;m<=E->sphere.caps_per_proc;m++)
       for(a=1;a<=ends;a++)   {
	  p = (a-1)*dims;
          for(e=1;e<=nel;e++) {
	    b = E->IEN[level][m][e].node[a];
	    j1= E->ID[level][m][b].doff[1];
            j2= E->ID[level][m][b].doff[2];
	    j3= E->ID[level][m][b].doff[3];
	    divU[m][e] += E->elt_del[level][m][e].g[p  ][0] * U[m][j1]
	                + E->elt_del[level][m][e].g[p+1][0] * U[m][j2]
	                + E->elt_del[level][m][e].g[p+2][0] * U[m][j3];
	    }
	 }

    return;
}


/* ==========================================
   Assemble a grad_P vector element by element
   ==========================================  */

void assemble_grad_p(E,P,gradP,lev)
     struct All_variables *E;
     double **P,**gradP;
     int lev;

{
  int m,e,i,j1,j2,j3,p,a,b,nel,neq;
  void strip_bcs_from_residual();

  const int ends=enodes[E->mesh.nsd];
  const int dims=E->mesh.nsd;

  for(m=1;m<=E->sphere.caps_per_proc;m++)  {

    nel=E->lmesh.NEL[lev];
    neq=E->lmesh.NEQ[lev];

    for(i=0;i<neq;i++)
      gradP[m][i] = 0.0;

    for(e=1;e<=nel;e++) {

	if(0.0==P[m][e])
	    continue;

	for(a=1;a<=ends;a++)       {
	     p = (a-1)*dims;
	     b = E->IEN[lev][m][e].node[a];
	     j1= E->ID[lev][m][b].doff[1];
	     j2= E->ID[lev][m][b].doff[2];
	     j3= E->ID[lev][m][b].doff[3];
		        /*for(b=0;b<ploc_mat_size[E->mesh.nsd];b++)  */
             gradP[m][j1] += E->elt_del[lev][m][e].g[p  ][0] * P[m][e];
             gradP[m][j2] += E->elt_del[lev][m][e].g[p+1][0] * P[m][e];
             gradP[m][j3] += E->elt_del[lev][m][e].g[p+2][0] * P[m][e];
	     }
        }       /* end for el */
     }       /* end for m */

  (E->solver.exchange_id_d)(E, gradP,  lev); /*  correct gradP   */


  strip_bcs_from_residual(E,gradP,lev);

return;
}


/* ==============================================================
   Assemble the strict-ALA pressure operator (D+C)^T * P.
   This is the exact discrete transpose of assemble_div_rho_u().
   ============================================================== */

static void assemble_grad_rho_p_local_terms(struct All_variables *E,
                                            double **P, double **gradP,
                                            int lev)
{
  int m,e,i,j1,j2,j3,p,a,b,nel,neq;

  const int ends=enodes[E->mesh.nsd];
  const int dims=E->mesh.nsd;

  for(m=1;m<=E->sphere.caps_per_proc;m++) {
    nel=E->lmesh.NEL[lev];
    neq=E->lmesh.NEQ[lev];

    for(i=0;i<neq;i++)
      gradP[m][i] = 0.0;

    for(e=1;e<=nel;e++) {
      if(0.0==P[m][e])
        continue;

      for(a=1;a<=ends;a++) {
        p = (a-1)*dims;
        b = E->IEN[lev][m][e].node[a];
        j1 = E->ID[lev][m][b].doff[1];
        j2 = E->ID[lev][m][b].doff[2];
        j3 = E->ID[lev][m][b].doff[3];
        gradP[m][j1] += (E->elt_del[lev][m][e].g[p][0]
                         + E->elt_c[lev][m][e].c[p][0]) * P[m][e];
        gradP[m][j2] += (E->elt_del[lev][m][e].g[p+1][0]
                         + E->elt_c[lev][m][e].c[p+1][0]) * P[m][e];
        gradP[m][j3] += (E->elt_del[lev][m][e].g[p+2][0]
                         + E->elt_c[lev][m][e].c[p+2][0]) * P[m][e];
      }
    }
  }
}


void assemble_grad_rho_p(struct All_variables *E,
                         double **P, double **gradP, int lev)
{
  void strip_bcs_from_residual();

  assemble_grad_rho_p_local_terms(E,P,gradP,lev);
  (E->solver.exchange_id_d)(E, gradP, lev);
  strip_bcs_from_residual(E,gradP,lev);

  return;
}


/* Exact diagnostic transpose of assemble_c_u(). */
void assemble_grad_c_p(struct All_variables *E,
                       double **P, double **gradP, int lev)
{
  int m,e,i,j1,j2,j3,p,a,b,nel,neq;
  const int ends=enodes[E->mesh.nsd];
  const int dims=E->mesh.nsd;
  void strip_bcs_from_residual();

  for(m=1;m<=E->sphere.caps_per_proc;m++) {
    nel=E->lmesh.NEL[lev];
    neq=E->lmesh.NEQ[lev];
    for(i=0;i<neq;i++)
      gradP[m][i]=0.0;
    for(e=1;e<=nel;e++) {
      if(P[m][e]==0.0)
        continue;
      for(a=1;a<=ends;a++) {
        p=(a-1)*dims;
        b=E->IEN[lev][m][e].node[a];
        j1=E->ID[lev][m][b].doff[1];
        j2=E->ID[lev][m][b].doff[2];
        j3=E->ID[lev][m][b].doff[3];
        gradP[m][j1] += E->elt_c[lev][m][e].c[p][0]*P[m][e];
        gradP[m][j2] += E->elt_c[lev][m][e].c[p+1][0]*P[m][e];
        gradP[m][j3] += E->elt_c[lev][m][e].c[p+2][0]*P[m][e];
      }
    }
  }
  (E->solver.exchange_id_d)(E,gradP,lev);
  strip_bcs_from_residual(E,gradP,lev);
}


double assemble_dAhatp_entry(E,e,level,m)
     struct All_variables *E;
     int e,level,m;

{
    int i,j,p,a,b,node,npno;
    void strip_bcs_from_residual();

    double gradP[81],divU,pressure_op;

    const int ends=enodes[E->mesh.nsd];
    const int dims=E->mesh.nsd;

    npno=E->lmesh.NPNO[level];

    for(i=0;i<81;i++)
	gradP[i] = 0.0;

    divU=0.0;

    for(a=1;a<=ends;a++) {
      p = (a-1)*dims;
      node = E->IEN[level][m][e].node[a];
      j=E->ID[level][m][node].doff[1];
      pressure_op = E->elt_del[level][m][e].g[p][0];
      if(E->control.ala_pressure_buoyancy)
        pressure_op += E->elt_c[level][m][e].c[p][0];
      gradP[p] += E->ALA_velocity_BI[level][m][j] * pressure_op;

      j=E->ID[level][m][node].doff[2];
      pressure_op = E->elt_del[level][m][e].g[p+1][0];
      if(E->control.ala_pressure_buoyancy)
        pressure_op += E->elt_c[level][m][e].c[p+1][0];
      gradP[p+1] += E->ALA_velocity_BI[level][m][j] * pressure_op;

      j=E->ID[level][m][node].doff[3];
      pressure_op = E->elt_del[level][m][e].g[p+2][0];
      if(E->control.ala_pressure_buoyancy)
        pressure_op += E->elt_c[level][m][e].c[p+2][0];
      gradP[p+2] += E->ALA_velocity_BI[level][m][j] * pressure_op;
      }


    /* calculate div U from the same thing .... */

    /* only need to run over nodes with non-zero grad P, i.e. the ones in
       the element accessed above, BUT it is only necessary to update the
       value in the original element, because the diagonal is all we use at
       the end ... */

    for(b=1;b<=ends;b++) {
      p = (b-1)*dims;
      pressure_op = E->elt_del[level][m][e].g[p][0];
      if(E->control.ala_pressure_buoyancy)
        pressure_op += E->elt_c[level][m][e].c[p][0];
      divU += pressure_op * gradP[p];

      pressure_op = E->elt_del[level][m][e].g[p+1][0];
      if(E->control.ala_pressure_buoyancy)
        pressure_op += E->elt_c[level][m][e].c[p+1][0];
      divU += pressure_op * gradP[p+1];

      pressure_op = E->elt_del[level][m][e].g[p+2][0];
      if(E->control.ala_pressure_buoyancy)
        pressure_op += E->elt_c[level][m][e].c[p+2][0];
      divU += pressure_op * gradP[p+2];
      }

return(divU);  }


double assemble_leng_zhong_dAhatp_entry(E,e,level,m)
     struct All_variables *E;
     int e,level,m;
{
    int a,b,j,node,p;
    double divU,gradP[81];
    const int ends=enodes[E->mesh.nsd];
    const int dims=E->mesh.nsd;

    for(p=0;p<81;p++)
      gradP[p]=0.0;

    for(a=1;a<=ends;a++) {
      p=(a-1)*dims;
      node=E->IEN[level][m][e].node[a];

      j=E->ID[level][m][node].doff[1];
      gradP[p]=E->ALA_velocity_BI[level][m][j]
          *E->elt_del[level][m][e].g[p][0];
      j=E->ID[level][m][node].doff[2];
      gradP[p+1]=E->ALA_velocity_BI[level][m][j]
          *E->elt_del[level][m][e].g[p+1][0];
      j=E->ID[level][m][node].doff[3];
      gradP[p+2]=E->ALA_velocity_BI[level][m][j]
          *E->elt_del[level][m][e].g[p+2][0];
    }

    divU=0.0;
    for(b=1;b<=ends;b++) {
      p=(b-1)*dims;
      divU += E->elt_del[level][m][e].g[p][0]*gradP[p];
      divU += E->elt_del[level][m][e].g[p+1][0]*gradP[p+1];
      divU += E->elt_del[level][m][e].g[p+2][0]*gradP[p+2];
    }

    return(divU);
}


/*==============================================================
  Function to supply the element c matrix for a given element e.
  ==============================================================  */

void get_elt_c(struct All_variables *E, int el,
               higher_precision elt_c[24][1], int lev, int m)
{
    void get_global_shape_fn();
    void construct_c3x3matrix_el();
    int p, a, i, j, nz, fine_first, fine_last, fine_nz;
    double temp, beta, dr, dr_total, rho_avg, x[4], rho[9];

    struct Shape_function GN;
    struct Shape_function_dx GNx;
    struct Shape_function_dA dOmega;
    double rtf[4][9];

    const int dims = E->mesh.nsd;
    const int ends = enodes[dims];
    const int sphere_key = 1;

    get_global_shape_fn(E,el,&GN,&GNx,&dOmega,2,sphere_key,rtf,lev,m);

    if ((el-1)%E->lmesh.ELZ[lev]==0)
        construct_c3x3matrix_el(E,el,&E->element_Cc,&E->element_Ccx,lev,m,1);

    temp = p_point[1].weight[dims-1] * dOmega.ppt[1];

    if(E->control.ala_pressure_buoyancy) {
        /* Strict ALA: restrict the single selected finest-grid element beta
         * by radial length.  For source=interval, these are the exact
         * serialized-density log secants.  On levmax this is identical to
         * the beta used by pressure buoyancy. */
        nz = ((el-1) % E->lmesh.ELZ[lev]) + 1;
        fine_first = ((nz-1) * E->lmesh.elz) / E->lmesh.ELZ[lev] + 1;
        fine_last = (nz * E->lmesh.elz) / E->lmesh.ELZ[lev];
        beta = 0.0;
        dr_total = 0.0;
        for(fine_nz=fine_first; fine_nz<=fine_last; fine_nz++) {
            dr = E->sx[1][3][fine_nz+1] - E->sx[1][3][fine_nz];
            beta += E->refstate.ala_beta[fine_nz] * dr;
            dr_total += dr;
        }
        beta /= dr_total;

        /* ala_beta=-d(ln rho)/dr, hence div(u)=ala_beta*u_r. */
        for(a=1;a<=ends;a++) {
            for(i=1;i<=dims;i++) {
                x[i] = E->N.ppt[GNPINDEX(a,1)]
                    * E->element_Cc.ppt[BPINDEX(3,i,a,1)];
            }
            p=dims*(a-1);
            elt_c[p  ][0] = x[1] * temp * beta;
            elt_c[p+1][0] = x[2] * temp * beta;
            elt_c[p+2][0] = x[3] * temp * beta;
        }
    }
    else if(E->refstate.choice == 1) {
        /* Explicit legacy analytic TALA benchmark policy. */
        beta = -E->control.disptn_number * E->control.inv_gruneisen;
        for(a=1;a<=ends;a++) {
            for(i=1;i<=dims;i++) {
                x[i] = E->N.ppt[GNPINDEX(a,1)]
                    * E->element_Cc.ppt[BPINDEX(3,i,a,1)];
            }
            p=dims*(a-1);
            elt_c[p  ][0] = -x[1] * temp * beta;
            elt_c[p+1][0] = -x[2] * temp * beta;
            elt_c[p+2][0] = -x[3] * temp * beta;
        }
    }
    else {
        /* Historical file-based TALA density-gradient assembly, unchanged. */
        for(a=1;a<=ends;a++) {
            j = E->IEN[lev][m][el].node[a];
            nz = (j - 1) % E->lmesh.noz + 1;
            rho[a] = E->refstate.rho[nz];
        }
        rho_avg = 0.0;
        for(a=1;a<=ends;a++)
            rho_avg += rho[a];
        rho_avg /= ends;
        for(a=1;a<=ends;a++) {
            for(i=1;i<=dims;i++) {
                x[i] = rho[a] * GNx.ppt[GNPXINDEX(2,a,1)]
                    * E->N.ppt[GNPINDEX(a,1)]
                    * E->element_Cc.ppt[BPINDEX(3,i,a,1)];
            }
            p=dims*(a-1);
            elt_c[p  ][0] = -x[1] * temp / rho_avg;
            elt_c[p+1][0] = -x[2] * temp / rho_avg;
            elt_c[p+2][0] = -x[3] * temp / rho_avg;
        }
    }

    return;
}


/*==============================================================
  Function to supply the element g matrix for a given element e.
  ==============================================================  */

void get_elt_g(E,el,elt_del,lev,m)
     struct All_variables *E;
     int el,m;
     higher_precision elt_del[24][1];
     int lev;

{  void get_global_shape_fn();
   void construct_c3x3matrix_el();
   int p,a,i;
   double ra,ct,si,x[4],rtf[4][9];
   double temp;

   struct Shape_function GN;
   struct Shape_function_dA dOmega;
   struct Shape_function_dx GNx;

   const int dims=E->mesh.nsd;
   const int ends=enodes[dims];
   const int sphere_key = 1;

   /* Special case, 4/8 node bilinear cartesian square/cube element -> 1 pressure point */

/*   es = (el-1)/E->lmesh.ELZ[lev]+1;  */

   if ((el-1)%E->lmesh.ELZ[lev]==0)
      construct_c3x3matrix_el(E,el,&E->element_Cc,&E->element_Ccx,lev,m,1);

   get_global_shape_fn(E,el,&GN,&GNx,&dOmega,2,sphere_key,rtf,lev,m);


   temp=p_point[1].weight[dims-1] * dOmega.ppt[1];

   ra = rtf[3][1];
   si = 1.0/sin(rtf[1][1]);
   ct = cos(rtf[1][1])*si;

   for(a=1;a<=ends;a++)      {
     for (i=1;i<=dims;i++)
       x[i]=GNx.ppt[GNPXINDEX(2,a,1)]*E->element_Cc.ppt[BPINDEX(3,i,a,1)]
        + 2.0*ra*E->N.ppt[GNPINDEX(a,1)]*E->element_Cc.ppt[BPINDEX(3,i,a,1)]
        + ra*(GNx.ppt[GNPXINDEX(0,a,1)]*E->element_Cc.ppt[BPINDEX(1,i,a,1)]
        +E->N.ppt[GNPINDEX(a,1)]*E->element_Ccx.ppt[BPXINDEX(1,i,1,a,1)]
        +ct*E->N.ppt[GNPINDEX(a,1)]*E->element_Cc.ppt[BPINDEX(1,i,a,1)]
        +si*(GNx.ppt[GNPXINDEX(1,a,1)]*E->element_Cc.ppt[BPINDEX(2,i,a,1)]
        +E->N.ppt[GNPINDEX(a,1)]*E->element_Ccx.ppt[BPXINDEX(2,i,2,a,1)]));

     p=dims*(a-1);
     elt_del[p  ][0] = -x[1] * temp;
     elt_del[p+1][0] = -x[2] * temp;
     elt_del[p+2][0] = -x[3] * temp;

      /* fprintf (E->fp,"B= %d %d %g %g %g %g %g\n",el,a,dOmega.ppt[1],GNx.ppt[GNPXINDEX(0,a,1)],GNx.ppt[GNPXINDEX(1,a,1)],elt_del[p][0],elt_del[p+1][0]);
      */
     }

return;
}


/*=================================================================
  Function to create the element force vector (allowing for velocity b.c.'s)
  ================================================================= */

void get_elt_f(E,el,elt_f,bcs,m)
     struct All_variables *E;
     int el,m;
     double elt_f[24];
     int bcs;

{

  int i,p,a,b,j,k,q,es;
  int got_elt_k,nodea,nodeb;
  unsigned int type;
  const unsigned int vbc_flag[] = {0, VBX, VBY, VBZ};

  double force[9],force_at_gs[9],elt_k[24*24];
  double rtf[4][9];

  void get_global_shape_fn();
  void construct_c3x3matrix_el();
  void get_ala_aug_k();

  struct Shape_function GN;
  struct Shape_function_dA dOmega;
  struct Shape_function_dx GNx;

  const int dims=E->mesh.nsd;
  const int n=loc_mat_size[dims];
  const int ends=enodes[dims];
  const int vpts=vpoints[dims];
  const int sphere_key=1;

  es = (el-1)/E->lmesh.elz + 1;

  get_global_shape_fn(E,el,&GN,&GNx,&dOmega,0,sphere_key,rtf,E->mesh.levmax,m);

  if ((el-1)%E->lmesh.elz==0)
      construct_c3x3matrix_el(E,el,&E->element_Cc,&E->element_Ccx,E->mesh.levmax,m,0);

  for(p=0;p<n;p++) elt_f[p] = 0.0;

  for(p=1;p<=ends;p++)
    force[p] = E->buoyancy[m][E->ien[m][el].node[p]];

  for(j=1;j<=vpts;j++)       {   /*compute force at each int point */
    force_at_gs[j] = 0.0;
    for(k=1;k<=ends;k++)
      force_at_gs[j] += force[k] * E->N.vpt[GNVINDEX(k,j)] ;
    }

  for(i=1;i<=dims;i++)  {
    for(a=1;a<=ends;a++)  {
      nodea=E->ien[m][el].node[a];
      p= dims*(a-1)+i-1;

      for(j=1;j<=vpts;j++)     /*compute sum(Na(j)*F(j)*det(j)) */
        elt_f[p] += force_at_gs[j] * E->N.vpt[GNVINDEX(a,j)]
           *dOmega.vpt[j]*g_point[j].weight[dims-1]
           *E->element_Cc.vpt[BVINDEX(3,i,a,j)];

      /* Strict ALA pressure buoyancy is -C^T*p.  Use the same
       * element matrix as continuity so the two operators are adjoints. */
      if(E->control.ala_pressure_buoyancy)
        elt_f[p] -= E->elt_c[E->mesh.levmax][m][el].c[p][0]
                    * E->P[m][el];

	  /* imposed velocity terms */

      if(bcs)  {
        got_elt_k = 0;
        for(j=1;j<=dims;j++) {
	  type=vbc_flag[j];
          for(b=1;b<=ends;b++) {
            nodeb=E->ien[m][el].node[b];
            if ((E->node[m][nodeb]&type)&&(E->sphere.cap[m].VB[j][nodeb]!=0.0)){
              if(!got_elt_k) {
                get_elt_k(E,el,elt_k,E->mesh.levmax,m,1);
                if(E->control.ala_augmented_lagrangian_gamma>0.0)
                  get_ala_aug_k(E,el,elt_k,E->mesh.levmax,m);
                got_elt_k = 1;
                }
              q = dims*(b-1)+j-1;
              if(p!=q) {
                elt_f[p] -= elt_k[p*n+q] * E->sphere.cap[m].VB[j][nodeb];
                }
              }
            }  /* end for b */
          }      /* end for j */
        }      /* end if for if bcs */

      }
    } /*  Complete the loops for a,i  	*/



  return;
}


/*=================================================================
  Function to create the element force vector due to stress b.c.
  ================================================================= */

void get_elt_tr(struct All_variables *E, int bel, int side, double elt_tr[24], int m)
{
	void construct_side_c3x3matrix_el();

	const int dims=E->mesh.nsd;
	const int ends1=enodes[dims-1];
	const int oned = onedvpoints[dims];

	struct CC Cc;
	struct CCX Ccx;

	const unsigned sbc_flag[4] = {0,SBX,SBY,SBZ};

	double traction[4][5],traction_at_gs[4][5], value, tmp;
	int j, b, p, k, a, nodea, d;
	int el = E->boundary.element[m][bel];
	int flagged;
	int found = 0;

	const float rho = E->data.density;
	const float g = E->data.grav_acc;
	const float R = 6371000.0;
	const float eta = E->data.ref_viscosity;
	const float kappa = E->data.kappa0;
	const float factor = 1.0e+00;
	int nodeas;

	if(E->control.side_sbcs)
		for(a=1;a<=ends1;a++)  {
			nodea = E->ien[m][el].node[ sidenodes[side][a] ];
			for(d=1;d<=dims;d++) {
				value = E->sbc.SB[m][side][d][ E->sbc.node[m][nodea] ];
				flagged = (E->node[m][nodea] & sbc_flag[d]) && (value);
				found |= flagged;
				traction[d][a] = ( flagged ? value : 0.0 );
			}
		}
	else {
		/* if side_sbcs is false, only apply sbc on top and bottom surfaces */
		if(side == SIDE_BOTTOM || side == SIDE_TOP) {
			for(a=1;a<=ends1;a++)  {
				nodea = E->ien[m][el].node[ sidenodes[side][a] ];
				for(d=1;d<=dims;d++) {
					value = E->sphere.cap[m].VB[d][nodea];
					flagged = (E->node[m][nodea] & sbc_flag[d]) && (value);
					found |= flagged;
					traction[d][a] = ( flagged ? value : 0.0 );
				}
			}
		}
	}

	/* skip the following computation if no sbc_flag is set
	   or value of sbcs are zero */
	if(!found) return;

	/* compute traction at each int point */
	construct_side_c3x3matrix_el(E,el,&Cc,&Ccx,
				     E->mesh.levmax,m,0,side);

	for(k=1;k<=oned;k++)
		for(d=1;d<=dims;d++) {
			traction_at_gs[d][k] = 0.0;
			for(j=1;j<=ends1;j++)
				traction_at_gs[d][k] += traction[d][j] * E->M.vpt[GMVINDEX(j,k)] ;
		}

	for(j=1;j<=ends1;j++) {
		a = sidenodes[side][j];
		for(d=1;d<=dims;d++) {
			p = dims*(a-1)+d-1;
			for(k=1;k<=oned;k++) {
				tmp = 0.0;
				for(b=1;b<=dims;b++)
					tmp += traction_at_gs[b][k] * Cc.vpt[BVINDEX(b,d,a,k)];

				elt_tr[p] += tmp * E->M.vpt[GMVINDEX(j,k)]
					* E->boundary.det[m][side][k][bel] * g_1d[k].weight[dims-1];

			}
		}
	}
}

void get_elt_tr_pseudo_surf(struct All_variables *E, int bel, int side, double elt_tr[24], int m)
{
	void construct_side_c3x3matrix_el();

	const int dims=E->mesh.nsd;
	const int ends1=enodes[dims-1];
	const int oned = onedvpoints[dims];

	struct CC Cc;
	struct CCX Ccx;

	const unsigned sbc_flag[4] = {0,SBX,SBY,SBZ};

	double traction[4][5],traction_at_gs[4][5], value, tmp;
	int j, b, p, k, a, nodea, d;
	int el = E->boundary.element[m][bel];
	int flagged;
	int found = 0;

	const float rho = E->data.density;
	const float g = E->data.grav_acc;
	const float R = 6371000.0;
	const float eta = E->data.ref_viscosity;
	const float kappa = E->data.kappa0;
	const float factor = 1.0e+00;
	int nodeas;

	if(E->control.side_sbcs)
		for(a=1;a<=ends1;a++)  {
			nodea = E->ien[m][el].node[ sidenodes[side][a] ];
			for(d=1;d<=dims;d++) {
				value = E->sbc.SB[m][side][d][ E->sbc.node[m][nodea] ];
				flagged = (E->node[m][nodea] & sbc_flag[d]) && (value);
				found |= flagged;
				traction[d][a] = ( flagged ? value : 0.0 );
			}
		}
	else {
		if( side == SIDE_TOP && E->parallel.me_loc[3]==E->parallel.nprocz-1 && (el%E->lmesh.elz==0)) {
			for(a=1;a<=ends1;a++)  {
				nodea = E->ien[m][el].node[ sidenodes[side][a] ];
				nodeas = E->ien[m][el].node[ sidenodes[side][a] ]/E->lmesh.noz;
				traction[1][a] = 0.0;
				traction[2][a] = 0.0;
				traction[3][a] = -1.0*factor*rho*g*(R*R*R)/(eta*kappa)
					*(E->slice.freesurf[m][nodeas]+E->sphere.cap[m].V[3][nodea]*E->advection.timestep);
				if(E->parallel.me==11 && nodea==3328)
					fprintf(stderr,"traction=%e vnew=%e timestep=%e coeff=%e\n",traction[3][a],E->sphere.cap[m].V[3][nodea],E->advection.timestep,-1.0*factor*rho*g*(R*R*R)/(eta*kappa));
				found = 1;
#if 0
				if(found && E->parallel.me==1)
					fprintf(stderr,"me=%d bel=%d el=%d side=%d TOP=%d a=%d sidenodes=%d ien=%d noz=%d nodea=%d traction=%e %e %e\n",
						E->parallel.me,bel,el,side,SIDE_TOP,a,sidenodes[side][a],
						E->ien[m][el].node[ sidenodes[side][a] ],E->lmesh.noz,
						nodea,traction[1][a],traction[2][a],traction[3][a]);

#endif
			}
		}
		else {
			for(a=1;a<=ends1;a++)  {
				nodea = E->ien[m][el].node[ sidenodes[side][a] ];
				for(d=1;d<=dims;d++) {
					value = E->sphere.cap[m].VB[d][nodea];
					flagged = (E->node[m][nodea] & sbc_flag[d]) && (value);
					found |= flagged;
					traction[d][a] = ( flagged ? value : 0.0 );
				}
			}
		}
	}

	/* skip the following computation if no sbc_flag is set
	   or value of sbcs are zero */
	if(!found) return;

	/* compute traction at each int point */
	construct_side_c3x3matrix_el(E,el,&Cc,&Ccx,
				     E->mesh.levmax,m,0,side);

	for(k=1;k<=oned;k++)
		for(d=1;d<=dims;d++) {
			traction_at_gs[d][k] = 0.0;
			for(j=1;j<=ends1;j++)
				traction_at_gs[d][k] += traction[d][j] * E->M.vpt[GMVINDEX(j,k)] ;
		}

	for(j=1;j<=ends1;j++) {
		a = sidenodes[side][j];
		for(d=1;d<=dims;d++) {
			p = dims*(a-1)+d-1;
			for(k=1;k<=oned;k++) {
				tmp = 0.0;
				for(b=1;b<=dims;b++)
					tmp += traction_at_gs[b][k] * Cc.vpt[BVINDEX(b,d,a,k)];

				elt_tr[p] += tmp * E->M.vpt[GMVINDEX(j,k)]
					* E->boundary.det[m][side][k][bel] * g_1d[k].weight[dims-1];

			}
		}
	}
}


/* =================================================================
 subroutine to get augmented lagrange part of stiffness matrix
================================================================== */

void get_aug_k(E,el,elt_k,level,m)
     struct All_variables *E;
     int el,m;
     double elt_k[24*24];
     int level;
{
     int i,p[9],a,b,nodea,nodeb;
     double Visc;

     const int n=loc_mat_size[E->mesh.nsd];
     const int ends=enodes[E->mesh.nsd];
     const int vpts=vpoints[E->mesh.nsd];
     const int dims=E->mesh.nsd;

     Visc = 0.0;
     for(a=1;a<=vpts;a++) {
	  p[a] = (a-1)*dims;
	  Visc += E->EVI[level][m][(el-1)*vpts+a];
       }
     Visc = Visc/vpts;

     for(a=1;a<=ends;a++) {
        nodea=E->IEN[level][m][el].node[a];
        for(b=1;b<=ends;b++) {
           nodeb=E->IEN[level][m][el].node[b];      /* for Kab dims*dims  */
	   i = (a-1)*n*dims+(b-1)*dims;
	   elt_k[i  ] += Visc*E->control.augmented*
	              E->elt_del[level][m][el].g[p[a]][0]*
		      E->elt_del[level][m][el].g[p[b]][0];   /*for 11 */
	   elt_k[i+1] += Visc*E->control.augmented*
	              E->elt_del[level][m][el].g[p[a]][0]*
		      E->elt_del[level][m][el].g[p[b]+1][0];  /* for 12 */
	   elt_k[i+n] += Visc*E->control.augmented*
	              E->elt_del[level][m][el].g[p[a]+1][0]*
		      E->elt_del[level][m][el].g[p[b]][0];    /* for 21 */
	   elt_k[i+n+1] += Visc*E->control.augmented*
	              E->elt_del[level][m][el].g[p[a]+1][0]*
		      E->elt_del[level][m][el].g[p[b]+1][0];  /* for 22 */

           if(3==dims) {
	       elt_k[i+2] += Visc*E->control.augmented*
	              E->elt_del[level][m][el].g[p[a]][0]*
		      E->elt_del[level][m][el].g[p[b]+2][0];  /* for 13 */
	       elt_k[i+n+2] += Visc*E->control.augmented*
	              E->elt_del[level][m][el].g[p[a]+1][0]*
		      E->elt_del[level][m][el].g[p[b]+2][0];  /* for 23 */
	       elt_k[i+n+n] += Visc*E->control.augmented*
	              E->elt_del[level][m][el].g[p[a]+2][0]*
		      E->elt_del[level][m][el].g[p[b]][0];    /* for 31 */
	       elt_k[i+n+n+1] += Visc*E->control.augmented*
	              E->elt_del[level][m][el].g[p[a]+2][0]*
		      E->elt_del[level][m][el].g[p[b]+1][0];  /* for 32 */
	       elt_k[i+n+n+2] += Visc*E->control.augmented*
	              E->elt_del[level][m][el].g[p[a]+2][0]*
		      E->elt_del[level][m][el].g[p[b]+2][0];  /* for 33 */
               }
           }
       }

   return;
   }


/* Add gamma * G^T * M_p^-1 * G for strict ALA, where the elementwise
   P0 pressure mass is M_p=V_e and G=D+C. */
void get_ala_aug_k(struct All_variables *E, int el,
                   double elt_k[24*24], int level, int m)
{
    int a,b,da,db,ia,ib;
    double ga,gb,scale,volume;
    const int dims=E->mesh.nsd;
    const int ends=enodes[dims];
    const int n=loc_mat_size[dims];

    volume=E->ECO[level][m][el].area;
    if(!isfinite(volume) || volume<=0.0)
        myerror(E,"Invalid P0 pressure mass in strict-ALA augmentation");
    scale=E->control.ala_augmented_lagrangian_gamma/volume;

    for(a=1;a<=ends;a++)
        for(da=0;da<dims;da++) {
            ia=(a-1)*dims+da;
            ga=E->elt_del[level][m][el].g[ia][0]
              +E->elt_c[level][m][el].c[ia][0];
            for(b=1;b<=ends;b++)
                for(db=0;db<dims;db++) {
                    ib=(b-1)*dims+db;
                    gb=E->elt_del[level][m][el].g[ib][0]
                      +E->elt_c[level][m][el].c[ib][0];
                    elt_k[ia*n+ib] += scale*ga*gb;
                }
        }
}


/* Apply only the strict-ALA augmentation.  Keeping this operation separate
   permits an exact audit of the original, unaugmented momentum equation. */
void assemble_ala_augmented_u(struct All_variables *E, double **u,
                              double **Au, int level, int strip_bcs)
{
    int m,e,a,d,node,eq,ia;
    double gu,scale,volume;
    const int dims=E->mesh.nsd;
    const int ends=enodes[dims];
    const int nel=E->lmesh.NEL[level];
    const int neq=E->lmesh.NEQ[level];
    void strip_bcs_from_residual();

    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        for(eq=0;eq<neq;eq++)
            Au[m][eq]=0.0;
        for(e=1;e<=nel;e++) {
            gu=0.0;
            for(a=1;a<=ends;a++) {
                node=E->IEN[level][m][e].node[a];
                for(d=0;d<dims;d++) {
                    ia=(a-1)*dims+d;
                    eq=E->ID[level][m][node].doff[d+1];
                    gu += (E->elt_del[level][m][e].g[ia][0]
                          +E->elt_c[level][m][e].c[ia][0])*u[m][eq];
                }
            }
            volume=E->ECO[level][m][e].area;
            if(!isfinite(volume) || volume<=0.0)
                myerror(E,"Invalid P0 pressure mass in strict-ALA action");
            scale=E->control.ala_augmented_lagrangian_gamma*gu/volume;
            for(a=1;a<=ends;a++) {
                node=E->IEN[level][m][e].node[a];
                for(d=0;d<dims;d++) {
                    ia=(a-1)*dims+d;
                    eq=E->ID[level][m][node].doff[d+1];
                    Au[m][eq] += scale*(E->elt_del[level][m][e].g[ia][0]
                                      +E->elt_c[level][m][e].c[ia][0]);
                }
            }
        }
    }
    (E->solver.exchange_id_d)(E,Au,level);
    if(strip_bcs)
        strip_bcs_from_residual(E,Au,level);
}


void assemble_unaugmented_del2_u(struct All_variables *E, double **u,
                                 double **Au, int level, int strip_bcs)
{
    int m,i;
    double *aug[NCS];
    const int neq=E->lmesh.NEQ[level];
    void assemble_del2_u();
    void assemble_ala_augmented_u();

    assemble_del2_u(E,u,Au,level,strip_bcs);
    if(E->control.ala_augmented_lagrangian_gamma==0.0)
        return;

    for(m=1;m<=E->sphere.caps_per_proc;m++)
        aug[m]=(double *)malloc(neq*sizeof(double));
    assemble_ala_augmented_u(E,u,aug,level,strip_bcs);
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        for(i=0;i<neq;i++)
            Au[m][i] -= aug[m][i];
        free(aug[m]);
    }
}


/* Apply the complete homogeneous strict-ALA velocity-pressure block.  This
 * is deliberately independent of assemble_forces(): pressure is an operator
 * unknown here, never a cached buoyancy contribution.  Keeping this action
 * allocation-free makes it suitable for a future coupled Krylov method and
 * coupled multigrid on every existing level.
 *
 * The caller owns five distinct fields.  In particular, velocity_work cannot
 * alias velocity or momentum because the nodal K_gamma action uses the neq
 * sentinel in its input and assemble_grad_rho_p() overwrites its output. */
void apply_ala_coupled_operator(struct All_variables *E,
                                double **velocity, double **pressure,
                                double **momentum, double **continuity,
                                double **velocity_work, int level)
{
    int m,i;
    const int neq=E->lmesh.NEQ[level];
    void assemble_del2_u();

    if(velocity==momentum || velocity==velocity_work ||
       momentum==velocity_work || pressure==continuity)
        myerror(E,"Aliased field in strict-ALA coupled operator");

    assemble_del2_u(E,velocity,momentum,level,1);
    assemble_grad_rho_p(E,pressure,velocity_work,level);
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(i=0;i<neq;i++)
            momentum[m][i] += velocity_work[m][i];
    assemble_div_rho_u(E,velocity,continuity,level);
}


/* Recover the pressure-independent force needed by a monolithic block solve.
 * The legacy finest-level assembly satisfies
 *
 *     E->F = f_external - C^T p.
 *
 * Since (D+C)^T p-D^T p=C^T p, the expression below reconstructs
 * f_external without mutating the current pressure or changing the legacy
 * E->F contract.  A later coupled solver can assemble this once per operator
 * lifecycle and then keep pressure entirely on the left-hand side. */
void assemble_ala_pressure_independent_force(struct All_variables *E,
                                             double **force,
                                             double **velocity_work,
                                             int level)
{
    int m,i;
    const int neq=E->lmesh.NEQ[level];
    void assemble_forces();
    void assemble_grad_p();
    void strip_bcs_from_residual();

    if(level!=E->mesh.levmax)
        myerror(E,"ALA pressure-independent force requires finest level");
    if(force==velocity_work || force==E->F || velocity_work==E->F)
        myerror(E,"Aliased field in ALA pressure-independent force");

    assemble_forces(E,0);
    assemble_grad_rho_p(E,E->P,velocity_work,level);
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(i=0;i<neq;i++)
            force[m][i]=E->F[m][i]+velocity_work[m][i];

    assemble_grad_p(E,E->P,velocity_work,level);
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(i=0;i<neq;i++)
            force[m][i]-=velocity_work[m][i];
    strip_bcs_from_residual(E,force,level);
}


static void ala_require_adjacent_factor_two(struct All_variables *E,
                                             int coarse_level);


void ala_coupled_prolong_velocity(struct All_variables *E, int coarse_level,
                                  double **coarse, double **fine)
{
    void interp_vector();
    void strip_bcs_from_residual();

    if(coarse_level<E->mesh.levmin || coarse_level>=E->mesh.levmax)
        myerror(E,"Invalid coupled velocity prolongation level");
    interp_vector(E,coarse_level,coarse,fine);
    strip_bcs_from_residual(E,fine,coarse_level+1);
}


void ala_coupled_restrict_velocity(struct All_variables *E, int fine_level,
                                   double **fine, double **coarse)
{
    int m,i,j,k,d,node,fine_eq,coarse_node,coarse_eq;
    int ix,iy,iz,nx,ny,nz;
    int cx[2],cy[2],cz[2];
    double wx[2],wy[2],wz[2],weight,value,x1,x2;
    int coarse_level=fine_level-1;
    int fnox,fnoy,fnoz,cnox,cnoz,coarse_neq;
    void from_rtf_to_xyz();
    void from_xyz_to_rtf();
    void strip_bcs_from_residual();

    if(fine_level<=E->mesh.levmin || fine_level>E->mesh.levmax)
        myerror(E,"Invalid coupled velocity restriction level");
    ala_require_adjacent_factor_two(E,coarse_level);
    fnox=E->lmesh.NOX[fine_level];
    fnoy=E->lmesh.NOY[fine_level];
    fnoz=E->lmesh.NOZ[fine_level];
    cnox=E->lmesh.NOX[coarse_level];
    cnoz=E->lmesh.NOZ[coarse_level];
    coarse_neq=E->lmesh.NEQ[coarse_level];

    /* interp_vector is T_f^T L T_c: rotate coarse r-t-phi values to
     * Cartesian, trilinearly interpolate, then rotate back on the fine
     * mesh.  Its owned-DOF transpose is T_c^T L^T T_f.  Accumulate only
     * owned fine nodes, assemble the coarse Cartesian contributions, and
     * finally project the result into the coarse constrained space. */
    from_rtf_to_xyz(E,fine_level,fine,E->temp);
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(i=0;i<=coarse_neq;i++)
            E->temp1[m][i]=0.0;

    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(k=1;k<=fnoy;k++)
            for(i=1;i<=fnox;i++)
                for(j=1;j<=fnoz;j++) {
                    node=j+(i-1)*fnoz+(k-1)*fnox*fnoz;
                    if(E->NODE[fine_level][m][node] & SKIP)
                        continue;
                    if(i%2) {
                        nx=1; cx[0]=(i+1)/2; wx[0]=1.0;
                    }
                    else {
                        nx=2; cx[0]=i/2; cx[1]=cx[0]+1;
                        wx[0]=wx[1]=0.5;
                    }
                    if(k%2) {
                        ny=1; cy[0]=(k+1)/2; wy[0]=1.0;
                    }
                    else {
                        ny=2; cy[0]=k/2; cy[1]=cy[0]+1;
                        wy[0]=wy[1]=0.5;
                    }
                    if(j%2) {
                        nz=1; cz[0]=(j+1)/2; wz[0]=1.0;
                    }
                    else {
                        nz=2; cz[0]=j/2; cz[1]=cz[0]+1;
                        x1=E->sphere.R[fine_level][j]
                           -E->sphere.R[fine_level][j-1];
                        x2=E->sphere.R[fine_level][j+1]
                           -E->sphere.R[fine_level][j];
                        wz[0]=x2/(x1+x2);
                        wz[1]=1.0-wz[0];
                    }
                    for(d=1;d<=E->mesh.nsd;d++) {
                        fine_eq=E->ID[fine_level][m][node].doff[d];
                        value=E->temp[m][fine_eq];
                        for(iy=0;iy<ny;iy++)
                            for(ix=0;ix<nx;ix++)
                                for(iz=0;iz<nz;iz++) {
                                    coarse_node=cz[iz]
                                      +(cx[ix]-1)*cnoz
                                      +(cy[iy]-1)*cnox*cnoz;
                                    coarse_eq=E->ID[coarse_level][m]
                                      [coarse_node].doff[d];
                                    weight=wx[ix]*wy[iy]*wz[iz];
                                    E->temp1[m][coarse_eq] += weight*value;
                                }
                    }
                }

    (E->solver.exchange_id_d)(E,E->temp1,coarse_level);
    from_xyz_to_rtf(E,coarse_level,E->temp1,coarse);
    strip_bcs_from_residual(E,coarse,coarse_level);
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        coarse[m][coarse_neq]=0.0;
}


static void ala_require_adjacent_factor_two(struct All_variables *E,
                                             int coarse_level)
{
    int fine_level=coarse_level+1;

    if(coarse_level<E->mesh.levmin || fine_level>E->mesh.levmax ||
       E->lmesh.ELX[fine_level]!=2*E->lmesh.ELX[coarse_level] ||
       E->lmesh.ELY[fine_level]!=2*E->lmesh.ELY[coarse_level] ||
       E->lmesh.ELZ[fine_level]!=2*E->lmesh.ELZ[coarse_level])
        myerror(E,"Coupled P0 transfer requires adjacent factor-two levels");
}


void ala_coupled_prolong_pressure_p0(struct All_variables *E,
                                     int coarse_level,
                                     double **coarse, double **fine)
{
    int m,ex,ey,ez,e,cx,cy,cz,ce;
    int fine_level=coarse_level+1;
    int elx,ely,elz,celx,celz;

    ala_require_adjacent_factor_two(E,coarse_level);
    elx=E->lmesh.ELX[fine_level];
    ely=E->lmesh.ELY[fine_level];
    elz=E->lmesh.ELZ[fine_level];
    celx=E->lmesh.ELX[coarse_level];
    celz=E->lmesh.ELZ[coarse_level];
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        fine[m][0]=0.0;
        for(ey=1;ey<=ely;ey++)
            for(ex=1;ex<=elx;ex++)
                for(ez=1;ez<=elz;ez++) {
                    e=ez+(ex-1)*elz+(ey-1)*elz*elx;
                    cx=(ex-1)/2+1;
                    cy=(ey-1)/2+1;
                    cz=(ez-1)/2+1;
                    ce=cz+(cx-1)*celz+(cy-1)*celz*celx;
                    fine[m][e]=coarse[m][ce];
                }
    }
}


void ala_coupled_restrict_pressure_p0_transpose(struct All_variables *E,
                                                int fine_level,
                                                double **fine,
                                                double **coarse)
{
    int m,ex,ey,ez,e,cx,cy,cz,ce;
    int coarse_level=fine_level-1;
    int elx,ely,elz,celx,celz,cnpno;

    ala_require_adjacent_factor_two(E,coarse_level);
    elx=E->lmesh.ELX[fine_level];
    ely=E->lmesh.ELY[fine_level];
    elz=E->lmesh.ELZ[fine_level];
    celx=E->lmesh.ELX[coarse_level];
    celz=E->lmesh.ELZ[coarse_level];
    cnpno=E->lmesh.NPNO[coarse_level];
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        for(ce=0;ce<=cnpno;ce++)
            coarse[m][ce]=0.0;
        for(ey=1;ey<=ely;ey++)
            for(ex=1;ex<=elx;ex++)
                for(ez=1;ez<=elz;ez++) {
                    e=ez+(ex-1)*elz+(ey-1)*elz*elx;
                    cx=(ex-1)/2+1;
                    cy=(ey-1)/2+1;
                    cz=(ez-1)/2+1;
                    ce=cz+(cx-1)*celz+(cy-1)*celz*celx;
                    coarse[m][ce] += fine[m][e];
                }
    }
}


/* Matrix-free Galerkin action on one adjacent coupled level:
 *
 *     A_c^RAP x_c = P^T A_f P x_c.
 *
 * This audit path deliberately reuses the production velocity/P0 transfers
 * and the complete strict-ALA block action.  It observes the operator that a
 * true coupled Galerkin hierarchy would use without changing the current
 * rediscretized V-cycle. */
static void ala_apply_coupled_galerkin_rap(
    struct All_variables *E, const struct ala_block_vector *coarse,
    struct ala_block_vector *coarse_action,
    struct ala_block_vector *fine_prolong,
    struct ala_block_vector *fine_action,
    struct ala_block_vector *fine_work, int fine_level)
{
    int coarse_level=fine_level-1;

    if(coarse->level!=coarse_level ||
       coarse_action->level!=coarse_level ||
       fine_prolong->level!=fine_level ||
       fine_action->level!=fine_level || fine_work->level!=fine_level)
        myerror(E,"Mismatched level in coupled Galerkin RAP action");
    ala_block_vector_zero(E,fine_prolong);
    ala_block_vector_zero(E,fine_action);
    ala_block_vector_zero(E,fine_work);
    ala_block_vector_zero(E,coarse_action);
    ala_coupled_prolong_velocity(E,coarse_level,coarse->velocity,
                                 fine_prolong->velocity);
    ala_coupled_prolong_pressure_p0(E,coarse_level,coarse->pressure,
                                    fine_prolong->pressure);
    apply_ala_coupled_operator(
        E,fine_prolong->velocity,fine_prolong->pressure,
        fine_action->velocity,fine_action->pressure,
        fine_work->velocity,fine_level);
    ala_coupled_restrict_velocity(E,fine_level,fine_action->velocity,
                                  coarse_action->velocity);
    ala_coupled_restrict_pressure_p0_transpose(
        E,fine_level,fine_action->pressure,coarse_action->pressure);
}


static void ala_fill_level_probe(struct All_variables *E,
                                 struct ala_block_vector *probe,
                                 double phase)
{
    int level=probe->level;
    int m,node,e,d,eq,nz;
    double radius,value;
    const unsigned int vbc_flag[4]={0,VBX,VBY,VBZ};

    ala_block_vector_zero(E,probe);
    for(m=1;m<=E->sphere.caps_per_proc;m++) {
        for(node=1;node<=E->lmesh.NNO[level];node++) {
            /* This is a primal Krylov probe, not an assembled residual.
             * Prescribe the same analytic radial value independently on
             * every replica of a shared node.  Using additive exchange to
             * construct the probe can multiply corner values along several
             * communication routes and manufacture an adjoint defect before
             * either G or G^T is applied. */
            radius=E->SX[level][m][3][node];
            for(d=1;d<=E->mesh.nsd;d++) {
                eq=E->ID[level][m][node].doff[d];
                value=(d==3) ? sin(phase+2.0*radius) : 0.0;
                if(E->NODE[level][m][node] & vbc_flag[d])
                    value=0.0;
                probe->velocity[m][eq]=value;
            }
        }
        for(e=1;e<=E->lmesh.NPNO[level];e++) {
            nz=(e-1)%E->lmesh.ELZ[level]+1;
            probe->pressure[m][e]=cos(phase+0.17*nz);
        }
    }
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        probe->velocity[m][E->lmesh.NEQ[level]]=0.0;
}


static double ala_norm_scaled_adjoint_defect(double left, double right,
                                             double left_vector_norm,
                                             double left_action_norm,
                                             double right_vector_norm,
                                             double right_action_norm)
{
    double scale=left_vector_norm*left_action_norm
                +right_vector_norm*right_action_norm;

    return fabs(left-right)/max(scale,1.0e-300);
}


static double ala_restricted_beta(struct All_variables *E, int level, int nz,
                                  double *radial_width)
{
    int fine_first,fine_last,fine_nz;
    double beta=0.0;
    double dr,dr_total=0.0;

    fine_first=((nz-1)*E->lmesh.elz)/E->lmesh.ELZ[level]+1;
    fine_last=(nz*E->lmesh.elz)/E->lmesh.ELZ[level];
    for(fine_nz=fine_first;fine_nz<=fine_last;fine_nz++) {
        dr=E->sx[1][3][fine_nz+1]-E->sx[1][3][fine_nz];
        beta += E->refstate.ala_beta[fine_nz]*dr;
        dr_total += dr;
    }
    if(!isfinite(dr_total) || dr_total<=0.0)
        myerror(E,"Invalid radial width in coupled beta audit");
    if(radial_width!=NULL)
        *radial_width=dr_total;
    return beta/dr_total;
}


void audit_ala_coupled_multilevel_contracts(struct All_variables *E)
{
    static int completed=0;
    int level,m,e,nz,child;
    int local_invalid,global_invalid,local_skips,global_skips,output_index;
    double left,right,unused,k_defect,g_defect,v_transfer_defect,p_transfer_defect;
    double k_left,k_right,g_left,g_right,g_element_right;
    double g_exchange_right,g_stripped_right;
    double local_g_element_right,g_multiplicity_right;
    double local_g_multiplicity_right;
    double g_element_defect,g_exchange_defect,g_stripped_defect;
    double g_multiplicity_defect;
    double x_velocity_norm,x_pressure_norm,y_velocity_norm,y_pressure_norm;
    double Ax_velocity_norm,Ax_pressure_norm,Ay_velocity_norm,Ay_pressure_norm;
    double coarse_velocity_norm,coarse_pressure_norm;
    double coarse_action_velocity_norm,coarse_action_pressure_norm;
    double rap_left,rap_right,rap_symmetry_defect;
    double rap_x_norm,rap_y_norm,redis_norm,rap_norm,difference_norm;
    double redis_velocity_norm,redis_pressure_norm;
    double rap_velocity_norm,rap_pressure_norm;
    double difference_velocity_norm,difference_pressure_norm;
    double velocity_action_difference,pressure_action_difference;
    double block_action_difference;
    double beta,beta_min,beta_max,volume_min,volume_max;
    double child_beta,child_width,nested_beta,nested_width;
    double local_beta_nested_defect,global_beta_nested_defect;
    double local_bounds[4],global_min[2],global_max[2];
    struct ala_block_vector *x,*y,*Ax,*Ay,*coarse,*coarse_action;
    struct ala_block_vector *multiplicity;
    struct ala_block_vector *coarse_x,*coarse_y,*rap_x,*rap_y;
    struct ala_block_vector *redis_x,*redis_work,*difference;
    struct ala_block_vector *fine_prolong,*fine_action,*fine_work;
    FILE *output;
    void assemble_del2_u();
    void strip_bcs_from_residual();

    if(completed)
        return;
    completed=1;
    for(level=E->mesh.levmin;level<=E->mesh.levmax;level++) {
        x=ala_block_vector_create(E,level);
        y=ala_block_vector_create(E,level);
        Ax=ala_block_vector_create(E,level);
        Ay=ala_block_vector_create(E,level);
        multiplicity=ala_block_vector_create(E,level);
        ala_fill_level_probe(E,x,0.31);
        ala_fill_level_probe(E,y,0.79);
        ala_block_vector_component_norms(
            E,x,&x_velocity_norm,&x_pressure_norm);
        ala_block_vector_component_norms(
            E,y,&y_velocity_norm,&y_pressure_norm);

        assemble_del2_u(E,x->velocity,Ax->velocity,level,1);
        assemble_del2_u(E,y->velocity,Ay->velocity,level,1);
        ala_block_vector_component_products(E,x,Ay,&k_left,&unused);
        ala_block_vector_component_products(E,y,Ax,&k_right,&unused);
        ala_block_vector_component_norms(
            E,Ax,&Ax_velocity_norm,&Ax_pressure_norm);
        ala_block_vector_component_norms(
            E,Ay,&Ay_velocity_norm,&Ay_pressure_norm);
        k_defect=ala_norm_scaled_adjoint_defect(
            k_left,k_right,x_velocity_norm,Ay_velocity_norm,
            y_velocity_norm,Ax_velocity_norm);

        assemble_grad_rho_p(E,x->pressure,Ax->velocity,level);
        assemble_div_rho_u(E,y->velocity,Ay->pressure,level);
        ala_block_vector_component_products(E,y,Ax,&g_right,&unused);
        ala_block_vector_component_products(E,x,Ay,&unused,&g_left);
        ala_block_vector_component_norms(
            E,Ax,&Ax_velocity_norm,&Ax_pressure_norm);
        ala_block_vector_component_norms(
            E,Ay,&Ay_velocity_norm,&Ay_pressure_norm);
        g_defect=ala_norm_scaled_adjoint_defect(
            g_left,g_right,x_pressure_norm,Ay_pressure_norm,
            y_velocity_norm,Ax_velocity_norm);

        /* The legacy SKIP owner mask and the additive exchange routes were
         * designed independently.  Build the actual sharing multiplicity by
         * exchanging one contribution from every local nodal copy, then use
         * its reciprocal as a partition of unity over assembled copies. */
        for(m=1;m<=E->sphere.caps_per_proc;m++) {
            for(e=0;e<E->lmesh.NEQ[level];e++)
                multiplicity->velocity[m][e]=1.0;
            multiplicity->velocity[m][E->lmesh.NEQ[level]]=0.0;
        }
        (E->solver.exchange_id_d)(E,multiplicity->velocity,level);
        local_g_multiplicity_right=0.0;
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(e=0;e<E->lmesh.NEQ[level];e++) {
                if(!isfinite(multiplicity->velocity[m][e]) ||
                   multiplicity->velocity[m][e]<=0.0)
                    myerror(E,"Invalid velocity sharing multiplicity");
                local_g_multiplicity_right += y->velocity[m][e]
                    *Ax->velocity[m][e]/multiplicity->velocity[m][e];
            }
        MPI_Allreduce(&local_g_multiplicity_right,&g_multiplicity_right,1,
                      MPI_DOUBLE,MPI_SUM,E->parallel.world);
        g_multiplicity_defect=fabs(g_left-g_multiplicity_right)
            /max(fabs(g_left)+fabs(g_multiplicity_right),1.0e-300);

        /* Separate the element transpose from distributed additive assembly.
         * The unexchanged local G^T contributions dotted with the replicated
         * velocity must equal the elementwise p^T G u exactly. */
        assemble_grad_rho_p_local_terms(
            E,x->pressure,Ax->velocity,level);
        local_g_element_right=0.0;
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(e=0;e<E->lmesh.NEQ[level];e++)
                local_g_element_right += y->velocity[m][e]
                                       *Ax->velocity[m][e];
        MPI_Allreduce(&local_g_element_right,&g_element_right,1,
                      MPI_DOUBLE,MPI_SUM,E->parallel.world);
        g_element_defect=fabs(g_left-g_element_right)
            /max(fabs(g_left)+fabs(g_element_right),1.0e-300);

        /* The production transpose wrapper performs additive assembly and
         * essential-BC stripping consecutively.  Audit the intermediate
         * assembled vector before stripping so a communication defect cannot
         * be confused with a constrained-space adjoint contract. */
        (E->solver.exchange_id_d)(E,Ax->velocity,level);
        ala_block_vector_component_products(
            E,y,Ax,&g_exchange_right,&unused);
        ala_block_vector_component_norms(
            E,Ax,&Ax_velocity_norm,&Ax_pressure_norm);
        g_exchange_defect=ala_norm_scaled_adjoint_defect(
            g_left,g_exchange_right,x_pressure_norm,Ay_pressure_norm,
            y_velocity_norm,Ax_velocity_norm);
        strip_bcs_from_residual(E,Ax->velocity,level);
        ala_block_vector_component_products(
            E,y,Ax,&g_stripped_right,&unused);
        ala_block_vector_component_norms(
            E,Ax,&Ax_velocity_norm,&Ax_pressure_norm);
        g_stripped_defect=ala_norm_scaled_adjoint_defect(
            g_left,g_stripped_right,x_pressure_norm,Ay_pressure_norm,
            y_velocity_norm,Ax_velocity_norm);

        v_transfer_defect=0.0;
        p_transfer_defect=0.0;
        if(level>E->mesh.levmin) {
            coarse=ala_block_vector_create(E,level-1);
            coarse_action=ala_block_vector_create(E,level-1);
            ala_fill_level_probe(E,coarse,0.53);
            ala_coupled_prolong_velocity(E,level-1,coarse->velocity,
                                         Ax->velocity);
            ala_coupled_restrict_velocity(E,level,y->velocity,
                                          coarse_action->velocity);
            ala_block_vector_component_products(E,y,Ax,&left,&unused);
            ala_block_vector_component_products(E,coarse,coarse_action,
                                                &right,&unused);
            ala_block_vector_component_norms(
                E,Ax,&Ax_velocity_norm,&Ax_pressure_norm);
            ala_block_vector_component_norms(
                E,coarse,&coarse_velocity_norm,&coarse_pressure_norm);
            ala_block_vector_component_norms(
                E,coarse_action,&coarse_action_velocity_norm,
                &coarse_action_pressure_norm);
            v_transfer_defect=ala_norm_scaled_adjoint_defect(
                left,right,y_velocity_norm,Ax_velocity_norm,
                coarse_velocity_norm,coarse_action_velocity_norm);
            ala_coupled_prolong_pressure_p0(E,level-1,coarse->pressure,
                                            Ax->pressure);
            ala_coupled_restrict_pressure_p0_transpose(
                E,level,y->pressure,coarse_action->pressure);
            ala_block_vector_component_products(E,y,Ax,&unused,&left);
            ala_block_vector_component_products(E,coarse,coarse_action,
                                                &unused,&right);
            ala_block_vector_component_norms(
                E,Ax,&Ax_velocity_norm,&Ax_pressure_norm);
            ala_block_vector_component_norms(
                E,coarse_action,&coarse_action_velocity_norm,
                &coarse_action_pressure_norm);
            p_transfer_defect=ala_norm_scaled_adjoint_defect(
                left,right,y_pressure_norm,Ax_pressure_norm,
                coarse_pressure_norm,coarse_action_pressure_norm);
            ala_block_vector_destroy(E,coarse);
            ala_block_vector_destroy(E,coarse_action);

            /* Compare the true matrix-free Galerkin action P^T A_f P with
             * the operator currently rediscretized directly on this coarse
             * mesh.  The comparison is observe-only: the production V-cycle
             * remains unchanged until the RAP contract is measured. */
            coarse_x=ala_block_vector_create(E,level-1);
            coarse_y=ala_block_vector_create(E,level-1);
            rap_x=ala_block_vector_create(E,level-1);
            rap_y=ala_block_vector_create(E,level-1);
            redis_x=ala_block_vector_create(E,level-1);
            redis_work=ala_block_vector_create(E,level-1);
            difference=ala_block_vector_create(E,level-1);
            fine_prolong=ala_block_vector_create(E,level);
            fine_action=ala_block_vector_create(E,level);
            fine_work=ala_block_vector_create(E,level);
            ala_fill_level_probe(E,coarse_x,0.37);
            ala_fill_level_probe(E,coarse_y,0.83);
            ala_apply_coupled_galerkin_rap(
                E,coarse_x,rap_x,fine_prolong,fine_action,fine_work,level);
            ala_apply_coupled_galerkin_rap(
                E,coarse_y,rap_y,fine_prolong,fine_action,fine_work,level);
            rap_left=ala_block_vector_dot(E,coarse_x,rap_y,1.0,1.0);
            rap_right=ala_block_vector_dot(E,coarse_y,rap_x,1.0,1.0);
            rap_x_norm=ala_block_vector_norm(E,rap_x,1.0,1.0);
            rap_y_norm=ala_block_vector_norm(E,rap_y,1.0,1.0);
            rap_symmetry_defect=fabs(rap_left-rap_right)/max(
                ala_block_vector_norm(E,coarse_x,1.0,1.0)*rap_y_norm
               +ala_block_vector_norm(E,coarse_y,1.0,1.0)*rap_x_norm,
                1.0e-300);
            apply_ala_coupled_operator(
                E,coarse_x->velocity,coarse_x->pressure,
                redis_x->velocity,redis_x->pressure,
                redis_work->velocity,level-1);
            ala_block_vector_copy(E,redis_x,difference);
            ala_block_vector_axpy(E,-1.0,rap_x,difference);
            redis_norm=ala_block_vector_norm(E,redis_x,1.0,1.0);
            rap_norm=ala_block_vector_norm(E,rap_x,1.0,1.0);
            difference_norm=ala_block_vector_norm(E,difference,1.0,1.0);
            block_action_difference=difference_norm
                /max(redis_norm+rap_norm,1.0e-300);
            ala_block_vector_component_norms(
                E,redis_x,&redis_velocity_norm,&redis_pressure_norm);
            ala_block_vector_component_norms(
                E,rap_x,&rap_velocity_norm,&rap_pressure_norm);
            ala_block_vector_component_norms(
                E,difference,&difference_velocity_norm,
                &difference_pressure_norm);
            velocity_action_difference=difference_velocity_norm/max(
                redis_velocity_norm+rap_velocity_norm,1.0e-300);
            pressure_action_difference=difference_pressure_norm/max(
                redis_pressure_norm+rap_pressure_norm,1.0e-300);
            if(E->parallel.me==0)
              for(output_index=0;output_index<2;output_index++) {
                output=(output_index==0) ? E->fp : stderr;
                if(output==NULL || (output_index==1 && output==E->fp))
                    continue;
                fprintf(output,
                    "ALA COUPLED GALERKIN RAP AUDIT fine_level=%d "
                    "coarse_level=%d rap_symmetry_defect=%e "
                    "rediscretized_action_difference=(velocity:%e,"
                    "pressure:%e,block:%e) rap_bilinear=(%e,%e) "
                    "action=observe_only\n",
                    level,level-1,rap_symmetry_defect,
                    velocity_action_difference,pressure_action_difference,
                    block_action_difference,rap_left,rap_right);
                fflush(output);
              }
            ala_block_vector_destroy(E,coarse_x);
            ala_block_vector_destroy(E,coarse_y);
            ala_block_vector_destroy(E,rap_x);
            ala_block_vector_destroy(E,rap_y);
            ala_block_vector_destroy(E,redis_x);
            ala_block_vector_destroy(E,redis_work);
            ala_block_vector_destroy(E,difference);
            ala_block_vector_destroy(E,fine_prolong);
            ala_block_vector_destroy(E,fine_action);
            ala_block_vector_destroy(E,fine_work);
        }

        beta_min=1.0e300;
        beta_max=-1.0e300;
        volume_min=1.0e300;
        volume_max=0.0;
        local_beta_nested_defect=0.0;
        local_invalid=0;
        local_skips=0;
        for(m=1;m<=E->sphere.caps_per_proc;m++) {
            local_skips += E->parallel.Skip_neq[level][m];
            for(e=1;e<=E->lmesh.NEL[level];e++) {
                nz=(e-1)%E->lmesh.ELZ[level]+1;
                beta=ala_restricted_beta(E,level,nz,NULL);
                if(level<E->mesh.levmax && m==1 &&
                   e<=E->lmesh.ELZ[level]) {
                    nested_beta=0.0;
                    nested_width=0.0;
                    for(child=2*nz-1;child<=2*nz;child++) {
                        child_beta=ala_restricted_beta(
                            E,level+1,child,&child_width);
                        nested_beta += child_beta*child_width;
                        nested_width += child_width;
                    }
                    nested_beta /= nested_width;
                    local_beta_nested_defect=max(
                        local_beta_nested_defect,
                        fabs(beta-nested_beta)/max(fabs(beta),1.0e-300));
                }
                if(!isfinite(beta) || beta<=0.0 ||
                   !isfinite(E->ECO[level][m][e].area) ||
                   E->ECO[level][m][e].area<=0.0)
                    local_invalid++;
                beta_min=min(beta_min,beta);
                beta_max=max(beta_max,beta);
                volume_min=min(volume_min,E->ECO[level][m][e].area);
                volume_max=max(volume_max,E->ECO[level][m][e].area);
            }
        }
        local_bounds[0]=beta_min;
        local_bounds[1]=volume_min;
        local_bounds[2]=beta_max;
        local_bounds[3]=volume_max;
        MPI_Allreduce(local_bounds,global_min,2,MPI_DOUBLE,MPI_MIN,
                      E->parallel.world);
        MPI_Allreduce(local_bounds+2,global_max,2,MPI_DOUBLE,MPI_MAX,
                      E->parallel.world);
        MPI_Allreduce(&local_invalid,&global_invalid,1,MPI_INT,MPI_SUM,
                      E->parallel.world);
        MPI_Allreduce(&local_skips,&global_skips,1,MPI_INT,MPI_SUM,
                      E->parallel.world);
        MPI_Allreduce(&local_beta_nested_defect,&global_beta_nested_defect,
                      1,MPI_DOUBLE,MPI_MAX,E->parallel.world);
        if(E->parallel.me==0)
          for(output_index=0;output_index<2;output_index++) {
            output=(output_index==0) ? E->fp : stderr;
            if(output==NULL || (output_index==1 && output==E->fp))
                continue;
            fprintf(output,
                "ALA COUPLED LEVEL CONTRACT level=%d neq=%d np0=%d "
                "velocity_probe=direct_conforming_radial_bc_projected "
                "K_symmetry_defect=%e G_adjoint_defect=%e "
                "G_element_adjoint_defect=%e "
                "G_exchange_adjoint_defect=%e "
                "G_stripped_adjoint_defect=%e "
                "G_multiplicity_adjoint_defect=%e "
                "K_bilinear=(%e,%e) G_bilinear=(%e,%e) "
                "G_element_bilinear=(%e,%e) "
                "G_exchange_bilinear=(%e,%e) "
                "G_stripped_bilinear=(%e,%e) "
                "G_multiplicity_bilinear=(%e,%e) "
                "velocity_transfer_adjoint_defect=%e "
                "pressure_transfer_adjoint_defect=%e "
                "beta_range=[%e,%e] beta_nested_defect=%e "
                "pressure_mass_range=[%e,%e] "
                "duplicate_velocity_dofs=%d invalid=%d\n",
                level,E->lmesh.NEQ[level],E->lmesh.NPNO[level],
                k_defect,g_defect,g_element_defect,g_exchange_defect,
                g_stripped_defect,g_multiplicity_defect,
                k_left,k_right,g_left,g_right,g_left,g_element_right,
                g_left,g_exchange_right,g_left,g_stripped_right,
                g_left,g_multiplicity_right,
                v_transfer_defect,p_transfer_defect,
                global_min[0],global_max[0],global_beta_nested_defect,
                global_min[1],global_max[1],
                global_skips,global_invalid);
            fflush(output);
          }
        ala_block_vector_destroy(E,x);
        ala_block_vector_destroy(E,y);
        ala_block_vector_destroy(E,Ax);
        ala_block_vector_destroy(E,Ay);
        ala_block_vector_destroy(E,multiplicity);
    }
}


/* version */
/* $Id: Element_calculations.c 8250 2007-11-08 23:28:46Z tan2 $ */

/* End of file  */
