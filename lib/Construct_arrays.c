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
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"


/*========================================================
  Function to make the IEN array for a mesh of given
  dimension. IEN is an externally defined structure array

  NOTE: this is not really general enough for new elements:
  it should be done through a pre-calculated lookup table.
  ======================================================== */

void construct_ien(E)
     struct All_variables *E;

{
  int lev,p,q,r,rr,j;
  int element,start,nel,nno;
  int elz,elx,ely,nox,noy,noz;

  const int dims=E->mesh.nsd;
  const int ends=enodes[dims];

  for (lev=E->mesh.levmax;lev>=E->mesh.levmin;lev--)  {
    for (j=1;j<=E->sphere.caps_per_proc;j++)  {

      elx = E->lmesh.ELX[lev];
      elz = E->lmesh.ELZ[lev];
      ely = E->lmesh.ELY[lev];
      nox = E->lmesh.NOX[lev];
      noz = E->lmesh.NOZ[lev];
      noy = E->lmesh.NOY[lev];
      nel=E->lmesh.NEL[lev];
      nno=E->lmesh.NNO[lev];

      for(r=1;r<=ely;r++)
        for(q=1;q<=elx;q++)
          for(p=1;p<=elz;p++)     {
             element = (r-1)*elx*elz + (q-1)*elz  + p;
             start = (r-1)*noz*nox + (q-1)*noz + p;
             for(rr=1;rr<=ends;rr++)
               E->IEN[lev][j][element].node[rr]= start
                  + offset[rr].vector[0]
                  + offset[rr].vector[1]*noz
                  + offset[rr].vector[2]*noz*nox;
	     }

      }     /* end for cap j */
    }     /* end loop for lev */


/* if(E->control.verbose)  { */
/*   for (lev=E->mesh.levmax;lev>=E->mesh.levmin;lev--)  { */
/*     fprintf(E->fp_out,"output_IEN_arrays me=%d lev=%d \n",E->parallel.me,lev); */
/*   for (j=1;j<=E->sphere.caps_per_proc;j++) { */
/*     fprintf(E->fp_out,"output_IEN_arrays me=%d %d %d\n",E->parallel.me,j,E->sphere.capid[j]); */
/*     for (i=1;i<=E->lmesh.NEL[lev];i++) */
/*        fprintf(E->fp_out,"%d %d %d %d %d %d %d %d %d\n",i,E->IEN[lev][j][i].node[1],E->IEN[lev][j][i].node[2],E->IEN[lev][j][i].node[3],E->IEN[lev][j][i].node[4],E->IEN[lev][j][i].node[5],E->IEN[lev][j][i].node[6],E->IEN[lev][j][i].node[7],E->IEN[lev][j][i].node[8]); */
/*     } */
/*     } */
/*   fflush (E->fp_out); */
/*   } */

  return;
}


/*  determine surface things */

void construct_surface( struct All_variables *E)
{
  int i, j, e, element;

  for (j=1;j<=E->sphere.caps_per_proc;j++)  {
    e = 0;
    for(element=1;element<=E->lmesh.nel;element++)
      if ( element%E->lmesh.elz==0) { /* top */
        e ++;
        E->sien[j][e].node[1] = E->ien[j][element].node[5]/E->lmesh.noz;
        E->sien[j][e].node[2] = E->ien[j][element].node[6]/E->lmesh.noz;
        E->sien[j][e].node[3] = E->ien[j][element].node[7]/E->lmesh.noz;
        E->sien[j][e].node[4] = E->ien[j][element].node[8]/E->lmesh.noz;
        E->surf_element[j][e] = element;
        }

    E->lmesh.snel = e;
    for (i=1;i<=E->lmesh.nsf;i++)
      E->surf_node[j][i] = i*E->lmesh.noz;

  }     /* end for cap j */
}


/*============================================
  Function to make the ID array for above case
  ============================================ */

void construct_id(E)
     struct All_variables *E;
{
    int i,j,k;
    int eqn_count,node,nno;
    unsigned int type,doff;
    int lev;
    void get_bcs_id_for_residual();

    const int dims=E->mesh.nsd,dofs=E->mesh.dof;
    const int ends=enodes[dims];

  for(lev=E->mesh.gridmax;lev>=E->mesh.gridmin;lev--)  {
    for(j=1;j<=E->sphere.caps_per_proc;j++)  {
      eqn_count = 0;

      for(node=1;node<=E->lmesh.NNO[lev];node++)
        for(doff=1;doff<=dims;doff++)  {
          E->ID[lev][j][node].doff[doff] = eqn_count;
          eqn_count ++;
          }

      E->lmesh.NEQ[lev] = eqn_count;

      i = 0;
      for(node=1;node<=E->lmesh.NNO[lev];node++) {
        if (E->NODE[lev][j][node] & SKIP)
        for(doff=1;doff<=dims;doff++)  {
	  i++;
          E->parallel.Skip_id[lev][j][i] = E->ID[lev][j][node].doff[doff];
          }
        }

      E->parallel.Skip_neq[lev][j] = i;

      get_bcs_id_for_residual(E,lev,j);

      }       /* end for j */
    }      /* end for lev */

    E->lmesh.neq = E->lmesh.NEQ[E->mesh.levmax];

/*     if (E->control.verbose) { */
/*       fprintf(E->fp_out,"output_ID_arrays \n"); */
/*       for(j=1;j<=E->sphere.caps_per_proc;j++)    */
/*         for (i=1;i<=E->lmesh.nno;i++) */
/*           fprintf(E->fp_out,"%d %d %d %d %d\n",eqn_count,i,E->ID[lev][j][i].doff[1],E->ID[lev][j][i].doff[2],E->ID[lev][j][i].doff[3]); */
/*       fflush(E->fp_out); */
/*       } */


    return;
    }



void get_bcs_id_for_residual(E,level,m)
    struct All_variables *E;
    int level,m;
  {

    int i,j;

    const int nno=E->lmesh.NNO[level];

   j = 0;
   for(i=1;i<=nno;i++) {
      if ( (E->NODE[level][m][i] & VBX) != 0 )  {
	j++;
        E->zero_resid[level][m][j] = E->ID[level][m][i].doff[1];
	}
      if ( (E->NODE[level][m][i] & VBY) != 0 )  {
	j++;
        E->zero_resid[level][m][j] = E->ID[level][m][i].doff[2];
	}
      if ( (E->NODE[level][m][i] & VBZ) != 0 )  {
	j++;
        E->zero_resid[level][m][j] = E->ID[level][m][i].doff[3];
	}
      }

    E->num_zero_resid[level][m] = j;

    return;
}

/*==========================================================
  Function to construct  the LM array from the ID and IEN arrays
  ========================================================== */

void construct_lm(E)
     struct All_variables *E;
{
  int i,j,a,e;
  int lev,eqn_no;
  int nel, nel2;

  const int dims=E->mesh.nsd,dofs=E->mesh.dof;
  const int ends=enodes[dims];

  return;
}


/* =====================================================
   Function to build the local node matrix indexing maps
   ===================================================== */

void construct_node_maps(E)
    struct All_variables *E;
{
    double time1,CPU_time0();

    int ii,noz,noxz,m,n,nn,lev,i,j,k,jj,kk,ia,ja,is,ie,js,je,ks,ke,doff;
    int neq,nno,dims2,matrix,nox,noy;

    const int dims=E->mesh.nsd,dofs=E->mesh.dof;
    const int ends=enodes[dims];
    int max_eqn;

  dims2 = dims-1;
  for(lev=E->mesh.gridmax;lev>=E->mesh.gridmin;lev--)
    for (m=1;m<=E->sphere.caps_per_proc;m++)             {
       neq=E->lmesh.NEQ[lev];
       nno=E->lmesh.NNO[lev];
       noxz = E->lmesh.NOX[lev]*E->lmesh.NOZ[lev];
       noz = E->lmesh.NOZ[lev];
       noy = E->lmesh.NOY[lev];
       nox = E->lmesh.NOX[lev];
       max_eqn = 14*dims;
       matrix = max_eqn*nno;

       E->Node_map[lev][m]=(int *) malloc (matrix*sizeof(int));

       for(i=0;i<matrix;i++)
	   E->Node_map[lev][m][i] = neq;  /* neq indicates an invalid eqn # */

       for (ii=1;ii<=noy;ii++)
       for (jj=1;jj<=nox;jj++)
       for (kk=1;kk<=noz;kk++)  {
	 nn = kk + (jj-1)*noz+ (ii-1)*noxz;
	 for(doff=1;doff<=dims;doff++)
	   E->Node_map[lev][m][(nn-1)*max_eqn+doff-1] = E->ID[lev][m][nn].doff[doff];

         ia = 0;
	 is=1; ie=dims2;
	 js=1; je=dims;
	 ks=1; ke=dims;
	 if (kk==1  ) ks=2;
	 if (kk==noz) ke=2;
	 if (jj==1  ) js=2;
	 if (jj==nox) je=2;
	 if (ii==1  ) is=2;
	 if (ii==noy) ie=2;
         for (i=is;i<=ie;i++)
           for (j=js;j<=je;j++)
             for (k=ks;k<=ke;k++)  {
               ja = nn-((2-i)*noxz + (2-j)*noz + 2-k);
               if (ja<nn)   {
		 ia++;
                 for (doff=1;doff<=dims;doff++)
                   E->Node_map[lev][m][(nn-1)*max_eqn+ia*dims+doff-1]=E->ID[lev][m][ja].doff[doff];
                 }
               }
         }

       E->Eqn_k1[lev][m] = (higher_precision *)malloc(matrix*sizeof(higher_precision));
       E->Eqn_k2[lev][m] = (higher_precision *)malloc(matrix*sizeof(higher_precision));
       E->Eqn_k3[lev][m] = (higher_precision *)malloc(matrix*sizeof(higher_precision));

       E->mesh.matrix_size[lev] = matrix;

       if(E->control.verbose) {
           fprintf(E->fp_out, "output Node_map lev=%d m=%d\n", lev, m);
           fprintf(E->fp_out, "neq=%d nno=%d max_eqn=%d matrix=%d\n", neq, nno, max_eqn, matrix);
           for(i=0;i<matrix;i++)
               fprintf(E->fp_out, "%d %d\n", i, E->Node_map[lev][m][i]);
       }

    }         /* end for level and m */

    return;
}


void construct_node_ks(E)
     struct All_variables *E;
{
    int m,level,i,j,k,e,a,da,eq;
    int node,node1,eqn1,eqn2,eqn3,loc0,loc1,loc2,loc3,found,element,index,pp,qq;
    int neq,nno,nel,max_eqn;

    double elt_K[24*24];
    double w1,w2,w3,ww1,ww2,ww3,zero;

    higher_precision *B1,*B2,*B3;

    void get_elt_k();
    void get_aug_k();
    void get_ala_aug_k();
    void build_diagonal_of_K();
    void parallel_process_termination();
    void myerror();

    const int dims=E->mesh.nsd,dofs=E->mesh.dof;
    const int ends=enodes[dims];
    const int lms=loc_mat_size[E->mesh.nsd];

    zero = 0.0;
    max_eqn = 14*dims;

   for(level=E->mesh.gridmax;level>=E->mesh.gridmin;level--)   {

      for(m=1;m<=E->sphere.caps_per_proc;m++)     {

        neq=E->lmesh.NEQ[level];
        nel=E->lmesh.NEL[level];
        nno=E->lmesh.NNO[level];
	for(i=0;i<neq;i++)
	    E->BI[level][m][i] = zero;
        if(E->control.ala_element_vanka_smoother)
            for(i=0;i<neq;i++)
                E->ALA_vanka_overlap_BI[level][m][i] = zero;
        for(i=0;i<E->mesh.matrix_size[level];i++) {
            E->Eqn_k1[level][m][i] = zero;
            E->Eqn_k2[level][m][i] = zero;
            E->Eqn_k3[level][m][i] = zero;
            }

        for(element=1;element<=nel;element++) {

	    get_elt_k(E,element,elt_K,level,m,0);

            if(E->control.ala_element_vanka_smoother)
                for(a=1;a<=ends;a++) {
                    node=E->IEN[level][m][element].node[a];
                    for(da=0;da<dims;da++) {
                        if((da==0 && (E->NODE[level][m][node] & VBX)) ||
                           (da==1 && (E->NODE[level][m][node] & VBY)) ||
                           (da==2 && (E->NODE[level][m][node] & VBZ)))
                            continue;
                        eq=E->ID[level][m][node].doff[da+1];
                        E->ALA_vanka_overlap_BI[level][m][eq] += 1.0;
                    }
                }

	    if (E->control.augmented_Lagr)
	         get_aug_k(E,element,elt_K,level,m);

            if(E->control.ala_augmented_lagrangian_gamma > 0.0)
                get_ala_aug_k(E,element,elt_K,level,m);

            build_diagonal_of_K(E,element,elt_K,level,m);

	    for(i=1;i<=ends;i++) {  /* i, is the node we are storing to */
	       node=E->IEN[level][m][element].node[i];

	       pp=(i-1)*dims;
	       w1=w2=w3=1.0;

	       loc0=(node-1)*max_eqn;

	       if(E->NODE[level][m][node] & VBX) w1=0.0;
	       if(E->NODE[level][m][node] & VBZ) w3=0.0;
	       if(E->NODE[level][m][node] & VBY) w2=0.0;

	       for(j=1;j<=ends;j++) { /* j is the node we are receiving from */
	         node1=E->IEN[level][m][element].node[j];

                        /* only for half of the matrix ,because of the symmetry */
                 if (node1<=node)  {

		    ww1=ww2=ww3=1.0;
		    qq=(j-1)*dims;
		    eqn1=E->ID[level][m][node1].doff[1];
		    eqn2=E->ID[level][m][node1].doff[2];
		    eqn3=E->ID[level][m][node1].doff[3];

		    if(E->NODE[level][m][node1] & VBX) ww1=0.0;
		    if(E->NODE[level][m][node1] & VBZ) ww3=0.0;
		    if(E->NODE[level][m][node1] & VBY) ww2=0.0;

		    /* search for direction 1*/

		    found=0;
		    for(k=0;k<max_eqn;k++)
		      if(E->Node_map[level][m][loc0+k] == eqn1) { /* found, index next equation */
			    index=k;
			    found++;
			    break;
			}

		    assert(found /* direction 1 */);

		    E->Eqn_k1[level][m][loc0+index] +=  w1*ww1*elt_K[pp*lms+qq]; /* direction 1 */
		    E->Eqn_k2[level][m][loc0+index] +=  w2*ww1*elt_K[(pp+1)*lms+qq]; /* direction 1 */
		    E->Eqn_k3[level][m][loc0+index] +=  w3*ww1*elt_K[(pp+2)*lms+qq]; /* direction 1 */

		     /* search for direction 2*/

		    found=0;
		    for(k=0;k<max_eqn;k++)
			if(E->Node_map[level][m][loc0+k] == eqn2) { /* found, index next equation */
			    index=k;
			    found++;
			    break;
			}

		    assert(found /* direction 2 */);

		    E->Eqn_k1[level][m][loc0+index] += w1*ww2*elt_K[pp*lms+qq+1]; /* direction 1 */
		    E->Eqn_k2[level][m][loc0+index] += w2*ww2*elt_K[(pp+1)*lms+qq+1]; /* direction 2 */
		    E->Eqn_k3[level][m][loc0+index] += w3*ww2*elt_K[(pp+2)*lms+qq+1]; /* direction 3 */

		    /* search for direction 3*/

                    found=0;
		    for(k=0;k<max_eqn;k++)
		    if(E->Node_map[level][m][loc0+k] == eqn3) { /* found, index next equation */
			index=k;
			found++;
			break;
		        }

                    assert(found /* direction 3 */);

		    E->Eqn_k1[level][m][loc0+index] += w1*ww3*elt_K[pp*lms+qq+2]; /* direction 1 */
                    E->Eqn_k2[level][m][loc0+index] += w2*ww3*elt_K[(pp+1)*lms+qq+2]; /* direction 2 */
		    E->Eqn_k3[level][m][loc0+index] += w3*ww3*elt_K[(pp+2)*lms+qq+2]; /* direction 3 */

		    }   /* end for j */
		  }   /* end for node1<= node */
		}      /* end for i */
	    }            /* end for element */
	}           /* end for m */

     (E->solver.exchange_id_d)(E, E->BI[level], level);
     if(E->control.ala_element_vanka_smoother)
         (E->solver.exchange_id_d)(E,E->ALA_vanka_overlap_BI[level],level);

     for(m=1;m<=E->sphere.caps_per_proc;m++)     {
        neq=E->lmesh.NEQ[level];

        for(j=0;j<neq;j++)                 {
            if(E->BI[level][m][j] ==0.0)  fprintf(stderr,"me= %d level %d, equation %d/%d has zero diagonal term\n",E->parallel.me,level,j,neq);
	    assert( E->BI[level][m][j] != 0 /* diagonal of matrix = 0, not acceptable */);
	    E->BI[level][m][j]  = (double) 1.0/E->BI[level][m][j];
	    if(E->control.ala_element_vanka_smoother) {
                if(E->ALA_vanka_overlap_BI[level][m][j] > 0.0) {
                    E->ALA_vanka_overlap_BI[level][m][j] =
                        1.0/E->ALA_vanka_overlap_BI[level][m][j];
                }
                else {
                    E->ALA_vanka_overlap_BI[level][m][j] = 0.0;
                }
            }
	    }
	}           /* end for m */


    }     /* end for level */

    return;
}


/* Cache a Cholesky factor of the complete local augmented velocity block.
   A single element stiffness has rigid-body null modes, so its diagonal is
   completed with the positive contribution from all other assembled elements
   incident on the same velocity dofs. */
void build_ala_element_vanka_factors(struct All_variables *E)
{
    int level,m,e,a,da,i,j,k,node,eq,fixed[ALA_VANKA_DOF];
    int coarse_level,coarse_e,coarse_x,coarse_y,coarse_z;
    int fine_x,fine_y,fine_z,fine_e;
    int local_elements,global_elements;
    double matrix[ALA_VANKA_DOF*ALA_VANKA_DOF];
    double L[ALA_VANKA_DOF*ALA_VANKA_DOF];
    double gradient[ALA_VANKA_DOF],forward[ALA_VANKA_DOF];
    double velocity_pressure[ALA_VANKA_DOF];
    double sum,pivot,maxdiag,global_diag,external,external_weight,shift;
    double schur,diagonal;
    double local_min_ratio,global_min_ratio,local_megabytes,global_max_megabytes;
    double build_start,build_seconds,global_build_seconds;
    higher_precision *cache;
    void get_elt_k();
    void get_ala_aug_k();
    void myerror();
    double CPU_time0();
    const int dims=E->mesh.nsd;
    const int ends=enodes[dims];
    const int n=loc_mat_size[dims];

    if(n!=ALA_VANKA_DOF)
        myerror(E,"Full ALA element-Vanka requires 3-D 24-dof elements");

    build_start=CPU_time0();
    local_elements=0;
    local_megabytes=0.0;
    local_min_ratio=1.0e300;
    for(level=E->mesh.gridmin;level<=E->mesh.gridmax;level++)
        for(m=1;m<=E->sphere.caps_per_proc;m++) {
            /* The coupled Schur correction consumes the finest-level patch
             * factors.  Coarser caches are also the production velocity-MG
             * smoother and must retain their assembled-diagonal completion.
             * Keeping one cache avoids another 43 MB/rank allocation. */
            external_weight=(level==E->mesh.gridmax)
                ? E->control.ala_element_vanka_external_diagonal_weight
                : 1.0;
            local_elements += E->lmesh.NEL[level];
            local_megabytes += (E->lmesh.NEL[level]+1)
                *(ALA_VANKA_CHOL_SIZE*sizeof(higher_precision)
                  +sizeof(double))/(1024.0*1024.0);
            for(e=1;e<=E->lmesh.NEL[level];e++) {
                get_elt_k(E,e,matrix,level,m,0);
                get_ala_aug_k(E,e,matrix,level,m);
                for(i=0;i<n*n;i++)
                    L[i]=0.0;
                for(a=1;a<=ends;a++) {
                    node=E->IEN[level][m][e].node[a];
                    for(da=0;da<dims;da++) {
                        i=(a-1)*dims+da;
                        fixed[i]=(da==0 && (E->NODE[level][m][node]&VBX)) ||
                                 (da==1 && (E->NODE[level][m][node]&VBY)) ||
                                 (da==2 && (E->NODE[level][m][node]&VBZ));
                    }
                }
                maxdiag=0.0;
                for(i=0;i<n;i++) {
                    a=i/dims+1;
                    da=i%dims;
                    node=E->IEN[level][m][e].node[a];
                    eq=E->ID[level][m][node].doff[da+1];
                    if(!isfinite(E->BI[level][m][eq]) ||
                       E->BI[level][m][eq]<=0.0)
                        myerror(E,"ALA element-Vanka assembled diagonal is invalid");
                    gradient[i]=fixed[i] ? 0.0 :
                        ALA_COMBINED_PRESSURE_COEFFICIENT(
                            E->elt_del[level][m][e].g[i][0],
                            E->elt_c[level][m][e].c[i][0]);
                    global_diag=1.0/E->BI[level][m][eq];
                    if(!fixed[i]) {
                        external=global_diag-matrix[i*n+i];
                        if(external < -1.0e-8*max(global_diag,1.0e-300)) {
                            fprintf(stderr,"rank=%d level=%d element=%d dof=%d "
                                    "global_diag=%e local_diag=%e\n",
                                    E->parallel.me,level,e,i,global_diag,
                                    matrix[i*n+i]);
                            myerror(E,"ALA element-Vanka external diagonal is negative");
                        }
                        matrix[i*n+i] +=
                            external_weight*max(external,0.0);
                        maxdiag=max(maxdiag,matrix[i*n+i]);
                    }
                }
                if(!isfinite(maxdiag) || maxdiag<=0.0)
                    myerror(E,"ALA element-Vanka local diagonal is invalid");
                shift=E->control.ala_element_vanka_regularization*maxdiag;
                for(i=0;i<n;i++) {
                    if(fixed[i]) {
                        for(j=0;j<n;j++)
                            matrix[i*n+j]=matrix[j*n+i]=0.0;
                        matrix[i*n+i]=maxdiag;
                    }
                    else {
                        matrix[i*n+i] += shift;
                        for(j=0;j<i;j++) {
                            sum=0.5*(matrix[i*n+j]+matrix[j*n+i]);
                            matrix[i*n+j]=matrix[j*n+i]=sum;
                        }
                    }
                }
                for(i=0;i<n;i++) {
                    for(j=0;j<=i;j++) {
                        sum=matrix[i*n+j];
                        for(k=0;k<j;k++)
                            sum -= L[i*n+k]*L[j*n+k];
                        if(i==j) {
                            pivot=sum;
                            if(!isfinite(pivot) ||
                               pivot<=1.0e-14*maxdiag) {
                                fprintf(stderr,"rank=%d level=%d element=%d "
                                        "pivot_dof=%d pivot=%e maxdiag=%e shift=%e\n",
                                        E->parallel.me,level,e,i,pivot,maxdiag,shift);
                                myerror(E,"ALA element-Vanka Cholesky failed");
                            }
                            L[i*n+j]=sqrt(pivot);
                            local_min_ratio=min(local_min_ratio,pivot/maxdiag);
                        }
                        else
                            L[i*n+j]=sum/L[j*n+j];
                    }
                }
                cache=E->ALA_vanka_chol[level][m]
                    +e*ALA_VANKA_CHOL_SIZE;
                k=0;
                for(i=0;i<n;i++)
                    for(j=0;j<=i;j++)
                        cache[k++]=(higher_precision)L[i*n+j];
                for(i=0;i<n;i++) {
                    sum=gradient[i];
                    for(j=0;j<i;j++)
                        sum -= L[i*n+j]*forward[j];
                    diagonal=L[i*n+i];
                    forward[i]=sum/diagonal;
                }
                for(i=n-1;i>=0;i--) {
                    sum=forward[i];
                    for(j=i+1;j<n;j++)
                        sum -= L[j*n+i]*velocity_pressure[j];
                    velocity_pressure[i]=sum/L[i*n+i];
                }
                schur=E->control.ala_element_vanka_regularization
                    *E->ECO[level][m][e].area;
                for(i=0;i<n;i++)
                    schur += gradient[i]*velocity_pressure[i];
                if(!isfinite(schur) || schur<=1.0e-300)
                    myerror(E,"ALA element-Vanka cached Schur is invalid");
                E->ALA_vanka_schur[level][m][e]=schur;
                E->ALA_vanka_valid[level][m][e]=1;
            }
        }

    /* Optional Galerkin pressure-Schur scale.  The local mixed patch remains
     * the same K/G solve, but its P0 Schur denominator on a coarse element is
     * inherited as Pp^T S_f Pp: the sum of its eight child denominators.  This
     * is a targeted Stage 3 diagnostic for the large Schur-coarsening defect;
     * it does not alter the physical operator or the velocity factors. */
    if(E->control.ala_element_vanka_galerkin_schur)
        /* Build from fine to coarse so every child Schur is already the
         * Galerkin value before it is aggregated into its parent. */
        for(level=E->mesh.gridmax;level>=E->mesh.gridmin+1;level--) {
            coarse_level=level-1;
            for(m=1;m<=E->sphere.caps_per_proc;m++)
                for(coarse_y=1;coarse_y<=E->lmesh.ELY[coarse_level];
                    coarse_y++)
                    for(coarse_x=1;coarse_x<=E->lmesh.ELX[coarse_level];
                        coarse_x++)
                        for(coarse_z=1;
                            coarse_z<=E->lmesh.ELZ[coarse_level];
                            coarse_z++) {
                            coarse_e=coarse_z
                                +(coarse_x-1)*E->lmesh.ELZ[coarse_level]
                                +(coarse_y-1)*E->lmesh.ELZ[coarse_level]
                                  *E->lmesh.ELX[coarse_level];
                            schur=0.0;
                            for(fine_y=2*coarse_y-1;fine_y<=2*coarse_y;
                                fine_y++)
                                for(fine_x=2*coarse_x-1;
                                    fine_x<=2*coarse_x;fine_x++)
                                    for(fine_z=2*coarse_z-1;
                                        fine_z<=2*coarse_z;fine_z++) {
                                        fine_e=fine_z
                                            +(fine_x-1)*E->lmesh.ELZ[level]
                                            +(fine_y-1)*E->lmesh.ELZ[level]
                                              *E->lmesh.ELX[level];
                                        schur +=
                                            E->ALA_vanka_schur[level][m]
                                              [fine_e];
                                    }
                            if(!isfinite(schur) || schur<=1.0e-300)
                                myerror(E,"Invalid Galerkin pressure Schur");
                            E->ALA_vanka_schur[coarse_level][m][coarse_e]
                                =schur;
                        }
        }
        if(E->control.ala_element_vanka_galerkin_schur &&
           E->parallel.me==0) {
            fprintf(E->fp,"ALA GALERKIN PRESSURE SCHUR SCALE enabled "
                    "operator=PpT_Sfine_Pp action=coarse_child_sum "
                    "velocity_factors=rediscretized observe_velocity=on\n");
            fprintf(stderr,"ALA GALERKIN PRESSURE SCHUR SCALE enabled "
                    "operator=PpT_Sfine_Pp action=coarse_child_sum "
                    "velocity_factors=rediscretized observe_velocity=on\n");
        }
    MPI_Allreduce(&local_elements,&global_elements,1,MPI_INT,MPI_SUM,
                  E->parallel.world);
    MPI_Allreduce(&local_min_ratio,&global_min_ratio,1,MPI_DOUBLE,MPI_MIN,
                  E->parallel.world);
    MPI_Allreduce(&local_megabytes,&global_max_megabytes,1,MPI_DOUBLE,MPI_MAX,
                  E->parallel.world);
    build_seconds=CPU_time0()-build_start;
    MPI_Allreduce(&build_seconds,&global_build_seconds,1,MPI_DOUBLE,MPI_MAX,
                  E->parallel.world);
    if(E->parallel.me==0) {
        fprintf(E->fp,"ALA full element-Vanka factors levels=%d "
                "global_elements=%d max_cache_per_rank_mb=%g "
                "regularization=%e min_pivot_ratio=%e "
                "finest_external_diagonal_weight=%e "
                "coarse_external_diagonal_weight=1.000000e+00 "
                "build_seconds_max=%e cycle=%d\n",
                E->mesh.gridmax-E->mesh.gridmin+1,global_elements,
                global_max_megabytes,
                E->control.ala_element_vanka_regularization,global_min_ratio,
                E->control.ala_element_vanka_external_diagonal_weight,
                global_build_seconds,E->monitor.solution_cycles);
        fprintf(stderr,"ALA full element-Vanka factors levels=%d "
                "global_elements=%d max_cache_per_rank_mb=%g "
                "regularization=%e min_pivot_ratio=%e "
                "finest_external_diagonal_weight=%e "
                "coarse_external_diagonal_weight=1.000000e+00 "
                "build_seconds_max=%e cycle=%d\n",
                E->mesh.gridmax-E->mesh.gridmin+1,global_elements,
                global_max_megabytes,
                E->control.ala_element_vanka_regularization,global_min_ratio,
                E->control.ala_element_vanka_external_diagonal_weight,
                global_build_seconds,E->monitor.solution_cycles);
        fflush(E->fp);
        fflush(stderr);
    }

    if(E->control.ala_viscosity_spectrum_diagnostics &&
       E->monitor.solution_cycles%
         E->control.ala_viscosity_spectrum_interval==0) {
        int radial,global_radial,v,radial_count;
        int *local_eta_count,*global_eta_count;
        int *local_diag_count,*global_diag_count;
        int *local_schur_count,*global_schur_count;
        double eta,diag,depth;
        double *local_eta_min,*global_eta_min,*local_eta_max,*global_eta_max;
        double *local_eta_sum,*global_eta_sum;
        double *local_diag_min,*global_diag_min,*local_diag_max,*global_diag_max;
        double *local_diag_sum,*global_diag_sum;
        double *local_schur_min,*global_schur_min;
        double *local_schur_max,*global_schur_max;
        double *local_schur_sum,*global_schur_sum;
        double *local_depth_sum,*global_depth_sum;
        const int vpts=vpoints[dims];

        for(level=E->mesh.gridmin;level<=E->mesh.gridmax;level++) {
            radial_count=E->mesh.ELZ[level];
#define ALA_ALLOC_DIAG(name,type) \
            name=(type *)malloc(radial_count*sizeof(type)); \
            if(name==NULL) myerror(E,"Unable to allocate ALA viscosity diagnostics")
            ALA_ALLOC_DIAG(local_eta_count,int);
            ALA_ALLOC_DIAG(global_eta_count,int);
            ALA_ALLOC_DIAG(local_diag_count,int);
            ALA_ALLOC_DIAG(global_diag_count,int);
            ALA_ALLOC_DIAG(local_schur_count,int);
            ALA_ALLOC_DIAG(global_schur_count,int);
            ALA_ALLOC_DIAG(local_eta_min,double);
            ALA_ALLOC_DIAG(global_eta_min,double);
            ALA_ALLOC_DIAG(local_eta_max,double);
            ALA_ALLOC_DIAG(global_eta_max,double);
            ALA_ALLOC_DIAG(local_eta_sum,double);
            ALA_ALLOC_DIAG(global_eta_sum,double);
            ALA_ALLOC_DIAG(local_diag_min,double);
            ALA_ALLOC_DIAG(global_diag_min,double);
            ALA_ALLOC_DIAG(local_diag_max,double);
            ALA_ALLOC_DIAG(global_diag_max,double);
            ALA_ALLOC_DIAG(local_diag_sum,double);
            ALA_ALLOC_DIAG(global_diag_sum,double);
            ALA_ALLOC_DIAG(local_schur_min,double);
            ALA_ALLOC_DIAG(global_schur_min,double);
            ALA_ALLOC_DIAG(local_schur_max,double);
            ALA_ALLOC_DIAG(global_schur_max,double);
            ALA_ALLOC_DIAG(local_schur_sum,double);
            ALA_ALLOC_DIAG(global_schur_sum,double);
            ALA_ALLOC_DIAG(local_depth_sum,double);
            ALA_ALLOC_DIAG(global_depth_sum,double);
#undef ALA_ALLOC_DIAG
            for(radial=0;radial<radial_count;radial++) {
                local_eta_count[radial]=local_diag_count[radial]=0;
                local_schur_count[radial]=0;
                local_eta_min[radial]=local_diag_min[radial]=1.0e300;
                local_schur_min[radial]=1.0e300;
                local_eta_max[radial]=local_diag_max[radial]=0.0;
                local_schur_max[radial]=0.0;
                local_eta_sum[radial]=local_diag_sum[radial]=0.0;
                local_schur_sum[radial]=local_depth_sum[radial]=0.0;
            }
            for(m=1;m<=E->sphere.caps_per_proc;m++)
                for(e=1;e<=E->lmesh.NEL[level];e++) {
                    radial=(e-1)%E->lmesh.ELZ[level];
                    global_radial=E->lmesh.EZS[level]+radial;
                    depth=(1.0-E->ECO[level][m][e].centre[3])
                        *E->data.radius_km;
                    local_depth_sum[global_radial] += depth;
                    local_schur_count[global_radial]++;
                    schur=E->ALA_vanka_schur[level][m][e];
                    local_schur_min[global_radial]=
                        min(local_schur_min[global_radial],schur);
                    local_schur_max[global_radial]=
                        max(local_schur_max[global_radial],schur);
                    local_schur_sum[global_radial] += schur;
                    for(v=1;v<=vpts;v++) {
                        eta=E->EVI[level][m][(e-1)*vpts+v];
                        local_eta_min[global_radial]=
                            min(local_eta_min[global_radial],eta);
                        local_eta_max[global_radial]=
                            max(local_eta_max[global_radial],eta);
                        local_eta_sum[global_radial] += eta;
                        local_eta_count[global_radial]++;
                    }
                    for(i=0;i<n;i++) {
                        a=i/dims+1;
                        da=i%dims;
                        node=E->IEN[level][m][e].node[a];
                        eq=E->ID[level][m][node].doff[da+1];
                        diag=1.0/E->BI[level][m][eq];
                        local_diag_min[global_radial]=
                            min(local_diag_min[global_radial],diag);
                        local_diag_max[global_radial]=
                            max(local_diag_max[global_radial],diag);
                        local_diag_sum[global_radial] += diag;
                        local_diag_count[global_radial]++;
                    }
                }
#define ALA_REDUCE_DIAG(local,global,type,operation) \
            MPI_Allreduce(local,global,radial_count,type,operation,E->parallel.world)
            ALA_REDUCE_DIAG(local_eta_count,global_eta_count,MPI_INT,MPI_SUM);
            ALA_REDUCE_DIAG(local_diag_count,global_diag_count,MPI_INT,MPI_SUM);
            ALA_REDUCE_DIAG(local_schur_count,global_schur_count,MPI_INT,MPI_SUM);
            ALA_REDUCE_DIAG(local_eta_min,global_eta_min,MPI_DOUBLE,MPI_MIN);
            ALA_REDUCE_DIAG(local_eta_max,global_eta_max,MPI_DOUBLE,MPI_MAX);
            ALA_REDUCE_DIAG(local_eta_sum,global_eta_sum,MPI_DOUBLE,MPI_SUM);
            ALA_REDUCE_DIAG(local_diag_min,global_diag_min,MPI_DOUBLE,MPI_MIN);
            ALA_REDUCE_DIAG(local_diag_max,global_diag_max,MPI_DOUBLE,MPI_MAX);
            ALA_REDUCE_DIAG(local_diag_sum,global_diag_sum,MPI_DOUBLE,MPI_SUM);
            ALA_REDUCE_DIAG(local_schur_min,global_schur_min,MPI_DOUBLE,MPI_MIN);
            ALA_REDUCE_DIAG(local_schur_max,global_schur_max,MPI_DOUBLE,MPI_MAX);
            ALA_REDUCE_DIAG(local_schur_sum,global_schur_sum,MPI_DOUBLE,MPI_SUM);
            ALA_REDUCE_DIAG(local_depth_sum,global_depth_sum,MPI_DOUBLE,MPI_SUM);
#undef ALA_REDUCE_DIAG
            if(E->parallel.me==0)
                for(radial=radial_count-1;radial>=0;radial--) {
                    if(global_schur_count[radial]<=0 ||
                       global_eta_count[radial]<=0 ||
                       global_diag_count[radial]<=0)
                        myerror(E,"ALA viscosity diagnostic radial layer is empty");
                    fprintf(E->fp,"ALA VISCOSITY SPECTRUM level=%d "
                            "radial_element=%d depth_km=%e "
                            "eta_range=[%e,%e] eta_mean=%e "
                            "diagK_range=[%e,%e] diagK_mean=%e "
                            "local_schur_range=[%e,%e] "
                            "local_schur_mean=%e "
                            "eta_mean_times_schur_mean=%e "
                            "cycle=%d\n",
                            level,radial,
                            global_depth_sum[radial]
                              /global_schur_count[radial],
                            global_eta_min[radial],global_eta_max[radial],
                            global_eta_sum[radial]/global_eta_count[radial],
                            global_diag_min[radial],global_diag_max[radial],
                            global_diag_sum[radial]/global_diag_count[radial],
                            global_schur_min[radial],
                            global_schur_max[radial],
                            global_schur_sum[radial]
                              /global_schur_count[radial],
                            global_eta_sum[radial]/global_eta_count[radial]
                              *global_schur_sum[radial]
                              /global_schur_count[radial],
                            E->monitor.solution_cycles);
                }
#define ALA_FREE_DIAG(name) free(name)
            ALA_FREE_DIAG(local_eta_count); ALA_FREE_DIAG(global_eta_count);
            ALA_FREE_DIAG(local_diag_count); ALA_FREE_DIAG(global_diag_count);
            ALA_FREE_DIAG(local_schur_count); ALA_FREE_DIAG(global_schur_count);
            ALA_FREE_DIAG(local_eta_min); ALA_FREE_DIAG(global_eta_min);
            ALA_FREE_DIAG(local_eta_max); ALA_FREE_DIAG(global_eta_max);
            ALA_FREE_DIAG(local_eta_sum); ALA_FREE_DIAG(global_eta_sum);
            ALA_FREE_DIAG(local_diag_min); ALA_FREE_DIAG(global_diag_min);
            ALA_FREE_DIAG(local_diag_max); ALA_FREE_DIAG(global_diag_max);
            ALA_FREE_DIAG(local_diag_sum); ALA_FREE_DIAG(global_diag_sum);
            ALA_FREE_DIAG(local_schur_min); ALA_FREE_DIAG(global_schur_min);
            ALA_FREE_DIAG(local_schur_max); ALA_FREE_DIAG(global_schur_max);
            ALA_FREE_DIAG(local_schur_sum); ALA_FREE_DIAG(global_schur_sum);
            ALA_FREE_DIAG(local_depth_sum); ALA_FREE_DIAG(global_depth_sum);
#undef ALA_FREE_DIAG
        }
        if(E->parallel.me==0)
            fflush(E->fp);
    }
}

void rebuild_BI_on_boundary(E)
     struct All_variables *E;
{
    int m,level,i,j;
    int eqn1,eqn2,eqn3;

    higher_precision *B1,*B2,*B3;
    int *C;

    const int dims=E->mesh.nsd,dofs=E->mesh.dof;

    const int max_eqn = dims*14;

   for(level=E->mesh.gridmax;level>=E->mesh.gridmin;level--)   {
     for (m=1;m<=E->sphere.caps_per_proc;m++)  {
        for(j=0;j<=E->lmesh.NEQ[level];j++)
            E->temp[m][j]=0.0;

        for(i=1;i<=E->lmesh.NNO[level];i++)  {
            eqn1=E->ID[level][m][i].doff[1];
            eqn2=E->ID[level][m][i].doff[2];
            eqn3=E->ID[level][m][i].doff[3];

            C=E->Node_map[level][m] + (i-1)*max_eqn;
            B1=E->Eqn_k1[level][m]+(i-1)*max_eqn;
            B2=E->Eqn_k2[level][m]+(i-1)*max_eqn;
            B3=E->Eqn_k3[level][m]+(i-1)*max_eqn;

            for(j=3;j<max_eqn;j++) {
                E->temp[m][eqn1] += fabs(B1[j]);
                E->temp[m][eqn2] += fabs(B2[j]);
                E->temp[m][eqn3] += fabs(B3[j]);
                }

            for(j=0;j<max_eqn;j++)
                E->temp[m][C[j]] += fabs(B1[j]) + fabs(B2[j]) + fabs(B3[j]);

            }
        }

     (E->solver.exchange_id_d)(E, E->temp, level);

     for (m=1;m<=E->sphere.caps_per_proc;m++)  {
        for(i=0;i<E->lmesh.NEQ[level];i++)  {
            E->temp[m][i] = E->temp[m][i] - 1.0/E->BI[level][m][i];
            }
        for(i=1;i<=E->lmesh.NNO[level];i++)
          if (E->NODE[level][m][i] & OFFSIDE)   {
            eqn1=E->ID[level][m][i].doff[1];
            eqn2=E->ID[level][m][i].doff[2];
            eqn3=E->ID[level][m][i].doff[3];
            E->BI[level][m][eqn1] = (double) 1.0/E->temp[m][eqn1];
            E->BI[level][m][eqn2] = (double) 1.0/E->temp[m][eqn2];
            E->BI[level][m][eqn3] = (double) 1.0/E->temp[m][eqn3];
            }
        }


    }     /* end for level */

 return;
}


/* ============================================
   Function to set up the boundary condition
   masks and other indicators.
   ============================================  */

void construct_masks(E)		/* Add lid/edge masks/nodal weightings */
     struct All_variables *E;
{
  int i,j,k,l,node,el,elt;
  int lev,elx,elz,ely,nno,nox,noz,noy;

  for(lev=E->mesh.gridmax;lev>=E->mesh.gridmin;lev--)
    for (j=1;j<=E->sphere.caps_per_proc;j++)           {
      elz = E->lmesh.ELZ[lev];
      ely = E->lmesh.ELY[lev];
      noy = E->lmesh.NOY[lev];
      noz = E->lmesh.NOZ[lev];
      nno = E->lmesh.NNO[lev];

        if (E->parallel.me_loc[3]==0 )
          for (i=1;i<=E->parallel.NUM_NNO[lev][j].bound[5];i++)   {
            node = E->parallel.NODE[lev][j][i].bound[5];
 	    E->NODE[lev][j][node] = E->NODE[lev][j][node] | TZEDGE;
	    }
        if ( E->parallel.me_loc[3]==E->parallel.nprocz-1 )
          for (i=1;i<=E->parallel.NUM_NNO[lev][j].bound[6];i++)   {
  	    node = E->parallel.NODE[lev][j][i].bound[6];
	    E->NODE[lev][j][node] = E->NODE[lev][j][node] | TZEDGE;
	    }

      }    /* end for j & lev */

/*   if (E->control.verbose) { */
/*     for(lev=E->mesh.gridmax;lev>=E->mesh.gridmin;lev--)  */
/*       for (j=1;j<=E->sphere.caps_per_proc;j++)           { */
/*         for (i=1;i<=E->parallel.NUM_NNO[lev][j].bound[5];i++)   {  */
/* 	  node = E->parallel.NODE[lev][j][i].bound[5]; */
/* 	  fprintf(E->fp_out,"bound=5  NODE[lev=%1d][node=%3d]=%d\n",lev,node,E->NODE[lev][j][node]); */
/* 	} */
/*         for (i=1;i<=E->parallel.NUM_NNO[lev][j].bound[6];i++)   {  */
/* 	  node = E->parallel.NODE[lev][j][i].bound[6]; */
/* 	  fprintf(E->fp_out,"bound=6  NODE[lev=%1d][node=%3d]=%d\n",lev,node,E->NODE[lev][j][node]); */
/* 	} */
/*       } */
/*     fflush(E->fp_out); */
/*   } */

  return;
  }


/*   ==========================================
     build the sub-element reference matrices
     ==========================================   */

void construct_sub_element(E)
     struct All_variables *E;

{    int i,j,k,l,m;
     int lev,nox,noy,noz,nnn,elx,elz,ely,elzu,elxu,elt,eltu;


  for(lev=E->mesh.levmax-1;lev>=E->mesh.levmin;lev--)
     for (m=1;m<=E->sphere.caps_per_proc;m++)       {
          elx = E->lmesh.ELX[lev];
	  elz = E->lmesh.ELZ[lev];
	  ely = E->lmesh.ELY[lev];
          nox = E->lmesh.NOX[lev];
          noy = E->lmesh.NOY[lev];
          noz = E->lmesh.NOZ[lev];
	  elz = E->lmesh.ELZ[lev];
	  ely = E->lmesh.ELY[lev];
	  elxu = 2 * elx;
	  elzu = 2 * elz;
          if (!(E->control.NMULTIGRID||E->control.EMULTIGRID))  {
             elzu = 1;
             if (lev == E->mesh.levmax-1)
                 elzu = E->lmesh.ELZ[E->mesh.levmax];
             }

	  for(i=1;i<=elx;i++)
	    for(j=1;j<=elz;j++)
	      for(k=1;k<=ely;k++)    {
		  elt = j + (i-1)*elz +(k-1)*elz*elx;
		  eltu = (j*2-1) + elzu *2*(i-1) + elxu*elzu*2*(k-1);

		  for(l=1;l<=enodes[E->mesh.nsd];l++)   {
		      E->EL[lev][m][elt].sub[l] = eltu
                                 + offset[l].vector[0]
                                 + offset[l].vector[1] * elzu
                                 + offset[l].vector[2] * elzu * elxu;
		      }
		  }

	  }


   return;
   }


void construct_elt_ks(E)
     struct All_variables *E;
{
    int e,el,lev,j,k,ii,m;
    void get_elt_k();
    void get_aug_k();
    void get_ala_aug_k();
    void build_diagonal_of_K();

    const int dims=E->mesh.nsd;
    const int n=loc_mat_size[E->mesh.nsd];

/*     if(E->parallel.me==0) */
/* 	fprintf(stderr,"storing elt k matrices\n"); */

    for(lev=E->mesh.gridmin;lev<=E->mesh.gridmax;lev++)  {

      for(m=1;m<=E->sphere.caps_per_proc;m++)     {

	for(el=1;el<=E->lmesh.NEL[lev];el++)    {

	    get_elt_k(E,el,E->elt_k[lev][m][el].k,lev,m,0);

	    if (E->control.augmented_Lagr)
	        get_aug_k(E,el,E->elt_k[lev][m][el].k,lev,m);

            if(E->control.ala_augmented_lagrangian_gamma > 0.0)
                get_ala_aug_k(E,el,E->elt_k[lev][m][el].k,lev,m);

            build_diagonal_of_K(E,el,E->elt_k[lev][m][el].k,lev,m);

	    }
	}        /* end for m */

      (E->solver.exchange_id_d)(E, E->BI[lev], lev);    /*correct BI   */

      for(m=1;m<=E->sphere.caps_per_proc;m++)

            for(j=0;j<E->lmesh.NEQ[lev];j++) {
	       if(E->BI[lev][m][j] ==0.0)  fprintf(stderr,"me= %d level %d, equation %d/%d has zero diagonal term\n",E->parallel.me,lev,j,E->mesh.NEQ[lev]);
               assert( E->BI[lev][m][j] != 0 /* diagonal of matrix = 0, not acceptable */);
               E->BI[lev][m][j]  = (float) 1.0/E->BI[lev][m][j];
	       }

    }       /* end for level */

  return;
}



void construct_elt_gs(E)
     struct All_variables *E;
{ int m,el,lev,a;
  void get_elt_g();

  const int dims=E->mesh.nsd,dofs=E->mesh.dof;
  const int ends=enodes[dims];

/*   if(E->control.verbose && E->parallel.me==0) */
/*       fprintf(stderr,"storing elt g matrices\n"); */

  for(lev=E->mesh.gridmin;lev<=E->mesh.gridmax;lev++)
    for(m=1;m<=E->sphere.caps_per_proc;m++)
      for(el=1;el<=E->lmesh.NEL[lev];el++)
        get_elt_g(E,el,E->elt_del[lev][m][el].g,lev,m);


  return;
}


/*==============================================
  For compressible cases, construct c matrix,
  where  c = \frac{d rho_r}{dr} / rho_r * u_r
  ==============================================*/

void construct_elt_cs(struct All_variables *E)
{
    int m, el, lev;
    void get_elt_c();

/*     if(E->control.verbose && E->parallel.me==0) */
/*         fprintf(stderr,"storing elt c matrices\n"); */

    for(lev=E->mesh.gridmin;lev<=E->mesh.gridmax;lev++)
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(el=1;el<=E->lmesh.NEL[lev];el++) {
                get_elt_c(E,el,E->elt_c[lev][m][el].c,lev,m);
            }


    return;
}


/* ==============================================================
 routine for constructing stiffness and node_maps
 ============================================================== */

void construct_stiffness_B_matrix(E)
  struct All_variables *E;
{
  void build_diagonal_of_K();
  void build_diagonal_of_Ahat();
  void build_radial_line_Ahat_preconditioner();
  void project_viscosity();
  void construct_node_maps();
  void construct_node_ks();
  void build_ala_element_vanka_factors(struct All_variables *);
  void construct_elt_ks();
  void rebuild_BI_on_boundary();
  int lev,m,j;

  if (E->control.NMULTIGRID)
    project_viscosity(E);

  if (E->control.NMULTIGRID || E->control.NASSEMBLE) {
    construct_node_ks(E);
    if(E->control.ala_element_vanka_smoother) {
      int rebuild_vanka=1;
      int factor_age=E->monitor.solution_cycles
          -E->control.ala_element_vanka_last_build_cycle;

      /* These factors are a preconditioner, not the physical operator.
         For temperature-dependent viscosity the matrix changes every step,
         but FGMRES still audits and converges the current assembled operator.
         Reusing a recent factorization therefore changes cost and convergence
         only.  The default interval of one preserves the historical path. */
      if(E->control.ala_coupled_element_vanka &&
         E->control.ala_element_vanka_rebuild_interval>1 &&
         E->control.ala_element_vanka_last_build_cycle>=0 &&
         factor_age>=0 &&
         factor_age<E->control.ala_element_vanka_rebuild_interval)
        rebuild_vanka=0;
      if(rebuild_vanka) {
        build_ala_element_vanka_factors(E);
        E->control.ala_element_vanka_last_build_cycle=
            E->monitor.solution_cycles;
      }
      else if(E->parallel.me==0) {
        fprintf(E->fp,"ALA full element-Vanka factors action=reuse "
                "source_cycle=%d current_cycle=%d age=%d interval=%d "
                "scope=preconditioner_only current_operator=reassembled\n",
                E->control.ala_element_vanka_last_build_cycle,
                E->monitor.solution_cycles,factor_age,
                E->control.ala_element_vanka_rebuild_interval);
        fprintf(stderr,"ALA full element-Vanka factors action=reuse "
                "source_cycle=%d current_cycle=%d age=%d interval=%d "
                "scope=preconditioner_only current_operator=reassembled\n",
                E->control.ala_element_vanka_last_build_cycle,
                E->monitor.solution_cycles,factor_age,
                E->control.ala_element_vanka_rebuild_interval);
        fflush(E->fp);
        fflush(stderr);
      }
    }
  }
  else {
    construct_elt_ks(E);
  }

  /* Preserve the positive assembled inverse velocity diagonal before the
     multigrid boundary reconstruction mutates BI ghost entries.  The strict
     ALA coarse pressure operator needs this fixed SPD diagonal. */
  for(lev=E->mesh.gridmin;lev<=E->mesh.gridmax;lev++)
    for(m=1;m<=E->sphere.caps_per_proc;m++)
      for(j=0;j<E->lmesh.NEQ[lev];j++)
        E->ALA_velocity_BI[lev][m][j]=E->BI[lev][m][j];

  /* BPI requires the positive assembled Jacobi inverse built above.  The
     subsequent offside-node BI reconstruction can contain nonpositive ghost
     entries and is only for the velocity multigrid path. */
  build_diagonal_of_Ahat(E);
  if(E->control.ala_radial_line_preconditioner)
    build_radial_line_Ahat_preconditioner(E);

  if (E->control.NMULTIGRID || (E->control.NASSEMBLE && !E->control.CONJ_GRAD))
    rebuild_BI_on_boundary(E);


  return;
}

/* took this apart to allow call from other subroutines */

/* 


determine layer number based on radial coordinate r

if E->viscosity.z... set to Earth values, then

1: lithosphere
2: 100-410
3: 410-660
4: lower mantle

*/
int layers_r(E,r)
     struct All_variables *E;
     float r;
{
  float rlith, r410, rlm;

  int llayers = 0;
  /* 
     the z-values, as read in, are non-dimensionalized depth
     convert to radii

  */
  rlith = E->sphere.ro - E->viscosity.zlith; /* lithosphere */
  r410  = E->sphere.ro - E->viscosity.z410;
  rlm   = E->sphere.ro - E->viscosity.zlm;

  if (r > rlith)		/* in lithospherre */
    llayers = 1;
  else if ((r > r410)&& (r  <= rlith)) /* in asthenosphere 100...410 km */
    llayers = 2;
  else if ((r > rlm) && (r <= r410)) /* in transition zone, 410 - 660 km */
    llayers = 3;
  else				/* lower mantle */
    llayers = 4;
  return (llayers);
}

/* determine layer number of node "node" of cap "m" */
int layers(E,m,node)
     struct All_variables *E;
     int m,node;
{
  return(layers_r(E,E->sx[m][3][node]));
}


/* ==============================================================
 construct array mat




 ============================================================== */
void construct_mat_group(E)
     struct All_variables *E;
{
  int m,i,j,k,kk,el,lev,a,nodea,els,llayer;

  const int dims=E->mesh.nsd,dofs=E->mesh.dof;
  const int ends=enodes[dims];

  for (m=1;m<=E->sphere.caps_per_proc;m++)   {
    for(el=1;el<=E->lmesh.nel;el++) {
      E->mat[m][el] = 1;
      nodea = E->ien[m][el].node[2];
      llayer = layers(E,m,nodea);
      if (llayer)  {
	E->mat[m][el] = llayer;
      }
    }
  }

  return;
}
