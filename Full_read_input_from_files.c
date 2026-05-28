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

void create_new_trench(struct All_variables *E, char *file_name, float *flag_depth, double maxdist);
void create_new_trench2(struct All_variables *E, char *file_name, float *flag_depth, double maxdist);
void mindst(double **sz, int *trench, int num, double x, double y, double *min, int *index);
void mindst2(double **sz, int num, double x, double y, double *min, int *index);
double cross_product(double x1, double y1, double x2, double y2, double x, double y);
double dot_product(double x1, double y1, double x2, double y2, double x, double y);

/*=======================================================================
  Calculate ages (MY) for opening input files -> material, ages, velocities
  Open these files, read in results, and average if necessary
=========================================================================*/

void full_read_input_files_for_timesteps(E,action,output)
    struct All_variables *E;
    int action, output;
{
    float find_age_in_MY();

    FILE *fp1, *fp2, *fp3,*fp4,*fp5,*fp6,*fp7;
    float age, newage1, newage2, inputdepth1,inputdepth2,lon,lat,flagdepth1,flagdepth2,tfdepth;
    char output_file1[255],output_file2[255],output_file3[255],output_file4[255],output_file5[255],output_file6[255],output_file7[255];
    float *VB1[4],*VB2[4], inputage1, inputage2, rad, theta, phi, mindist;
    int nox,noz,noy,nnn,nox1,noz1,noy1,lev;
    int i,ii,jj,kk,ll,m,mm,j,k,n,nodeg,nodel,node,cap;
    int intage, pos_age;
    int nodea;
    int nn,el;
    double temp;
    void create_new_trench();
    void create_new_trench2();

    const int dims=E->mesh.nsd;

    int elx,ely,elz,elg,emax;
    float *VIP1,*VIP2;
    int *LL1, *LL2;

    int llayer;
    int layers();

    nox=E->mesh.nox;
    noy=E->mesh.noy;
    noz=E->mesh.noz;
    nox1=E->lmesh.nox;
    noz1=E->lmesh.noz;
    noy1=E->lmesh.noy;
    lev=E->mesh.levmax;

    elx=E->lmesh.elx;
    elz=E->lmesh.elz;
    ely=E->lmesh.ely;
   
    age=find_age_in_MY(E);

    if (E->monitor.solution_cycles == 0){
	E->velo_1=(float*) malloc((nox*noy+1)*sizeof(float));
        E->velo_2=(float*) malloc((nox*noy+1)*sizeof(float));
        E->flag_depth=(float*) malloc((nox*noy+1)*sizeof(float));
	E->flag_depth1=(float*) malloc((nox*noy+1)*sizeof(float));
        E->flag_depth2=(float*) malloc((nox*noy+1)*sizeof(float));
        E->new_flag_depth=(float*) malloc((nox*noy+1)*sizeof(float));
        E->tf_depth=(float*) malloc((nox*noy+1)*sizeof(float));
	E->trench_visit_age=age+1;
    }

    emax=E->mesh.elx*E->mesh.elz*E->mesh.ely;

    if (age < 0.0) { /* age is negative -> use age=0 for input files */
      intage = 0;
      newage2 = newage1 = 0.0;
      pos_age = 0;
      exit(8);
    }
    else {
      intage = age;
      newage1 = 1.0*intage;
      newage2 = 1.0*intage + 1.0;
      pos_age = 1;
    }

    for (m=1;m<=E->sphere.caps_per_proc;m++)  {
      cap = E->sphere.capid[m] - 1;  /* capid: 1-12 */

      switch (action) { /* set up files to open */

      case 1:  /* read velocity boundary conditions */
	sprintf(output_file1,"%s%0.0f.%d",E->control.velocity_boundary_file,newage1,cap);
	sprintf(output_file2,"%s%0.0f.%d",E->control.velocity_boundary_file,newage2,cap);
	fp1=fopen(output_file1,"r");
	if (fp1 == NULL) {
          fprintf(E->fp,"(Problem_related #4) Cannot open %s\n",output_file1);
          exit(8);
	}
	if (pos_age) {
	  fp2=fopen(output_file2,"r");
	  if (fp2 == NULL) {
	    fprintf(E->fp,"(Problem_related #5) Cannot open %s\n",output_file2);
	    exit(8);
	  }
	}
	if((E->parallel.me==0) && (output==1))   {
	  fprintf(E->fp,"Velocity: Starting Age = %g, Elapsed time = %g, Current Age = %g\n",E->control.start_age,E->monitor.elapsed_time,age);
	  fprintf(E->fp,"Velocity: File1 = %s\n",output_file1);
	  if (pos_age)
	    fprintf(E->fp,"Velocity: File2 = %s\n",output_file2);
	  else
	    fprintf(E->fp,"Velocity: File2 = No file inputted (negative age)\n");
	}
	break;

      case 2:  /* read ages for lithosphere temperature assimilation */
	sprintf(output_file1,"%s%0.0f.%d",E->control.lith_age_file,newage1,cap);
	sprintf(output_file2,"%s%0.0f.%d",E->control.lith_age_file,newage2,cap);
	sprintf(output_file3,"%s%0.0f.xyz",E->control.flag_depth_new_file,newage1);
        sprintf(output_file4,"%s%0.0f.xyz",E->control.flag_depth_file,newage1);
        sprintf(output_file5,"%s%0.0f.xyz",E->control.flag_depth_file,newage2);
        sprintf(output_file6,"%s%0.0f.xyz",E->control.tf_file,newage1);
	fp1=fopen(output_file1,"r");
	if (fp1 == NULL) {
          fprintf(E->fp,"(Problem_related #6) Cannot open %s\n",output_file1);
          exit(8);
	}
	if (pos_age) {
	  fp2=fopen(output_file2,"r");
	  if (fp2 == NULL) {
	    fprintf(E->fp,"(Problem_related #7) Cannot open %s\n",output_file2);
	    exit(8);
	  }
	}
	/*fp6=fopen(output_file6,"r");
        if (fp6==NULL){
          fprintf(E->fp,"(Problem_related #12) Cannot open %s\n",output_file6);
          exit(8);
        }*/
	if((E->parallel.me==0) && (output==1))   {
	//if (output==1){
	  fprintf(E->fp,"Age: Starting Age = %g, Elapsed time = %g, Current Age = %g\n",E->control.start_age,E->monitor.elapsed_time,age);
	  fprintf(E->fp,"Age: File1 = %s\n",output_file1);
	  fprintf(E->fp,"Age: cap=%d,File1 = %s\n",E->sphere.capid[m],output_file1);
          fprintf(E->fp,"Age: cap=%d,File3 = %s\n",E->sphere.capid[m],output_file3);
          fprintf(E->fp,"Age: cap=%d,File5 = %s\n",E->sphere.capid[m],output_file5);
          //fprintf(E->fp,"Age: cap=%d,File7 = %s\n",E->sphere.capid[m],output_file7);
	  if (pos_age){
	    fprintf(E->fp,"Age: File2 = %s\n",output_file2);
	    fprintf(E->fp,"Age: cap=%d,File4 = %s\n",E->sphere.capid[m],output_file4);
            fprintf(E->fp,"Age: cap=%d,File6 = %s\n",E->sphere.capid[m],output_file6);
	  }
	  else
	    fprintf(E->fp,"Age: File2 = No file inputted (negative age)\n");
	}
	break;

      case 3:  /* read element materials */

	sprintf(output_file1,"%s%0.0f.%d",E->control.mat_file,newage1,cap);
	sprintf(output_file2,"%s%0.0f.%d",E->control.mat_file,newage2,cap);
	fp1=fopen(output_file1,"r");
	if (fp1 == NULL) {
          fprintf(E->fp,"(Problem_related #8) Cannot open %s\n",output_file1);
          exit(8);
	}
	if (pos_age) {
	  fp2=fopen(output_file2,"r");
	  if (fp2 == NULL) {
	    fprintf(E->fp,"(Problem_related #9) Cannot open %s\n",output_file2);
	    exit(8);
	  }
	}
	if((E->parallel.me==0) && (output==1))   {
	  fprintf(E->fp,"Mat: Starting Age = %g, Elapsed time = %g, Current Age = %g\n",E->control.start_age,E->monitor.elapsed_time,age);
	  fprintf(E->fp,"Mat: File1 = %s\n",output_file1);
	  if (pos_age)
	    fprintf(E->fp,"Mat: File2 = %s\n",output_file2);
	  else
	    fprintf(E->fp,"Mat: File2 = No file inputted (negative age)\n");
	}

	break;

      } /* end switch */



      switch (action) { /* Read the contents of files and average */

      case 1:  /* velocity boundary conditions */
	nnn=nox*noy;
	for(i=1;i<=dims;i++)  {
	  VB1[i]=(float*) malloc ((nnn+1)*sizeof(float));
	  VB2[i]=(float*) malloc ((nnn+1)*sizeof(float));
	}
	for(i=1;i<=nnn;i++)   {
	  fscanf(fp1,"%f %f",&(VB1[1][i]),&(VB1[2][i]));
	  VB1[1][i] *= E->data.timedir;
	  VB1[2][i] *= E->data.timedir;
	  E->velo_1[i] = VB1[1][i]*E->data.scalev;
	  E->velo_2[i] = VB1[1][i]*E->data.scalev;
	  if (pos_age) {
	    fscanf(fp2,"%f %f",&(VB2[1][i]),&(VB2[2][i]));
	    VB2[1][i] *= E->data.timedir;
	    VB2[2][i] *= E->data.timedir;
	  }
	  /* if( E->parallel.me ==0)
	     fprintf(stderr,"%d %f  %f  %f  %f\n",i,VB1[1][i],VB1[2][i],VB2[1][i],VB2[2][i]); */
	}
	fclose(fp1);
	if (pos_age) fclose(fp2);

	if(E->parallel.me_loc[3]==E->parallel.nprocz-1 )  {
          for(k=1;k<=noy1;k++)
	    for(i=1;i<=nox1;i++)    {
	      nodeg = E->lmesh.nxs+i-1 + (E->lmesh.nys+k-2)*nox;
	      nodel = (k-1)*nox1*noz1 + (i-1)*noz1+noz1;
	      if (pos_age) { /* positive ages - we must interpolate */
		E->sphere.cap[m].VB[1][nodel] = (VB1[1][nodeg] + (VB2[1][nodeg]-VB1[1][nodeg])/(newage2-newage1)*(age-newage1))*E->data.scalev;
		E->sphere.cap[m].VB[2][nodel] = (VB1[2][nodeg] + (VB2[2][nodeg]-VB1[2][nodeg])/(newage2-newage1)*(age-newage1))*E->data.scalev;
		E->sphere.cap[m].VB[3][nodel] = 0.0;
	      }
	      else { /* negative ages - don't do the interpolation */
		E->sphere.cap[m].VB[1][nodel] = VB1[1][nodeg] * E->data.scalev;
		E->sphere.cap[m].VB[2][nodel] = VB1[2][nodeg] * E->data.scalev;
		E->sphere.cap[m].VB[3][nodel] = 0.0;
	      }
	    }
	}   /* end of E->parallel.me_loc[3]==E->parallel.nproczl-1   */
	for(i=1;i<=dims;i++) {
          free ((void *) VB1[i]);
          free ((void *) VB2[i]);
	}
	break;

      case 2:  /* ages for lithosphere temperature assimilation */
	for(i=1;i<=noy;i++)
	  for(j=1;j<=nox;j++) {
	    node=j+(i-1)*nox;
	    fscanf(fp1,"%f",&inputage1);
            /*fscanf(fp6,"%f %f %f",&lon, &lat, &tfdepth);*/
	    if (pos_age) { /* positive ages - we must interpolate */
              fscanf(fp2,"%f",&inputage2);
              E->age_t[node] = (inputage1 + (inputage2-inputage1)/(newage2-newage1)*(age-newage1))/E->data.scalet;
              /*E->tf_depth[node] = tfdepth;*/
	    }
	    else { /* negative ages - don't do the interpolation */
              E->age_t[node] = inputage1/E->data.scalet;
              /*E->tf_depth[node] = tfdepth;*/
	    }
	  }
	fclose(fp1);
        /*fclose(fp6);*/
	if (pos_age){
		 fclose(fp2);
	}
	if(intage<E->trench_visit_age && E->monitor.solution_cycles>0) {
            if(E->parallel.me==0) fprintf(stderr,"age=%f\n",age);
            create_new_trench2(E,output_file3,E->new_flag_depth, 500.0);
            create_new_trench2(E,output_file4,E->flag_depth, 1000.0);
            create_new_trench2(E,output_file5,E->flag_depth1, 1000.0);
	    for(k=1;k<=noy1;k++)
                for(i=1;i<=nox1;i++)    {
                    nodeg = E->lmesh.nxs+i-1 + (E->lmesh.nys+k-2)*nox;
		    temp=E->flag_depth1[nodeg]-E->flag_depth[nodeg];
		    if(pos_age && fabs(temp)<0.5)
		        E->flag_depth2[nodeg]=E->flag_depth[nodeg]+temp/(newage2-newage1)*(age-newage1);
		    else
			E->flag_depth2[nodeg]=E->flag_depth[nodeg];
		}
            create_new_trench2(E,output_file6,E->tf_depth, 100.0);
            if(E->parallel.me==0) fprintf(stderr,"check after, age=%f\n",age);
            /*for (j=1; j<=E->sphere.caps_per_proc; j++) {
                for(ii=1;ii<=E->lmesh.noy;ii++)
                for(jj=1;jj<=E->lmesh.nox;jj++) {
                nodeg=E->lmesh.nxs-1+jj+(E->lmesh.nys+ii-2)*E->mesh.nox;
                mindist=E->flag_depth[nodeg];
                if(mindist<=0.99)
                    for(kk=1;kk<=E->lmesh.noz;kk++) {
                        node=kk+(jj-1)*E->lmesh.noz+(ii-1)*E->lmesh.nox*E->lmesh.noz;
                        rad=E->sx[j][3][node];
                        if(rad>0.975)
                            if(mindist<0.0)
                                E->T[j][node]=E->control.lith_age_mantle_temp;
                            else if(mindist<0.12) {
                                if(rad<=1-0.003*0.02)
                                    E->T[j][node]=E->control.lith_age_mantle_temp;
                            }
                            else if(mindist<0.6) {
                                if(rad<=1-(0.003+(mindist-0.12)/(0.6007-0.12))*0.02)
                                    E->T[j][node]=E->control.lith_age_mantle_temp;
                            }
                    }
                }
            }*/
        }
        else if(E->monitor.solution_cycles==0) {
            if(E->parallel.me==0) fprintf(stderr,"age=%f\n",age);
            create_new_trench2(E,output_file4,E->flag_depth, 1000.0);
            create_new_trench2(E,output_file5,E->flag_depth1, 1000.0);
	    for(k=1;k<=noy1;k++)
                for(i=1;i<=nox1;i++)    {
                    nodeg = E->lmesh.nxs+i-1 + (E->lmesh.nys+k-2)*nox;
		    temp=E->flag_depth1[nodeg]-E->flag_depth[nodeg];
		    if(pos_age && fabs(temp)<0.5)
		        E->flag_depth2[nodeg]=E->flag_depth[nodeg]+temp/(newage2-newage1)*(age-newage1);
		    else
			E->flag_depth2[nodeg]=E->flag_depth[nodeg];
		}
            create_new_trench2(E,output_file6,E->tf_depth, 100.0);
            if(E->parallel.me==0) fprintf(stderr,"check after, age2=%f\n",age);
	    fflush(stderr);
        }
	break;

      case 3:  /* read element materials */

        VIP1 = (float*) malloc ((emax+1)*sizeof(float));
        VIP2 = (float*) malloc ((emax+1)*sizeof(float));
        LL1 = (int*) malloc ((emax+1)*sizeof(int));
        LL2 = (int*) malloc ((emax+1)*sizeof(int));


        /* probably can be safely removed */
        for(m=1;m<=E->sphere.caps_per_proc;m++)
          for (el=1; el<=elx*ely*elz; el++)  {
            nodea = E->ien[m][el].node[2];
            llayer = layers(E,m,nodea);
            if (llayer)  { /* for layers:1-lithosphere,2-upper, 3-trans, and 4-lower mantle */
              E->mat[m][el] = llayer;
            }
          }
          for(i=1;i<=emax;i++)  {
               fscanf(fp1,"%d %d %f", &nn,&(LL1[i]),&(VIP1[i]));
               fscanf(fp2,"%d %d %f", &nn,&(LL2[i]),&(VIP2[i]));
          }

          fclose(fp1);
          fclose(fp2);

          for (m=1;m<=E->sphere.caps_per_proc;m++) {
            for (k=1;k<=ely;k++)   {
              for (i=1;i<=elx;i++)   {
                for (j=1;j<=elz;j++)  {
                  el = j + (i-1)*E->lmesh.elz + (k-1)*E->lmesh.elz*E->lmesh.elx;
                  elg = E->lmesh.ezs+j + (E->lmesh.exs+i-1)*E->mesh.elz + (E->lmesh.eys+k-1)*E->mesh.elz*E->mesh.elx;

                  E->VIP[m][el] = VIP1[elg]+(VIP2[elg]-VIP1[elg])/(newage2-newage1)*(age-newage1);
                  E->mat[m][el] = LL1[elg];

                }     /* end for j  */
              }     /*  end for i */
            }     /*  end for k  */
          }     /*  end for m  */

         free ((void *) VIP1);
         free ((void *) VIP2);
         free ((void *) LL1);
         free ((void *) LL2);

	break;

      } /* end switch */
    } /* end for m */

    fflush(E->fp);

    return;
}

/* Jiashun defined to create new trench */
void create_new_trench(struct All_variables *E, char *file_name, float *flag_depth, double maxdist) {
    int i,j,k,ii,jj,kk,node,nodeg,index,trenchid,numtrench,trench[200],*side,side_recorder,isfirst,islast,subzon;
    long seconds;
    double x,y,lat,lon,rad,**sz;
    double lon_rad,colat_rad,dist,mindist,dotprod,left,left1,left2,left3;
    char depthf[255],flagf[255],buf[500],cmd[200];
    FILE *Fsz,*Fdepth,*Fflag;
    void mindst();
    double cross_product();
    double dot_product();
    int Length;
    double PI=3.14159265359;
    double sin20=0.99;

    if(flag_depth==E->new_flag_depth)
	Length=2000;
    else if(flag_depth==E->tf_depth)
	Length=35000;
    else
	Length=12000;
    sz=(double **)malloc(2*sizeof(double *));
    sz[0]=(double *)malloc(Length*sizeof(double));
    sz[1]=(double *)malloc(Length*sizeof(double));
    side=(int *)malloc(Length*sizeof(int));
    Fsz=fopen(file_name,"r");
    i=0; /*record line number */
    numtrench=0; /*record trench number */
    while(fgets(buf,500,Fsz)!=NULL) {
        if(buf[0] == '>') {
            trench[numtrench]=i;
            numtrench++;
            if(buf[2] == 'L')
                side_recorder=1;
            else
                side_recorder=-1;
        }
        else {
            sscanf(buf,"%lf %lf",&lat,&lon);
            sz[0][i]=lon/180.0*PI;
            sz[1][i]=PI/2-lat/180.0*PI;
            side[i]=side_recorder;
            i++;
        }
        /*printf("%s\n",buf);*/
    }
    fclose(Fsz);

    for (j=1; j<=E->sphere.caps_per_proc; j++) {
    for(ii=1;ii<=E->lmesh.noy;ii++)
        for(jj=1;jj<=E->lmesh.nox;jj++) {
	    node = 1 + (jj-1)*E->lmesh.noz+(ii-1)*E->lmesh.noz*E->lmesh.nox;
            nodeg=E->lmesh.nxs-1+jj+(E->lmesh.nys+ii-2)*E->mesh.nox;
            if(i==0)
                flag_depth[nodeg]=1.0;
            else {
                mindist=1.0;
                lon_rad=E->sx[j][2][node];
                colat_rad=E->sx[j][1][node];
                for(trenchid=0; trenchid<numtrench-1;trenchid++) {
                    mindst(sz,trench,trenchid,lon_rad,colat_rad,&dist,&index);
                    if(dist>=maxdist)
                        continue;

                    /* this part of the code determines where the closest point is the first point or not */
                    isfirst=0;
                    islast=0;
                    if(trench[trenchid]==index)
                        isfirst=1;
                    if(trench[trenchid+1]==index+1)
                        islast=1;

                    if(isfirst==0 && islast==0) {
                        left1=cross_product(sz[0][index-1],sz[1][index-1],sz[0][index],sz[1][index],lon_rad,colat_rad);
                        left2=cross_product(sz[0][index],sz[1][index],sz[0][index+1],sz[1][index+1],lon_rad,colat_rad);
                        subzon=1;
                        if((left1>0 && left2>0 && side[index]>0) || (left1<0 && left2<0 && side[index]<0))
                            subzon=1;
                        else if((left1>0 && left2>0 && side[index]<0) || (left1<0 && left2<0 && side[index]>0))
                            subzon=-1;
                        else {
                            left3=cross_product(sz[0][index-1],sz[1][index-1],sz[0][index],sz[1][index],sz[0][index+1],sz[1][index+1]);
                            if(left3*left1>0) {
                                if(left2*side[index]>0)
                                    subzon=1;
                                else
                                    subzon=-1;
                            }
                            if(left3*left1<0) {
                                if(left1*side[index]>0)
                                    subzon=1;
                                else
                                    subzon=-1;
                            }
                        }
                        if(subzon==1) {
                            if(mindist>=0)
                                mindist=-dist/maxdist;
			    else if(mindist<0 && mindist<-dist/maxdist)
                                mindist=-dist/maxdist;
                        }
                        if(subzon==-1) {
                            if(mindist>dist/maxdist)
                                mindist=dist/maxdist;
                        }
                    }
                    else if(isfirst==1) {
                        dotprod=dot_product(sz[0][index],sz[1][index],sz[0][index+1],sz[1][index+1],lon_rad,colat_rad);
                        if(dotprod<-sin20 || dist>maxdist) {
                            continue;
                        }
                        if(dotprod<0 && dist>1.0*maxdist) {
                            continue;
                        }
                        else {
                            subzon=1;
                            left=cross_product(sz[0][index],sz[1][index],sz[0][index+1],sz[1][index+1],lon_rad,colat_rad);
                            if(left*side[index]<0)
                                subzon=-1;
                            if(subzon==1) {
                                if(mindist>=0)
                                    mindist=-dist/maxdist;
                                else if(mindist<0 && mindist<-dist/maxdist)
                                    mindist=-dist/maxdist;
                            }
                            if(subzon==-1) {
                                if(mindist>dist/maxdist)
                                    mindist=dist/maxdist;
                            }
                        }
                    }
                    else if(islast==1) {
                        dotprod=dot_product(sz[0][index],sz[1][index],sz[0][index-1],sz[1][index-1],lon_rad,colat_rad);
                        if(dotprod<-sin20 || dist>maxdist) {
                            continue;
                        }
                        if(dotprod<0 && dist>1.0*maxdist) {
                            continue;
                        }
                        else {
                            subzon=1;
			    left=cross_product(sz[0][index-1],sz[1][index-1],sz[0][index],sz[1][index],lon_rad,colat_rad);
                            if(left*side[index]<0)
                                subzon=-1;
                            if(subzon==1) {
                                if(mindist>=0)
                                    mindist=-dist/maxdist;
                                else if(mindist<0 && mindist<-dist/maxdist)
                                    mindist=-dist/maxdist;
                            }
                            if(subzon==-1) {
                                if(mindist>dist/maxdist)
                                    mindist=dist/maxdist;
                            }
                        }
                    }
                }
                flag_depth[nodeg]=mindist;
            }
        }
    }
    free(sz[0]);
    free(sz[1]);
    free(sz);
    free(side);
}


/* Jiashun defined to create new trench */
void create_new_trench2(struct All_variables *E, char *file_name, float *flag_depth, double maxdist) {
    int i,j,k,ii,jj,jjj,kk,node,nodeg,index,trenchid,numtrench,trench[300],*side,side_recorder,isfirst,islast,subzon;
    long seconds;
    double x,y,lat,lon,rad,**sz;
    double lon_rad,colat_rad,dist,mindist,dotprod,left,left1,left2,left3;
    char depthf[255],flagf[255],buf[500],cmd[200];
    FILE *Fsz,*Fdepth,*Fflag;
    void mindst2();
    double cross_product();
    double dot_product();
    int Length;
    double PI=3.14159265359;
    double sin20=0.99;

    if(flag_depth==E->new_flag_depth)
	Length=2000;
    else if(flag_depth==E->tf_depth)
	Length=35000;
    else
	Length=13500;
    sz=(double **)malloc(2*sizeof(double *));
    sz[0]=(double *)malloc(Length*sizeof(double));
    sz[1]=(double *)malloc(Length*sizeof(double));
    side=(int *)malloc(Length*sizeof(int));
    Fsz=fopen(file_name,"r");
    i=0; /*record line number */
    numtrench=0; /*record trench number */
    while(fgets(buf,500,Fsz)!=NULL) {
        if(buf[0] == '>') {
            trench[numtrench]=i;
	    /* remove trenches with a legnth of 3 */
            if(numtrench>0 && trench[numtrench]-trench[numtrench-1]<=3) {
                i=trench[numtrench-1];
                numtrench--;
            }
            numtrench++;
            if(buf[2] == 'L')
                side_recorder=1;
            else
                side_recorder=-1;
        }
        else {
            sscanf(buf,"%lf %lf",&lat,&lon);
            sz[0][i]=lon/180.0*PI;
            sz[1][i]=PI/2-lat/180.0*PI;
            side[i]=side_recorder;
            i++;
        }
        /*printf("%s\n",buf);*/
    }
    fclose(Fsz);

    for (j=1; j<=E->sphere.caps_per_proc; j++) {
    for(ii=1;ii<=E->lmesh.noy;ii++)
        for(jj=1;jj<=E->lmesh.nox;jj++) {
	    node = 1 + (jj-1)*E->lmesh.noz+(ii-1)*E->lmesh.noz*E->lmesh.nox;
            nodeg=E->lmesh.nxs-1+jj+(E->lmesh.nys+ii-2)*E->mesh.nox;
            if(i==0)
                flag_depth[nodeg]=1.0;
            else {
                mindist=1.0;
                lon_rad=E->sx[j][2][node];
                colat_rad=E->sx[j][1][node];
                
                mindst2(sz,i,lon_rad,colat_rad,&dist,&index);
                if(dist>=maxdist) {
		    flag_depth[nodeg]=1.0;
                    continue;
		}
		else if(flag_depth==E->tf_depth) {
                    flag_depth[nodeg]=dist/maxdist;
                    continue;
                }

                /* this part of the code determines where the closest point is the first point or not */
		isfirst=0;
		islast=0;
		for(jjj=0;jjj<numtrench;jjj++)
		    if(trench[jjj]==index)
			isfirst=1;
		    /* this requres the file has > in the last line */
		    else if(trench[jjj]==index+1)
			islast=1;
		    else if(trench[jjj]>index)
			break;

                    if(isfirst==0 && islast==0) {
                        left1=cross_product(sz[0][index-1],sz[1][index-1],sz[0][index],sz[1][index],lon_rad,colat_rad);
                        left2=cross_product(sz[0][index],sz[1][index],sz[0][index+1],sz[1][index+1],lon_rad,colat_rad);
                        subzon=1;
                        if((left1>0 && left2>0 && side[index]>0) || (left1<0 && left2<0 && side[index]<0))
                            subzon=1;
                        else if((left1>0 && left2>0 && side[index]<0) || (left1<0 && left2<0 && side[index]>0))
                            subzon=-1;
                        else {
                            left3=cross_product(sz[0][index-1],sz[1][index-1],sz[0][index],sz[1][index],sz[0][index+1],sz[1][index+1]);
                            if(left3*left1>0) {
                                if(left2*side[index]>0)
                                    subzon=1;
                                else
                                    subzon=-1;
                            }
                            if(left3*left1<0) {
                                if(left1*side[index]>0)
                                    subzon=1;
                                else
                                    subzon=-1;
                            }
                        }
                        if(subzon==1) {
                            if(mindist>=0)
                                mindist=-dist/maxdist;
			    else if(mindist<0 && mindist<-dist/maxdist)
                                mindist=-dist/maxdist;
                        }
                        if(subzon==-1) {
                            if(mindist>dist/maxdist)
                                mindist=dist/maxdist;
                        }
                    }
                    else if(isfirst==1) {
                        dotprod=dot_product(sz[0][index],sz[1][index],sz[0][index+1],sz[1][index+1],lon_rad,colat_rad);
                        if(dotprod<-sin20 || dist>maxdist) {
			    flag_depth[nodeg]=1.0;
                            continue;
                        }
                        if(dotprod<0 && dist>1.0*maxdist) {
			    flag_depth[nodeg]=1.0;
                            continue;
                        }
                        else {
                            subzon=1;
                            left=cross_product(sz[0][index],sz[1][index],sz[0][index+1],sz[1][index+1],lon_rad,colat_rad);
                            if(left*side[index]<0)
                                subzon=-1;
                            if(subzon==1) {
                                if(mindist>=0)
                                    mindist=-dist/maxdist;
                                else if(mindist<0 && mindist<-dist/maxdist)
                                    mindist=-dist/maxdist;
                            }
                            if(subzon==-1) {
                                if(mindist>dist/maxdist)
                                    mindist=dist/maxdist;
                            }
                        }
                    }
                    else if(islast==1) {
                        dotprod=dot_product(sz[0][index],sz[1][index],sz[0][index-1],sz[1][index-1],lon_rad,colat_rad);
                        if(dotprod<-sin20 || dist>maxdist) {
			    flag_depth[nodeg]=1.0;
                            continue;
                        }
                        if(dotprod<0 && dist>1.0*maxdist) {
			    flag_depth[nodeg]=1.0;
                            continue;
                        }
                        else {
                            subzon=1;
			    left=cross_product(sz[0][index-1],sz[1][index-1],sz[0][index],sz[1][index],lon_rad,colat_rad);
                            if(left*side[index]<0)
                                subzon=-1;
                            if(subzon==1) {
                                if(mindist>=0)
                                    mindist=-dist/maxdist;
                                else if(mindist<0 && mindist<-dist/maxdist)
                                    mindist=-dist/maxdist;
                            }
                            if(subzon==-1) {
                                if(mindist>dist/maxdist)
                                    mindist=dist/maxdist;
                            }
                        }
                    }
                
                flag_depth[nodeg]=mindist;
            }
        }
    }
    free(sz[0]);
    free(sz[1]);
    free(sz);
    free(side);
}

/* Jiashun defined to find out the minimum distance */
/* x is longitude, y is colat in radius */
void mindst(double **sz, int *trench, int num, double x, double y, double *min, int *index) {
    int i,tmp_index;
    double dst,tmp,R=6371.0;
    //double delx,dely,delxTmp,delyTmp; /* This filters out some points */

    dst=20000.0;
    tmp_index=0;
    for(i=trench[num];i<trench[num+1];i++) {
            tmp=R*acos(sin(y)*sin(sz[1][i])*cos(x-sz[0][i])+cos(y)*cos(sz[1][i]));
            if(tmp<dst) {
                dst=tmp;
                tmp_index=i;
            }
    }
    *min=dst;
    *index=tmp_index;
}

/* Jiashun defined to find out the minimum distance */
/* x is longitude, y is colat in radius */
void mindst2(double **sz, int num, double x, double y, double *min, int *index) {
    int i,tmp_index;
    double dst,tmp,R=6371.0;
    //double delx,dely,delxTmp,delyTmp; /* This filters out some points */

    dst=20000.0;
    tmp_index=0;
    //delx=10.0;
    //dely=5.0;
    for(i=0;i<num;i++) {
	//delxTmp=abs(x-sz[0][i]);
	//if(delxTmp>PI) delxTmp=2*PI-delxTmp;
	//delyTmp=abs(y-sz[1][i]);
	//if(delxTmp<delx || delyTmp<dely) {
            tmp=R*acos(sin(y)*sin(sz[1][i])*cos(x-sz[0][i])+cos(y)*cos(sz[1][i]));
	    //printf("i=%d,x=%lf,y=%lf,gx=%lf,gy=%lf,dist=%lf\n",i,x,y,sz[0][i],sz[1][i],tmp);
            if(tmp<dst) {
        	dst=tmp;
		tmp_index=i;
		//delx=delxTmp;
		//dely=delyTmp;
	    }
	//}
    }   
    *min=dst;
    *index=tmp_index;
}

/* Jiashun defined the cross produce */
/* (x, y) is test point, (x1,y1) is the point before the close point, (x2, y2) is the cloest point
   return > 0 if the test point is on the left side of the segment */
double cross_product(double x1, double y1, double x2, double y2, double x, double y) {
    double PI=3.14159265359;
    if(x-x1>PI) x=x-2*PI;
    if(x1-x>PI) x=x+2*PI;
    return((x2-x1)*(y2-y)-(x2-x)*(y2-y1));
}

/* Jiashun defined the dot product */
/* (x, y) is the test point, (x1, y1) is the closest point */
double dot_product(double x1, double y1, double x2, double y2, double x, double y) {
    double PI=3.14159265359;
    if(x-x1>PI) x=x-2*PI;
    if(x1-x>PI) x=x+2*PI;
    return(((x1-x2)*(x1-x)+(y1-y2)*(y1-y))/(sqrt(pow(x1-x2,2)+pow(y1-y2,2))*sqrt(pow(x1-x,2)+pow(y1-y,2))));
}
