#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define PI 3.14159265359
#define res 0.1
#define asm_type 1
#define new_trench 2
#define tf_type 3

int main(int argc, char **argv) {
    int type=new_trench;
    double maxdist;
    char szf[255];
    void create_new_trench2();

    maxdist=500.0;
    sprintf(szf,"%s",argv[1]);
    create_new_trench2(szf,maxdist,type);

}

/* Jiashun defined to create new trench */
void create_new_trench2(char *file_name, double maxdist, int type) {
    int i,j,k,jjj,kk,node,nodeg,index,trenchid,numtrench,trench[300],*side,side_recorder,isfirst,islast,subzon;
    float ii,jj;
    long seconds;
    double x,y,lat,lon,rad,**sz;
    double lon_rad,colat_rad,dist,mindist,dotprod,left,left1,left2,left3;
    char depthf[255],flagf[255],buf[500],cmd[200], *left_polarity="sL";
    FILE *Fsz,*Fdepth,*Flag_depth;
    void mindst2();
    double cross_product();
    double dot_product();
    int Length;
    double sin20=0.5;

    if(type==asm_type)
	Length=15000;
    else if(type==tf_type)
	Length=35000;
    else
	Length=1500;
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
	    if(type==tf_type)
                printf("%s %d %s\n",file_name, i, buf);
            if(strstr(buf,left_polarity) != NULL)
                side_recorder=1;
            else
                side_recorder=-1;
        }
        else {
            sscanf(buf,"%lf %lf",&lon,&lat);
            sz[0][i]=lon/180.0*PI;
            sz[1][i]=PI/2-lat/180.0*PI;
            side[i]=side_recorder;
            i++;
        }
    }
    fclose(Fsz);

    sprintf(flagf,"asm_depth.xy");
    Flag_depth=fopen(flagf,"w");
    for(ii=0;ii<=360;ii=ii+res)
        for(jj=-90;jj<=90;jj=jj+res) {
	    /*printf("%f, %f\n",ii,jj); */
            if(i==0)
		fprintf(Flag_depth,"%f %f %f\n",ii,jj,1.0);
            else {
                mindist=1.0;
                lon_rad=ii/180.0*PI;
                colat_rad=(90.0-jj)/180.0*PI;
                
                mindst2(sz,i,lon_rad,colat_rad,&dist,&index);
                if(dist>=maxdist) {
		    fprintf(Flag_depth,"%f %f %f\n",ii,jj,1.0);
                    continue;
		}
		else if(type==tf_type) {
		    fprintf(Flag_depth,"%f %f %f\n",ii,jj,1.0);
		    continue;
		}

                /* this part of the code determines where the closest point is the first point or not */
		isfirst=0;
		islast=0;
		for(jjj=0;jjj<numtrench;jjj++)
		    if(trench[jjj]==index)
			isfirst=1;
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
		            fprintf(Flag_depth,"%f %f %f\n",ii,jj,1.0);
                            continue;
                        }
                        if(dotprod<0 && dist>1.0*maxdist) {
		            fprintf(Flag_depth,"%f %f %f\n",ii,jj,1.0);
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
		            fprintf(Flag_depth,"%f %f %f\n",ii,jj,1.0);
                            continue;
                        }
                        if(dotprod<0 && dist>1.0*maxdist) {
		            fprintf(Flag_depth,"%f %f %f\n",ii,jj,1.0);
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
                
		fprintf(Flag_depth,"%f %f %f\n",ii,jj,mindist);
            }
        }
    fclose(Flag_depth);
    free((void *)sz[0]);
    free((void *)sz[1]);
    free((void *)sz);
    free((void *)side);
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
    if(x-x1>PI) x=x-2*PI;
    if(x1-x>PI) x=x+2*PI;
    return((x2-x1)*(y2-y)-(x2-x)*(y2-y1));
}

/* Jiashun defined the dot product */
/* (x, y) is the test point, (x1, y1) is the closest point */
double dot_product(double x1, double y1, double x2, double y2, double x, double y) {
    if(x-x1>PI) x=x-2*PI;
    if(x1-x>PI) x=x+2*PI;
    return(((x1-x2)*(x1-x)+(y1-y2)*(y1-y))/(sqrt(pow(x1-x2,2)+pow(y1-y2,2))*sqrt(pow(x1-x,2)+pow(y1-y,2))));
}
