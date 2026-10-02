/* Opt-in consistent particle projection. Absolute temperatures solve a bounded
 * least-squares problem; signed increments use the unbounded linear operator.
 * This is not a physical-particle-mass conservation scheme. */
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>
#include "global_defs.h"
#include "pices.h"
struct PICMatrix {struct All_variables *E;int nn,np;int *nodes;double *weights,*copies,*rows,*diag;};
static double *alloc(struct All_variables *E,size_t n) {double *x=calloc(n,sizeof(double));if(!x)pices_fail(E,"projection allocation failed");return x;}
static double global(struct All_variables *E,double x,MPI_Op op) {double y;MPI_Allreduce(&x,&y,1,MPI_DOUBLE,op,E->parallel.world);return y;}
static void share(struct All_variables *E,double *x) {double *a[NCS];a[1]=x;E->exchange_node_d(E,a,E->mesh.levmax);}
static int fixed_node(struct All_variables *E,int n) {return (E->node[1][n]&(TBX|TBY|TBZ))!=0;}
static void apply(struct PICMatrix *A,const double *x,double *y) {
 int p,a;memset(y,0,(A->nn+1)*sizeof(double));
 for(p=0;p<A->np;p++) {double q=0;for(a=0;a<8;a++)q+=A->weights[p*8+a]*x[A->nodes[p*8+a]];
  for(a=0;a<8;a++)y[A->nodes[p*8+a]]+=A->weights[p*8+a]*q;}
 share(A->E,y);
}
static double dot(struct PICMatrix *A,const double *x,const double *y) {
 int n;double s=0;for(n=1;n<=A->nn;n++)s+=x[n]*y[n]/A->copies[n];return global(A->E,s,MPI_SUM);
}
void pices_consistent_project(struct All_variables *E,const double *values,double *out,int boundary) {
 struct PICMatrix A;int n,p,a,cg=0,pg=0,empty=0,*covered;double lo=DBL_MAX,hi=-DBL_MAX,scale=0,tol,error=0,rz=0;
 double *b,*r,*z,*direction,*work;double active=0;
 A.E=E;A.nn=E->lmesh.nno;A.np=E->trace.ntracers[1];
 A.nodes=malloc((size_t)A.np*8*sizeof(int));if(!A.nodes)pices_fail(E,"projection index allocation failed");
 A.weights=alloc(E,(size_t)A.np*8);A.copies=alloc(E,A.nn+1);A.rows=alloc(E,A.nn+1);A.diag=alloc(E,A.nn+1);
 b=alloc(E,A.nn+1);r=alloc(E,A.nn+1);z=alloc(E,A.nn+1);direction=alloc(E,A.nn+1);work=alloc(E,A.nn+1);
 covered=calloc(E->lmesh.nel+1,sizeof(int));if(!covered)pices_fail(E,"projection coverage allocation failed");
 for(n=1;n<=A.nn;n++)A.copies[n]=1;share(E,A.copies);
 for(p=0;p<A.np;p++) {
  int nodes[9],e=E->trace.ielement[1][p+1];double w[9],v=values[p+1];
  if(!isfinite(v) || (boundary && E->data.Ttop+E->data.ref_temperature*v<0))pices_fail(E,"invalid projection particle value");
  if(e<1 || e>E->lmesh.nel)pices_fail(E,"invalid projection host element");covered[e]++;
  tracer_temperature_weights(E,1,p+1,nodes,w);lo=fmin(lo,v);hi=fmax(hi,v);scale=fmax(scale,fabs(v));
  for(a=0;a<8;a++) {n=nodes[a+1];A.nodes[p*8+a]=n;A.weights[p*8+a]=w[a+1];
   b[n]+=w[a+1]*v;A.rows[n]+=w[a+1];A.diag[n]+=w[a+1]*w[a+1];}
 }
 for(n=1;n<=E->lmesh.nel;n++)if(!covered[n])empty++;
 share(E,b);share(E,A.rows);share(E,A.diag);
 for(n=1;n<=A.nn;n++) {
  if(!(A.diag[n]>0) || !(A.rows[n]>0) || !isfinite(A.diag[n]))pices_fail(E,"consistent projection requires covered nodes");
  out[n]=boundary?E->T[1][n]:0;
  if(boundary) {if(!isfinite(out[n]) || E->data.Ttop+E->data.ref_temperature*out[n]<0)pices_fail(E,"invalid previous grid temperature");lo=fmin(lo,out[n]);hi=fmax(hi,out[n]);scale=fmax(scale,fabs(out[n]));}
 }
 lo=global(E,lo,MPI_MIN);hi=global(E,hi,MPI_MAX);scale=global(E,scale,MPI_MAX);tol=fmax(1e-14,1e-12*scale);
 apply(&A,out,work);
 for(n=1;n<=A.nn;n++) {r[n]=(boundary&&fixed_node(E,n))?0:b[n]-work[n];z[n]=r[n]/A.diag[n];direction[n]=z[n];error=fmax(error,fabs(r[n])/A.rows[n]);}
 error=global(E,error,MPI_MAX);rz=dot(&A,r,z);
 while(error>tol && cg<1000) {
  double pap,alpha,next,beta;apply(&A,direction,work);pap=dot(&A,direction,work);
  if(!(pap>0) || !isfinite(pap) || !isfinite(rz))pices_fail(E,"projection CG lost positive curvature");alpha=rz/pap;
  for(n=1;n<=A.nn;n++)if(!boundary || !fixed_node(E,n))out[n]+=alpha*direction[n];
  /* Recompute true residual periodically to avoid accepting recursive drift. */
  cg++;if(cg%25==0) {apply(&A,out,work);for(n=1;n<=A.nn;n++)r[n]=(boundary&&fixed_node(E,n))?0:b[n]-work[n];}
  else for(n=1;n<=A.nn;n++)r[n]=(boundary&&fixed_node(E,n))?0:r[n]-alpha*work[n];
  error=0;for(n=1;n<=A.nn;n++){z[n]=r[n]/A.diag[n];error=fmax(error,fabs(r[n])/A.rows[n]);}
  error=global(E,error,MPI_MAX);next=dot(&A,r,z);beta=next/rz;rz=next;
  for(n=1;n<=A.nn;n++)direction[n]=z[n]+beta*direction[n];
 }
 if(error>tol || !isfinite(error))pices_fail(E,"consistent projection CG did not converge");
 /* Row sums bound the scaled Hessian's largest eigenvalue by one. Projected
  * Jacobi is therefore a valid convex-QP descent iteration, not post-hoc clipping. */
 if(boundary)for(n=1;n<=A.nn;n++)if(!fixed_node(E,n))out[n]=fmax(lo,fmin(hi,out[n]));
 for(;;) {
  apply(&A,out,work);error=0;
  for(n=1;n<=A.nn;n++) {
   double next=out[n]+(b[n]-work[n])/A.rows[n];
   if(boundary)next=fixed_node(E,n)?E->T[1][n]:fmax(lo,fmin(hi,next));
   z[n]=next;error=fmax(error,fabs(next-out[n]));
  }
  error=global(E,error,MPI_MAX);
  if(!isfinite(error))pices_fail(E,"nonfinite projection KKT residual");
  if(error<=tol)break;
  if(!boundary || pg>=5000)pices_fail(E,"projection true residual/KKT check failed");
  for(n=1;n<=A.nn;n++)out[n]=z[n];pg++;
 }
 for(n=1;n<=A.nn;n++) {
  if(!isfinite(out[n]))pices_fail(E,"nonfinite projected field");
  if(boundary && !fixed_node(E,n) && (out[n]<=lo+tol || out[n]>=hi-tol))active+=1/A.copies[n];
 }
 active=global(E,active,MPI_SUM);
 if(E->parallel.me==0)fprintf(E->fp,"PICES_PROJECTION step=%d method=bounded_consistent_v1 kind=%s cg=%d qp=%d residual=%.17g tolerance=%.17g lower=%.17g upper=%.17g active=%.0f\n",E->monitor.solution_cycles,boundary?"absolute":"signed",cg,pg,error,tol,lo,hi,active);
 if(boundary)fprintf(E->fp,"PICES_COVERAGE step=%d empty_elements=%d zero_boundary_nodes=0\n",E->monitor.solution_cycles,empty);
 free(A.nodes);free(A.weights);free(A.copies);free(A.rows);free(A.diag);free(b);free(r);free(z);free(direction);free(work);free(covered);
}
