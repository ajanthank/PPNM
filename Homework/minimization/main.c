#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <float.h>
#include <gsl/gsl_blas.h>
#include <assert.h>

double *E;
double *cross;
double *sigma;
double Ndat=30;

void num_grad(double fun(gsl_vector* x), gsl_vector* x, gsl_vector* gradient){
  int n=x->size;
  double funval=fun(x);

  for(int i=0; i<n;i++){
    double step;
    double xi=gsl_vector_get(x,i);

    if(fabs(xi)<DBL_EPSILON){
      step=DBL_EPSILON;
    }
    else{
      step=fabs(xi)*DBL_EPSILON;
    }

    gsl_vector_set(x,i,xi+step);
    gsl_vector_set(gradient, i, (fun(x)-funval)/step);
    gsl_vector_set(x,i,xi-step);
  }
}

void qnewton(double fun(gsl_vector*), gsl_vector* x, double tol){
  int n=x->size;
  double DELTA=sqrt(DBL_EPSILON);

  int nsteps=0, nrests=0, nscale=0;

  gsl_matrix* Ihesse=gsl_matrix_alloc(n,n);
  gsl_matrix_set_identity(Ihesse);

  gsl_matrix* I=gsl_matrix_alloc(n,n);
  gsl_matrix_set_identity(I);

  gsl_vector* grad=gsl_vector_alloc(n);
  gsl_vector* newtonstep=gsl_vector_alloc(n);
  gsl_vector* dx=gsl_vector_alloc(n);
  gsl_vector* dgrad=gsl_vector_alloc(n);
  gsl_vector* sol=gsl_vector_alloc(n);
  gsl_vector* sol2=gsl_vector_alloc(n);
  gsl_vector* B=gsl_vector_alloc(n);

  num_grad(fun,x,grad);
  double funval=fun(x);
  double dfunval;

  while(nsteps<1e4){
    nsteps++;

    gsl_blas_dgemv(CblasNoTrans, -1., Ihesse, grad,0,newtonstep);
    if(gsl_blas_dnrm2(newtonstep)<DELTA*gsl_blas_dnrm2(x)){
      fprintf(stderr,"qnewton: |DX|<DELTA*|x|\n");
      break;
    }
    if(gsl_blas_dnrm2(grad)<tol){
      fprintf(stderr,"qnewton: |grad|<tolorance\n");
      break;
    }

    double scale=1;

    while(1){
      gsl_vector_memcpy(dx,x);
      gsl_vector_add(dx, newtonstep);

      dfunval=fun(dx);

      double sTg;
      gsl_blas_ddot(newtonstep, grad, &sTg);

      if(dfunval<funval+0.01*sTg){
	nscale++;
	break;
      }
      if(scale<DELTA){
	nrests++;
	gsl_matrix_set_identity(Ihesse);
	break;
      }
      scale*=0.5;
      gsl_vector_scale(newtonstep,0.5);
    }
    num_grad(fun,dx,dgrad);

    gsl_vector_memcpy(sol,dgrad);
    gsl_blas_daxpy(-1,grad,sol);
    gsl_vector_memcpy(sol2,newtonstep);

    gsl_blas_dgemv(CblasNoTrans,-1,Ihesse,sol,1,sol2);

    gsl_matrix* sol2sol2T=gsl_matrix_calloc(n,n);
    gsl_blas_dsyr(CblasUpper, 1.0, sol2, sol2sol2T);
    double s2T2;
    gsl_blas_ddot(sol2,sol,&s2T2);
    if(fabs(s2T2)>1e-12){
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0/s2T2, sol2sol2T,I,1.0, Ihesse);
    }

    gsl_vector_memcpy(x,dx);
    gsl_vector_memcpy(grad,dgrad);
    funval=dfunval;
  }

  

  gsl_matrix_free(Ihesse);
  gsl_matrix_free(I);
  gsl_vector_free(grad);
  gsl_vector_free(newtonstep);
  gsl_vector_free(dx);
  gsl_vector_free(dgrad);
  gsl_vector_free(sol);
  gsl_vector_free(sol2);
  gsl_vector_free(B);

  fprintf(stderr,"qnewton: nsteps=%i nscale=%i nresets=%i fx=%.1e\n",nsteps,nscale,nrests,funval);
}

double Rosen(gsl_vector* x){
  double y=gsl_vector_get(x,0);
  double z=gsl_vector_get(x,1);

  double funval=pow(1-y,2)+100*pow(z-y*y,2);
  return funval;
}

double Himb(gsl_vector* x){
  double y=gsl_vector_get(x,0);
  double z=gsl_vector_get(x,0);
  
  double funval=pow(y*y+z-11,2)+pow(y+z*z-7,2);
  return funval;
}

double bw(double m, double width, double scale, double Energy){
  return scale/((Energy-m)*(Energy-m)+(width*width)/4);
}

double dif(gsl_vector* x){
  double mass=gsl_vector_get(x,0);
  double width=gsl_vector_get(x,1);
  double scale=gsl_vector_get(x,2);

  double tot=0;
  for(int i=0; i<Ndat;i++){
    tot+=pow(bw(mass,width,scale,E[i])-cross[i],2)/pow(sigma[i],2);
  }
  return 0;
}
  
int main(){

  printf("------------------ Exercise A ------------------\n");

  double xmin=1, ymin=1;
  double ixmin=xmin-1, iymin=ymin-1;

  int n=2;
  double tol=1e-5;

  gsl_vector* x=gsl_vector_alloc(n);
  gsl_vector_set(x,0,ixmin);
  gsl_vector_set(x,0,iymin);

  qnewton(Rosen,x,tol);

  printf("The Rosen.-function, with initial guess, %10g and %10g\n", ixmin, iymin);

  printf("Calculated minima at (x,y)=(%10g,%10g)\n",gsl_vector_get(x,0),gsl_vector_get(x,1));

  xmin=3;
  ymin=2;

  ixmin=xmin-0.2;
  iymin=ymin-0.2;

  gsl_vector_set(x,0,ixmin);
  gsl_vector_set(x,0,iymin);


  qnewton(Himb,x,tol);

  printf("The Himmelb.-function, with initial guess, %10g and %10g\n",ixmin, iymin);

  printf("Minima at (x,y)=(%10g,%10g)\n",gsl_vector_get(x,0),gsl_vector_get(x,1));

  printf("------------------- Exercise B ------------------\n");

  double mass=125.3;

  double E[]={101,103,105,107,109,111,113,115,117,119,121,123,125,127,129,131,133,135,137,139,141,143,145,147,149,151,153,155,157,159};
  double cross[]={-0.25,-0.30,-0.15,-1.71,0.81,0.65,-0.91,0.91,0.96,-2.52,-1.01,2.01,4.83,4.58,1.26,1.01,-1.26,0.45,0.15,-0.91,-0.81,
    -1.41,1.36,0.50,-0.45,1.61,-2.21,-1.86, 1.76,-0.50};

  double sigma[]={2.0,2.0,1.9,1.9,1.9,1.9,1.9,1.9,1.6,1.6,1.6,1.6,1.6,1.6,1.3,1.3,1.3,1.3,1.3,1.3,1.1,1.1,1.1,1.1,1.1,1.1,1.1,0.9,0.9,0.9};

  n=3;
  gsl_vector* xhiggs=gsl_vector_alloc(n);
  gsl_vector_set(xhiggs,0,mass+1);
  gsl_vector_set(xhiggs,1,2.5);
  gsl_vector_set(xhiggs,2,9);

  qnewton(dif,xhiggs,tol);

  printf("\n\n\n");
  printf("# Higgsfit\n");

  for(int i=0;i<Ndat;i++){
    double funval=bw(gsl_vector_get(xhiggs,0),gsl_vector_get(xhiggs,1), gsl_vector_get(xhiggs,2),E[i]);
    printf("%10g %10g %10g %10g\n",E[i],cross[i],sigma[i], funval);
  }

  gsl_vector_free(x);
  gsl_vector_free(xhiggs);

  return 0;
}
