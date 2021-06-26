#include <stdio.h>
#include <stdlib.h> 
#include <math.h> 
#include <gsl/gsl_vector.h> 
#include <gsl/gsl_matrix.h> 
#include <gsl/gsl_blas.h> 
#include <gsl/gsl_linalg.h>
#include <assert.h>


void GS_decomp(gsl_matrix* A, gsl_matrix* R){
  int m=A->size2;
  
  double DotP;

  for(int i=0; i<m; i++){
    gsl_vector_view col=gsl_matrix_column(A,i);
    gsl_matrix_set(R, i, i, gsl_blas_dnrm2(&col.vector));
    gsl_vector_scale(&col.vector, 1/gsl_matrix_get(R, i, i));

    for(int j=i+1; j<m; j++){
      gsl_vector_view col2=gsl_matrix_column(A,j);
      gsl_blas_ddot(&col.vector, &col2.vector, &DotP);
      gsl_matrix_set(R,i,j,DotP);
      gsl_blas_daxpy(-gsl_matrix_get(R,i,j),&col.vector, &col2.vector);
    }
  }
}

void backsub(gsl_matrix* M, gsl_vector* v){
  for(int i=v->size-1;i>=0;i--){
    double Sum=gsl_vector_get(v,i);
    for(int j=i+1; j<v->size; j++){
      Sum-=gsl_matrix_get(M,i,j)*gsl_vector_get(v,j);
    }
    gsl_vector_set(v,i,Sum/gsl_matrix_get(M,i,i));
  }
}

void GS_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* v, gsl_vector* w){
  gsl_blas_dgemv(CblasTrans, 1, Q, v, 0, w);
  backsub(R,w);
}

void GS_inverse(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* B){
  gsl_vector* basevec=gsl_vector_alloc(Q->size2);
  for(int i=0; i<Q-> size2; i++){
    gsl_vector_set_basis(basevec,i);
    gsl_vector_view col=gsl_matrix_column(B,i);
    GS_solve(Q,R,basevec, &col.vector);
  }
  gsl_vector_free(basevec);
}

void lsfit(int m,  double f(int i, double x), gsl_vector* x, gsl_vector* y, gsl_vector* dy, gsl_vector* c, gsl_matrix* S){
  int n= x->size;

  gsl_matrix* A=gsl_matrix_alloc(n,m);
  gsl_vector* b=gsl_vector_alloc(n);
  gsl_matrix* R=gsl_matrix_alloc(m,m);
  gsl_matrix* invR=gsl_matrix_alloc(m,m);
  gsl_matrix* I=gsl_matrix_alloc(m,m);

  for(int i=0;i<n;i++){
    double xi=gsl_vector_get(x,i);
    double yi=gsl_vector_get(y,i);
    double dyi=gsl_vector_get(dy,i);
    assert(dyi>0);
    gsl_vector_set(b,i,yi/dyi);
    for(int j=0;j<m;j++){
      gsl_matrix_set(A,i,j,f(j,xi)/dyi);
    }
  }

  GS_decomp(A,R);
  GS_solve(A,R,b,c);

  gsl_matrix_set_identity(I);
  GS_inverse(I,R,invR);
  gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,invR,invR,0,S);

  gsl_matrix_free(A);
  gsl_matrix_free(R);
  gsl_matrix_free(invR);
  gsl_matrix_free(I);
  gsl_vector_free(b);
}

double fun(int i, double t){
    switch(i){
    case 0: return 1; break;
    case 1: return -t; break;
    default: return NAN;
    }
}

double fit(double x,int m,gsl_vector* c){
  double s=0;
  for(int i=0;i<m;i++){
    s+=gsl_vector_get(c,i)*fun(i,x);
  }
  return s;
}

double fit_plus(int i, double x, gsl_vector* dc,gsl_vector* c,int m){
  return fit(x,m,c)+gsl_vector_get(dc,i)*fun(i,x);
}
double fit_minus(int i, double x, gsl_vector* dc,gsl_vector* c, int m){
  return fit(x,m,c)-gsl_vector_get(dc,i)*fun(i,x);
}

int main(){
  double x[]={1,2,3,4,6,9,10,13,15};
  double y[]={117,100,88,72,53,29.5,25.2,15.2,11.1};
  int n=sizeof(x)/sizeof(x[0]);
  double dy[n];
  for(int i=0; i<n;i++){
    dy[i]=0.05*y[i];
  }
  for(int i=0;i<n;i++){
    dy[i]/=y[i];
    y[i]=log(y[i]);
    }
  gsl_vector* vx=gsl_vector_alloc(n);
  gsl_vector* vy=gsl_vector_alloc(n);
  gsl_vector* vdy=gsl_vector_alloc(n);

  for(int i=0;i<n;i++){
    gsl_vector_set(vx,i,x[i]);
    gsl_vector_set(vy,i,y[i]);
    gsl_vector_set(vdy,i,dy[i]);
  }

  int m=2;
  

  gsl_vector* c=gsl_vector_alloc(m);
  gsl_matrix* S=gsl_matrix_alloc(m,m);

  lsfit(m,fun,vx,vy,vdy,c,S);

  gsl_vector* dc=gsl_vector_alloc(m);

  for(int j=0;j<m;j++){
    double sjj=gsl_matrix_get(S,j,j);
    gsl_vector_set(dc,j,sqrt(sjj));
  }
  
  double cl=gsl_vector_get(c,1);
  double dcl=gsl_vector_get(dc,1);

  double T=1/cl*log(2.0);
  double dT=dcl/cl/cl;

  printf("The half-life of ThX is found to be, %g +/- %g days\n",T,dT);
  printf("The modern value is 3.63 days[Wikipedia], which seems to disagree with the values from the fit. This could be due to the contamination of other isotopes of radium.\n");

  FILE* fitdat=fopen("Fitdat.txt","w");

  fprintf(fitdat,"# time log(A) ∆(log(A))\n");
  for(int i=0;i<n;i++){
    fprintf(fitdat,"%g %g %g\n",x[i],y[i],dy[i]);
  }
  fprintf(fitdat,"\n\n");

  for(int i=0;i<m;i++){
    fprintf(fitdat,"# Time Fit Fit_plus Fit_minus; k=%i\n",i);
    for(double z=x[0],dz=(x[n-1]-x[0])/64;z<=x[n-1];z+=dz){
      fprintf(fitdat,"%g %g %g %g\n",z,fit(z,m,c),fit_plus(i,z,dc,c,m),fit_minus(i,z,dc,c,m));
    }
    fprintf(fitdat,"\n\n");
  }

  gsl_vector_free(vx);
  gsl_vector_free(vy);
  gsl_vector_free(vdy);
  gsl_vector_free(c);
  gsl_vector_free(dc);
  gsl_matrix_free(S);

  return 0;
}
