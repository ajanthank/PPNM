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
double fun(int i, double t){
    switch(i){
    case 0: return 1; break;
    case 1: return -t; break;
    default: return NAN;
    }
}

void printm(gsl_matrix* M, char* name){
  double m_ij;

  printf("%s=\n[",name);
  for(int i=0;i<M->size1; i++){
    for(int j=0;j<M->size2;j++){
      m_ij=gsl_matrix_get(M,i,j);
      printf("%10g,  ",m_ij);
    }
    if(i<M->size1){
      printf("\n");
    }
  }
  printf("]\n");
}

void printv(gsl_vector* v, char* name){
  double v_i;
  printf("%s=\n[",name);
  for(int i=0;i<v->size;i++){
    v_i=gsl_vector_get(v,i);
    printf("%10g",v_i);
    if(i<v->size){
      printf("\n");
    }
  }
  printf("]\n");
}

void lsfit(int m,  double fun(int i, double t), gsl_vector* x, gsl_matrix* data, gsl_matrix* Cov){
  int n= data->size1;

  gsl_matrix* A=gsl_matrix_alloc(n,m);
  gsl_vector* b=gsl_vector_alloc(n);
  gsl_matrix* R=gsl_matrix_alloc(m,m);
  gsl_matrix* B=gsl_matrix_alloc(m,m);
  gsl_matrix* invR=gsl_matrix_alloc(m,m);
  gsl_matrix* I=gsl_matrix_alloc(m,m);

  for(int i=0;i<n;i++){
    for(int j=0;j<m;j++){
      gsl_matrix_set(A,i,j,fun(j,gsl_matrix_get(data,i,0))/gsl_matrix_get(data,j,2));
    }
  }

  //printm(A);

  for(int i=0;i<m;i++){
    gsl_vector_set(b,i,gsl_matrix_get(data,i,1)/gsl_matrix_get(data,i,2));
  }

  gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,A,A,0,R);

  GS_decomp(R,B);
  GS_inverse(R,B,Cov);

  GS_decomp(A,R);
  gsl_blas_dgemv(CblasTrans,1,A,b,0,x);
  backsub(R,x);

  gsl_matrix_free(A);
  gsl_matrix_free(R);
  gsl_matrix_free(invR);
  gsl_matrix_free(B);
  gsl_vector_free(b);
  gsl_matrix_free(I);
}




int main(){

  int n=7;
  int m=4;

  gsl_matrix* A=gsl_matrix_alloc(n,m);
  gsl_matrix* R=gsl_matrix_alloc(m,m);
  gsl_matrix* B=gsl_matrix_alloc(m,m);
  unsigned int rand_seed=1;

  for (int i=0; i<n; i++){
    for (int j=0; j<m; j++){
      gsl_matrix_set(A,i,j,1.0*rand_r(&rand_seed)/RAND_MAX);
    }
  }

  double X[]={1,  2,  3, 4, 6, 9,   10,  13,  15};
  double Y[]={117,100,88,72,53,29.5,25.2,15.2,11.1};

  int N=sizeof(X)/sizeof(X[0]);

  gsl_matrix* data=gsl_matrix_alloc(N,3);

  for(int i=0;i<N;i++){
    gsl_matrix_set(data,i,0,X[i]);
    gsl_matrix_set(data,i,1,Y[i]);
    gsl_matrix_set(data,i,2,Y[i]/20);
  }
  printf("The data in its matrix form, with the last column being the uncertainty:\n");
  printm(data,"Data");
  
  printf("The matrix A is:\n");
  printm(A,"A");

  GS_decomp(A,R);
  printm(A,"Q");
  printm(R,"R");
  gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,A,A,0,B);

  printm(B,"Q^TQ");
  

  int Nf=2; //number of functions fitted to

  gsl_vector* c=gsl_vector_alloc(Nf);


  gsl_matrix* Cov=gsl_matrix_alloc(Nf,Nf);


  lsfit(Nf,&fun,c,data,Cov);
  printm(Cov,"Covariance matrix");

  printv(c,"c");

  printf("∆c=\n");

  for (int i=0; i<Nf; i++){
    printf("%10g   ",sqrt(gsl_matrix_get(Cov,i,i)));
  }

  double ht= -1*log(1./2)/gsl_vector_get(c,1);
  double deltaht=sqrt(gsl_matrix_get(Cov,1,1));
  

  printf("The calculated halftime is: %10g +/- %10g\n",ht,deltaht);
  printf("Online: 3.63 days[https://en.wikipedia.org/wiki/Isotopes_of_radium]. So they do not agree.\n");

  FILE* fit=fopen("fit.txt","w");

  for (int i=1; j<=15; j++){
    fprintf(fit_stream,"%10i %10g\n",j,exp(gsl_vector_get(c,0)-gsl_vector_get(c,1)*j));
  }
  fprintf(fit_stream,"\n\n\n\n\n\n\n\n");
  for (int j=1; j<=15; j++){
    fprintf(fit_stream,"%10i %10g\n",j,exp(gsl_vector_get(c,0)+sqrt(gsl_matrix_get(Covy,0,0))-(gsl_vector_get(c,1)+sqrt(gsl_matrix_get(Covy,1,1)))*j));
  }
  fprintf(fit_stream,"\n\n\n\n\n\n\n\n");

  for (int j=1; j<=15; j++){
    fprintf(fit_stream,"%10i %10g\n",j,exp(gsl_vector_get(c,0)-sqrt(gsl_matrix_get(Covy,0,0))-(gsl_vector_get(c,1)-sqrt(gsl_matrix_get(Covy,1,1)))*j));
  }
  printf("\n\n\n\n\n\n\n\n");*/


  gsl_matrix_free(Cov);
  gsl_vector_free(c);
  gsl_matrix_free(B);
  gsl_matrix_free(A);
  gsl_matrix_free(R);
  gsl_matrix_free(data);

  return 0;
}
