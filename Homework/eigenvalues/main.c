#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <time.h>

void timesJ(gsl_matrix* A, int p, int q, double theta){
  double c=cos(theta),s=sin(theta);
  for(int i=0;i<A->size1;i++){
    double new_aip=c*gsl_matrix_get(A,i,p)-s*gsl_matrix_get(A,i,q);
    double new_aiq=s*gsl_matrix_get(A,i,p)+c*gsl_matrix_get(A,i,q);
    gsl_matrix_set(A,i,p,new_aip);
    gsl_matrix_set(A,i,q,new_aiq);
  }
}

void Jtimes(gsl_matrix* A, int p, int q, double theta){
  double c=cos(theta),s=sin(theta);
  for(int j=0;j<A->size2;j++){
    double new_apj= c*gsl_matrix_get(A,p,j)+s*gsl_matrix_get(A,q,j);
    double new_aqj=-s*gsl_matrix_get(A,p,j)+c*gsl_matrix_get(A,q,j);
    gsl_matrix_set(A,p,j,new_apj);
    gsl_matrix_set(A,q,j,new_aqj);
  }
}

void JDiag(gsl_matrix* A, gsl_matrix* V){
  gsl_matrix_set_identity(V);

  int n=A->size1;
  int changed;
  do{
    changed=0;
    for(int p=0;p<n-1;p++){
      for(int q=p+1;q<n;q++){
	double apq=gsl_matrix_get(A,p,q);
	double app=gsl_matrix_get(A,p,p);
	double aqq=gsl_matrix_get(A,q,q);
	double theta=0.5*atan2(2*apq,aqq-app);
	double c=cos(theta),s=sin(theta);
	double new_app=c*c*app-2*s*c*apq+s*s*aqq;
	double new_aqq=s*s*app+2*s*c*apq+c*c*aqq;
	if(new_app!=app || new_aqq!=aqq){ // do rotation
	    changed=1;
	    timesJ(A,p,q, theta);
	    Jtimes(A,p,q,-theta); // A←J^T*A*J
	    timesJ(V,p,q, theta); // V←V*J
	}
      }
    }
  }while(changed!=0);
}

void printm(gsl_matrix* M, char* name){
  double m_ij;

  printf("%s=\n[",name);
  for(int i=0;i<M->size1; i++){
    for(int j=0;j<M->size2;j++){
      m_ij=gsl_matrix_get(M,i,j);
      if(m_ij<1e-10){
	m_ij=0.0;
      }
      printf("%10g,  ",m_ij);
    }
    if(i==(M->size1)-1){
      printf("]\n");
    }
    if(i<M->size1){
      printf("\n");
    }
  }
}

double randomval(unsigned int* seed){
  double ran=(double)rand_r(seed)/(double)RAND_MAX;
  return ran;
}

void symmatrix(gsl_matrix* M, unsigned int* seed){
  for(int i=0;i<M->size1; i++){
    gsl_matrix_set(M, i, i, randomval(seed));
    for(int j=i + 1;j<M->size2; j++){
      double MIJ=randomval(seed);
      gsl_matrix_set(M, i, j, MIJ);
      gsl_matrix_set(M, j, i, MIJ);
    }
  }
}

int main(){

  printf("-------------------------Exercise A --------------------\n");

  int n=5;
  unsigned int seed=time(NULL);
  
  gsl_matrix* A=gsl_matrix_alloc(n,n);
  gsl_matrix* V=gsl_matrix_alloc(n,n);
  gsl_matrix* lambdaM=gsl_matrix_alloc(n,n);

  symmatrix(A,&seed);
  gsl_matrix_memcpy(lambdaM,A);

  printm(A,"A");
  JDiag(lambdaM, V);

  gsl_matrix* test=gsl_matrix_alloc(n,n);
  gsl_matrix* testI=gsl_matrix_alloc(n,n);
  gsl_matrix* testdia=gsl_matrix_alloc(n,n);
  gsl_matrix* testM=gsl_matrix_alloc(n,n);

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A, V, 0.0, test);
  gsl_blas_dgemm(CblasTrans,   CblasNoTrans, 1.0, V, test, 0.0, testdia);
  gsl_blas_dgemm(CblasTrans,   CblasNoTrans, 1.0, V, V, 0.0, testI);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans,   1.0, lambdaM, V, 0.0, test);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, V, test, 0.0, testM);
    
  printm(V,"Eigenvector matrix, V");
  printm(testI,"V^T*V=I");
  printm(lambdaM,"Eigenvalues");
  printm(testdia,"V^T*A*V=D");
  printm(A,"A");
  printm(testM,"V*D*V^T=A");

  printf("-------------------------Exercise B --------------------\n");

  FILE* output=fopen("Energy.txt", "w");

  int N=30;

  double s=1.0/(N+1);

  gsl_matrix* H=gsl_matrix_alloc(N,N);
  for(int i=0;i<N-1;i++){
    gsl_matrix_set(H,i,i,-2.0);
    gsl_matrix_set(H,i,i+1,1.0);
    gsl_matrix_set(H,i+1,i,1.0);
  }
  gsl_matrix_set(H,N-1,N-1,-2.0);
  gsl_matrix_scale(H,-1/s/s);

  gsl_matrix* ketE=gsl_matrix_alloc(N,N);
  JDiag(H,ketE);

  //See if it is correct.

  fprintf(output,"# k    calculated     exact\n");

  for(int k=0;k<N/3;k++){
    double exact=M_PI*M_PI*(k+1)*(k+1);
    double calculated=gsl_matrix_get(H,k,k);
    fprintf(output,"%i \t %10g \t %10g\n",k,calculated,exact);
  }
  FILE* output2=fopen("Eigenfunc.txt","w");
  fprintf(output2, "# lowest eigenfunctions\n");
  for(int k=0;k<4;k++){
    fprintf(output2,"\n \n");
    fprintf(output2,"%10d \t %10d\n",0,0);
    for(int i=0;i<N;i++){
      fprintf(output2,"%10g \t %10g\n",(i+1.0)/(N+1), gsl_matrix_get(ketE,i,k));
    }
    fprintf(output2,"%10d \t %10d\n",1,0);
  }

  gsl_matrix_free(A);
  gsl_matrix_free(V);
  gsl_matrix_free(lambdaM);
  gsl_matrix_free(test);
  gsl_matrix_free(testI);
  gsl_matrix_free(testdia);
  gsl_matrix_free(testM);
  gsl_matrix_free(H);
  gsl_matrix_free(ketE);
}
