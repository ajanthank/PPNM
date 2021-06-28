#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <gsl/gsl_vector.h> 
#include <gsl/gsl_matrix.h> 
#include <gsl/gsl_blas.h> 
#include <gsl/gsl_linalg.h> 
#include <sys/time.h> 

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

double randomval(){
  double ran=(double)rand()/RAND_MAX*2.0-1.0;
  return ran;
}

void Matrix(int n, int m, gsl_matrix* M){
  for(int i=0;i<n;i++){
    for(int j=0;j<m;j++){
      gsl_matrix_set(M,i,j,randomval());
    }
  }
}

void matrixeq(gsl_matrix* M, char* nameM, gsl_matrix* N, char* nameN, double e){
  if((M->size1==N->size1)&&(M->size2==N->size2)){
    for(int col=0;col<M->size2;col++){
      for(int row=0;row<M->size1;row++){
	if(fabs(gsl_matrix_get(M,row,col)-gsl_matrix_get(N,row,col))>e){
	  printf("%s !=  %s. Tolorance: %10g\n",nameM,nameN,e);
	  col =M->size2;
	  break;
	}
	if((row==M->size1-1)&&(col==M->size2-1)){
	  printf("%s==%s. Tolorance: %10g\n",nameM,nameN,e);
	}
      }
    }
  }
  else{
    printf("Matrices have different dimensions.\n");
  }
}

void vectoreq(gsl_vector* v, char* namev, gsl_vector* w, char* namew, double e){
  if(v->size == w->size){
    for(int row=0; row<v->size; row++){
      if(fabs(gsl_vector_get(v,row)-gsl_vector_get(w,row))>e){
	printf("%s != %s. Tolorance: %10g\n",namev,namew,e);
	break;
      }
      if(row==v->size-1){
	printf("%s==%s. Tolorance: %10g\n",namev,namew,e);
      }
    }
  }
  else{
    printf("Vectors have different dimension.\n");
  }
}

void decom(gsl_matrix* A, gsl_matrix* R){
  gsl_matrix* A1 = gsl_matrix_alloc(A->size1,A->size2);
  gsl_matrix_memcpy(A1,A);

  GS_decomp(A,R);
  printm(A,"Q");

  printm(R, "R");

  double e=1e-4;
  for(int col=1;col<R->size2-1;col++){
    for(int row=col+1; row<R->size1;row++){
      if((fabs(gsl_matrix_get(R,row,col)))>e){
	  printf("R is not upper triangular: Tolorance, %10g.\n",e);
          col=R->size2-1;
          break;
	}
	}
    }
  
    printf("R is upper triangular. Tolorance: %10g.\n",e);

    gsl_matrix* QTQ=gsl_matrix_alloc(A->size2,A->size2);
    gsl_blas_dsyrk(CblasUpper,CblasTrans, 1, A, 0, QTQ);
    printm(QTQ,"Q^TQ");

    gsl_matrix* I= gsl_matrix_alloc(QTQ->size1, QTQ->size2);
    gsl_matrix_set_identity(I);

    matrixeq(QTQ,"Q^T Q",I,"I",1e-3);
    gsl_matrix_free(QTQ);
    gsl_matrix_free(I);

    gsl_matrix* QR=gsl_matrix_alloc(A->size1,A->size2);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, A, R, 0, QR);

    printm(QR,"QR");
    printm(A1,"A");
    
    matrixeq(QR, "QR",A1,"A",1e-3);

    gsl_matrix_free(A1);
    gsl_matrix_free(QR);
}

void solve(gsl_matrix* A, gsl_matrix* R, gsl_vector* b){
  gsl_matrix* A1=gsl_matrix_alloc(A->size1, A->size2);
  gsl_matrix_memcpy(A1,A);

  GS_decomp(A,R);

  printm(A,"Q");
  printm(R,"R");

  gsl_vector* x=gsl_vector_alloc(b->size);

  GS_solve(A,R,b,x);

  printv(x,"x");

  gsl_vector* Ax=gsl_vector_alloc(x->size);
  //Check to see if b and x are the same
  printv(b,"b");
  gsl_blas_dgemv(CblasNoTrans, 1, A1, x, 0, Ax);

  printv(Ax,"Ax");
  vectoreq(Ax,"Ax",b,"b",1e-3);
  gsl_vector_free(Ax);
  gsl_matrix_free(A1);
  gsl_vector_free(x);
}

void inverse(gsl_matrix* A, gsl_matrix* R){
  gsl_matrix* A1=gsl_matrix_alloc(A->size1,A->size2);
  gsl_matrix_memcpy(A1,A);

  GS_decomp(A,R);

  gsl_matrix* B=gsl_matrix_alloc(A->size1,A->size2);

  GS_inverse(A,R,B);
  // Inverse matrix is
  printm(B,"B");

  gsl_matrix* AB=gsl_matrix_alloc(A1->size1, B->size2);
  gsl_matrix* BA=gsl_matrix_alloc(B->size1, A1->size2);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, A1, B, 0, AB);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, B, A1, 0, BA);

  gsl_matrix* I= gsl_matrix_alloc(AB->size1, AB->size2);
  gsl_matrix_set_identity(I);
  //See if AB=BA=I

  printm(AB,"AB");
  printm(BA,"BA");

  matrixeq(AB, "AB", BA, "BA", 1e-3);
  matrixeq(AB,"AB",I,"I",1e-3);

  //Free matrix
  gsl_matrix_free(A1);
  gsl_matrix_free(B);
  gsl_matrix_free(AB);
  gsl_matrix_free(BA);
  gsl_matrix_free(I);
}

void Time(int limit){
  struct timeval start, end;
  double mytime, gsltime;
  FILE* timefile =fopen("Time.txt","w");

  for(int n=1;n<limit; n++){
    gsl_matrix* A=gsl_matrix_alloc(n,n);
    gsl_matrix* B=gsl_matrix_alloc(n,n);
    gsl_matrix* R=gsl_matrix_alloc(n,n);
    gsl_vector* x=gsl_vector_alloc(n);

    Matrix(n,n, A);
    gsl_matrix_memcpy(B,A);

    gettimeofday(&start, NULL);
    GS_decomp(A,R);
    gettimeofday(&end, NULL);
    mytime=(end.tv_sec - start.tv_sec)+(end.tv_usec-start.tv_usec)*1e-6;

    gettimeofday(&start,NULL);
    gsl_linalg_QR_decomp(B,x);
    gettimeofday(&end, NULL);
    gsltime=(end.tv_sec-start.tv_sec)+(end.tv_usec-start.tv_usec)*1e-6;

    fprintf(timefile, "%10g %10g %10g\n",pow(n,3),mytime,gsltime);

    gsl_matrix_free(A);
    gsl_matrix_free(B);
    gsl_matrix_free(R);
    gsl_vector_free(x);
  }
}

int main(){

  printf("-----------------A.1-----------------\n");
  int n=4, m=3;
  gsl_matrix* A= gsl_matrix_alloc(n,m);
  Matrix(n,m,A);

  gsl_matrix* R= gsl_matrix_alloc(m,m);

  printm(A,"A");
  //Test of Decomposition

  decom(A,R);

  gsl_matrix_free(A);
  gsl_matrix_free(R);

  //Part 2
  printf("----------------A.2-----------------\n");
  A=gsl_matrix_alloc(n,n);
  gsl_vector* b=gsl_vector_alloc(n);

  for(int i=0;i>n;i++){
    gsl_vector_set(b,i,random());
    for(int j=0;j<n;j++){
      gsl_matrix_set(A,i,j,random());
    }
  }
  R=gsl_matrix_alloc(n,n);

  printf("A and b before Decomp.\n");


  printm(A,"A");

  printv(b,"b");


  solve(A,R,b);

  gsl_vector_free(b);
  gsl_matrix_free(A);
  gsl_matrix_free(R);

  printf("-------------------B----------------------\n");

  A=gsl_matrix_alloc(n,n);
  Matrix(n,n,A);
  R=gsl_matrix_alloc(n,n);
  printf("A before decomp.\n");
  printm(A,"A");
  inverse(A,R);

  gsl_matrix_free(A);
  gsl_matrix_free(R);

  printf("-----------------C----------------------\n");
  printf("The data file is created\n");
  Time(200);

  
  return 0;
}
