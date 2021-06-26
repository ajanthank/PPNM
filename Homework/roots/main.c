#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <assert.h>
#include <stdlib.h>
#include <gsl/gsl_linalg.h>

double E, R_Max, E_B, R_Max_B;

void rkstep12(void f(double t, gsl_vector* y, gsl_vector* dydt), double t, gsl_vector* yt, double h, gsl_vector* yh, gsl_vector* err){
  int n=yt->size;
  gsl_vector* k0=gsl_vector_alloc(n);
  gsl_vector* k1=gsl_vector_alloc(n);
  gsl_vector* yx=gsl_vector_alloc(n);

  //First order
  f(t,yt,k0);
  for(int i=0;i<n;i++){
    double yti=gsl_vector_get(yt,i);
    double k0i=gsl_vector_get(k0,i);
    double yxi=yti+h/2*k0i;
    gsl_vector_set(yx, i, yxi);
  }

  //Second order
  f(t+h/2,yx,k1);
  for(int i=0;i<n;i++){
    double yti=gsl_vector_get(yt,i);
    double k1i=gsl_vector_get(k1,i);
    double yhi=yti+h*k1i;
    gsl_vector_set(yh,i,yhi);
  }

  //Error
  for(int i=0;i<n;i++){
    double k0i=gsl_vector_get(k0,i);
    double k1i=gsl_vector_get(k1,i);
    double erri=(k0i-k1i)*h/2;
    gsl_vector_set(err,i,erri);
  }

  gsl_vector_free(k0);
  gsl_vector_free(k1);
  gsl_vector_free(yx);
}

int driver2(void f(double t, gsl_vector* y, gsl_vector* dydt),
       double a, gsl_vector* ya, double **Yal, int steps, double **Xal,double b, double h, double acc, double eps){
  int k=0;
  double dy,normy,tol;
  int n=ya->size;

  gsl_vector* err=gsl_vector_alloc(n);
  gsl_matrix_view Y=gsl_matrix_view_array(*Yal, steps, n);
  gsl_vector_view X=gsl_vector_view_array(*Xal,steps);
  gsl_matrix_set_row(&Y.matrix,0,ya);
  gsl_vector* yplaceholder=gsl_vector_alloc(n);

  double xpos=a;
  while(xpos<b){
    if(xpos+h>b){
      h=b-xpos;
    }
    gsl_vector_view yx=gsl_matrix_row(&Y.matrix,k);

    rkstep12(f,xpos,&yx.vector,h,yplaceholder,err);

    dy=gsl_blas_dnrm2(err);
    normy=gsl_blas_dnrm2(yplaceholder);
    tol=(normy*eps+acc)*sqrt(h/(b-a));

    if(dy<tol){
      k++;
      if(k>=steps){
    return -k;
    Yal=realloc(*Yal, sizeof(double)*(n*(k+1)));
    Xal=realloc(*Xal, sizeof(double)*(k+1));
    Y=gsl_matrix_view_array(*Yal, k+1,n);

    X=gsl_vector_view_array(*Xal,k+1);
      }
      xpos+=h;
      gsl_vector_set(&X.vector,k,xpos);
      gsl_matrix_set_row(&Y.matrix,k,yplaceholder);
    }
    if(dy>0){
      h*=pow(tol/dy,0.25)*0.95;
    }
    else{
      h*=2;
    }
  }
  gsl_vector_free(err);
  gsl_vector_free(yplaceholder);
  return k+1;
}

int driver(void f(double x, gsl_vector* y, gsl_vector* dydx) ,double a, gsl_vector* ya, double b, double h, double acc, double eps){
  int k = 0;
  double dy,normy,tol;
  int n = ya-> size;
  gsl_vector* err = gsl_vector_alloc(n);
  gsl_vector* yplaceholder = gsl_vector_alloc(n);
  double xpos = a;
  while (xpos < b ){
    if (xpos+h>b) h=b-xpos;
    rkstep12(f,xpos,ya,h, yplaceholder,err);
    dy = gsl_blas_dnrm2(err);
    normy = gsl_blas_dnrm2(yplaceholder);
    tol = (normy*eps+acc)*sqrt(h/(b-a));
    if(dy<tol){
      xpos+=h;
      gsl_vector_memcpy(ya,yplaceholder);
    }
    if(dy>0) h*=pow(tol/dy,0.25)*0.95; else h*=2;
  }
  return k;
}

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

void newton(void fun(gsl_vector* x, gsl_vector* fx), gsl_vector* x, double eps){
  int n=x->size, N=0;
  double step=sqrt(DBL_EPSILON);
  gsl_matrix* J=gsl_matrix_alloc(n,n);
  gsl_matrix* R=gsl_matrix_alloc(n,n);
  gsl_vector* dx=gsl_vector_alloc(n);
  gsl_vector* sol=gsl_vector_alloc(n);
  gsl_vector* sol2=gsl_vector_alloc(n);
  gsl_vector* fx=gsl_vector_alloc(n);
  gsl_vector* fx2=gsl_vector_alloc(n);
  gsl_vector* dfx=gsl_vector_alloc(n);

  fun(x,fx);

  while(gsl_blas_dnrm2(fx)>eps){
    N++;
    assert(N<1e5);

    for(int i=0;i<n;i++){
      gsl_vector_set(x,i,gsl_vector_get(x,i)+step);
      fun(x,fx2);
      for(int j=0;j<n;j++){

	double funval2_i=gsl_vector_get(fx2, j);
	double funval_i=gsl_vector_get(fx, j);
	double Deltafunval=funval2_i-funval_i;
	double xij=Deltafunval/step;
	gsl_matrix_set(J,j,i,xij);
      }
      gsl_vector_set(x, i, gsl_vector_get(x,i)-step);
    }
    gsl_vector_scale(fx,-1.0);
    
    GS_decomp(J,R);
    GS_solve(J,R,fx,sol);
    gsl_vector_scale(fx, -1.0);

    double scale=2;

    while((gsl_blas_dnrm2(dfx) >= (1-scale/2)*gsl_blas_dnrm2(fx)) && scale>=0.02){
      scale/=2;
      gsl_vector_memcpy(sol2,sol);
      gsl_vector_scale(sol2,scale);
      gsl_vector_memcpy(dx, sol2);
      gsl_vector_add(dx,x);

      fun(dx,dfx);
    }
    gsl_vector_memcpy(x,dx);
    gsl_vector_memcpy(fx,dfx);

    if(gsl_blas_dnrm2(sol)<step){
      break;
    }
  }
    
  gsl_matrix_free(J);
  gsl_matrix_free(R);
  gsl_vector_free(dx);
  gsl_vector_free(sol);
  gsl_vector_free(sol2);
  gsl_vector_free(fx);
  gsl_vector_free(fx2);
  gsl_vector_free(dfx);
}

void Rosen(gsl_vector* x, gsl_vector* fx){
  double y=gsl_vector_get(x,0);
  double z=gsl_vector_get(x,1);
  gsl_vector_set(fx,0,(-1)*2*(1-y)+(-2*y)*2*100*(z-y*y));
  gsl_vector_set(fx,1,2*100*(z-y*y));
}

void SE(double x, gsl_vector* y, gsl_vector* dydx){
  gsl_vector_set(dydx, 0,gsl_vector_get(y,1));
  double fx=gsl_vector_get(y,0);
  double ddy=-2.*(E+1./x)*fx;
  gsl_vector_set(dydx,1,ddy);
}

void SE_B(double x, gsl_vector* y, gsl_vector* dydx){
  gsl_vector_set(dydx, 0,gsl_vector_get(y,1));
  double fx=gsl_vector_get(y,0);
  double ddy=-2.*(E_B+1./x)*fx;
  gsl_vector_set(dydx,1,ddy);
}



void wavefun(gsl_vector* x, gsl_vector* fx){
  double a = 1e-3, abs = 1e-3, eps = 1e-3;
  int n=2;
  E=gsl_vector_get(x,0);
  gsl_vector* ya = gsl_vector_alloc(n);
  gsl_vector_set(ya,0,(a-a*a));
  gsl_vector_set(ya,1,(1-2.*a));

  driver(SE,a,ya,R_Max,0.1,abs,eps);
  gsl_vector_set(fx,0,gsl_vector_get(ya,0));
}

void wavefun_B(gsl_vector* x, gsl_vector* fx){
  double a = 1e-3, abs = 1e-3, eps = 1e-3;
  int n=2;
  E_B=gsl_vector_get(x,0);
  gsl_vector* ya = gsl_vector_alloc(n);
  gsl_vector_set(ya,0,(a-a*a));
  gsl_vector_set(ya,1,(1-2.*a));

  driver(SE_B,a,ya,R_Max_B,0.1,abs,eps);
  gsl_vector_set(fx,0,gsl_vector_get(ya,0));
  double F1=gsl_vector_get(ya,0);
  double k=sqrt(-2.*E_B);
  double F2=R_Max_B*exp(-k*R_Max_B);
  gsl_vector_set(fx,0,F1-F2);
}

int main(){
  printf("---------------- Exercise A -----------------\n");
  gsl_vector* x=gsl_vector_alloc(2);

  gsl_vector_set(x,0,5);
  gsl_vector_set(x,1,5);

  newton(Rosen,x,1e-3);
  printf("Extremum of Rosen:\nx=%10g y=%10g\n",gsl_vector_get(x,0),gsl_vector_get(x,1));

  printf("---------------- Exercise B ----------------\n");

  gsl_vector* y=gsl_vector_alloc(1);
  gsl_vector* y_B=gsl_vector_alloc(1);
  gsl_vector_set(y,0,-3);
  gsl_vector_set(y_B,0,-1);

  R_Max=10;
  R_Max_B=0.5;

  newton(wavefun,y,1e-3);
  newton(wavefun_B,y_B,1e-3);

  printf("Lowest Energy: \nE0= %10g\nE0_B= %10g\n",E,E_B);

  printf("\n\n\n");
  printf("# Data for\n#xi      Theoratical          Calculated\n");

  double abs=1e-3, eps=1e-3, a=1e-3, b=R_Max;
  gsl_vector* ya=gsl_vector_calloc(2);
  gsl_vector_set(ya,0,(a-a*a));
  gsl_vector_set(ya,1,(1-2.*a));
  int n=2,m=500;
  double* Yal=malloc(sizeof(double)*(m*n));
  double* Xal=malloc(sizeof(double)*m);
  driver2(SE,a,ya,&Yal,m,&Xal,b,(b-a)/10,abs,eps);
  gsl_matrix_view Y=gsl_matrix_view_array(Yal,m,n);
  gsl_vector_view X=gsl_vector_view_array(Xal,m);

  for(int i=0;i<m;i++){
    double xi=gsl_vector_get(&X.vector,i);
    if(xi!=0){
      double theory=xi*exp(-xi);
      double cal=gsl_matrix_get(&Y.matrix,i,0);
      printf("%10g %10g %10g\n",xi,theory,cal);
    }
  }

  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_vector_free(y_B);
  gsl_vector_free(ya);
  free(Yal);
  free(Xal);

  
  return 0;
}
