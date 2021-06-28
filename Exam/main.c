#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <assert.h>

void rkstep12(void f(double t, gsl_vector* y, gsl_vector* dydt), double t, gsl_vector* yt, double h, gsl_vector* yh, gsl_vector* err){
  int n=yt->size;
  gsl_vector* k0=gsl_vector_alloc(n);
  gsl_vector* k1=gsl_vector_alloc(n);
  gsl_vector* yx=gsl_vector_alloc(n);

  f(t,yt,k0);
  for(int i=0;i<n;i++){
    double yti=gsl_vector_get(yt,i);
    double k0i=gsl_vector_get(k0,i);
    double yxi=yti+h/2*k0i;
    gsl_vector_set(yx, i, yxi);
  }

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

void abstep2(void f(double t, gsl_vector* y, gsl_vector* dydt), double t, gsl_vector* yt, double h, gsl_vector* yh, gsl_vector* err,
	     double t1){
  int n=yt->size;
  gsl_vector* yx=gsl_vector_alloc(n);
  gsl_vector* k=gsl_vector_alloc(n);

  f(t,yt,k);
  for(int i=1;i<n;i++){
    double yti=gsl_vector_get(yt,i);
    double yti1=gsl_vector_get(yt,i-1);
    double ki=gsl_vector_get(k,i);
    double ci=(yti1-yti+ki*(t-t1))/(pow(t-t1,2));
    double yxi=yti+ki*h+ci*pow(h,2);
    double erri=ci*h*h;
    gsl_vector_set(yx,i,yxi);
    gsl_vector_set(err,i,erri);
  }
  gsl_vector_free(yx);
}

int abdriver(void f(double t, gsl_vector* y, gsl_vector* dydt),
	     double a, gsl_vector* ya, double **Yal, int steps, double **Xal, double b, double h, double acc, double eps){
  int k=0;
  double dy,normy,tol;
  int n=ya->size;

  gsl_vector* err=gsl_vector_alloc(n);
  gsl_matrix_view Y=gsl_matrix_view_array(*Yal, steps, n);
  gsl_vector_view X=gsl_vector_view_array(*Xal,steps);
  gsl_matrix_set_row(&Y.matrix,0,ya);
  gsl_vector* yplaceholder=gsl_vector_alloc(n);

  double xpos=a,xpos1=0;
  int abrk=0;
  if(abrk==0){
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
	xpos1=xpos;
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
  }
  else{
    while(xpos<b){
      if(xpos+h>b){
	h=b-xpos;
      }
      gsl_vector_view yx=gsl_matrix_row(&Y.matrix,k);
      abstep2(f,xpos,&yx.vector,h,yplaceholder,err,xpos1);
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
	xpos1=xpos;
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
  }
  abrk++;
  gsl_vector_free(err);
  gsl_vector_free(yplaceholder);
  return k+1;
}


void harm(double t, gsl_vector* y, gsl_vector* dydt){
  double yi=gsl_vector_get(y,0);
  double dy=gsl_vector_get(y,1);
  gsl_vector_set(dydt,1,-yi);
  gsl_vector_set(dydt,0,dy);
}

int main(){
  //Harmonic oscillator
  FILE* plot=fopen("Plotdat.txt","w");

  double eps=5e-3, abs=5e-3, a=0, pi=M_PI;
  int k, n=2, m=100;
  gsl_vector* ya=gsl_vector_calloc(2);
  gsl_vector_set(ya,1,1.0);
  double* Yal=malloc(sizeof(double)*(m*n));
  double* Xal=malloc(sizeof(double)*(m));
  
  k=abdriver(harm,a,ya,&Yal,m,&Xal,pi,(pi-a)/30,abs,eps);
  gsl_matrix_view Y=gsl_matrix_view_array(Yal, k, n);
  gsl_vector_view X=gsl_vector_view_array(Xal, k);

  for(int i=0;i<k;i++){
    double xi=gsl_vector_get(&X.vector,i);
    fprintf(plot,"%10g %10g %10g\n",xi,gsl_matrix_get(&Y.matrix,i,0),sin(xi));
  }
  return 0;
}

