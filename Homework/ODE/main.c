#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <assert.h>

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

int driver(void f(double t, gsl_vector* y, gsl_vector* dydt),
	   double a, gsl_vector* ya, double **Yal, int steps, double **Xal, double b, double h, double acc, double eps){
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
	//	printf("Here %d\n",k);
	Yal=realloc(*Yal, sizeof(double)*(n*(k+1)));
	Xal=realloc(*Xal, sizeof(double)*(k+1));
	Y=gsl_matrix_view_array(*Yal, k+1,n);
	//	printf("Here %d\n",k);

	X=gsl_vector_view_array(*Xal,k+1);
	//	printf("Here %d",k);
	//	printf("%g",gsl_vector_get(&X.vector,0));
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

void harm(double t, gsl_vector* y, gsl_vector* dydt){
  double yi=gsl_vector_get(y,0);
  double dy=gsl_vector_get(y,1);
  gsl_vector_set(dydt,1,-yi);
  gsl_vector_set(dydt,0,dy);
}

void SIR15( double t, gsl_vector* y, gsl_vector* dydt){
  double con=15.0,time=10.0, Pop=5.5e6; //time is in units of days
  double Tc=time/con;
  double s=gsl_vector_get(y,0);
  double i=gsl_vector_get(y,1);
  double dsdt=-i*s/(Pop*Tc);
  assert(dsdt<0);
  double didt= i*s/(Pop*Tc)-i/time;
  double drdt=i/time;

  gsl_vector_set(dydt,0,dsdt);
  gsl_vector_set(dydt,1,didt);
  gsl_vector_set(dydt,2,drdt);
}
void SIR10( double t, gsl_vector* y, gsl_vector* dydt){
  double con=10.0,time=10.0, Pop=5.5e6; //time is in units of day
  double Tc=time/con;
  double s=gsl_vector_get(y,0);
  double i=gsl_vector_get(y,1);
  double dsdt=-i*s/(Pop*Tc);
  assert(dsdt<0);
  double didt= i*s/(Pop*Tc)-i/time;
  double drdt=i/time;

  gsl_vector_set(dydt,0,dsdt);
  gsl_vector_set(dydt,1,didt);
  gsl_vector_set(dydt,2,drdt);
}
void SIR5( double t, gsl_vector* y, gsl_vector* dydt){
  double con=5.0,time=10.0, Pop=5.5e6; //time is in units of days
  double Tc=time/con;
  double s=gsl_vector_get(y,0);
  double i=gsl_vector_get(y,1);
  double dsdt=-i*s/(Pop*Tc);
  assert(dsdt<0);
  double didt= i*s/(Pop*Tc)-i/time;
  double drdt=i/time;

  gsl_vector_set(dydt,0,dsdt);
  gsl_vector_set(dydt,1,didt);
  gsl_vector_set(dydt,2,drdt);
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
  
  k=driver(harm,a,ya,&Yal,m,&Xal,pi,(pi-a)/30,abs,eps);
  gsl_matrix_view Y=gsl_matrix_view_array(Yal, k, n);
  gsl_vector_view X=gsl_vector_view_array(Xal, k);

  for(int i=0;i<k;i++){
    double xi=gsl_vector_get(&X.vector,i);
    fprintf(plot,"%10g %10g %10g\n",xi,gsl_matrix_get(&Y.matrix,i,0),sin(xi));
  }

  //SIR Model
  int dk15, dk10, dk5;
  double Pop=5.5e6, da=0., db=100.;
  gsl_vector* dya=gsl_vector_calloc(3);
  double I=661., D0=234318.;
  double vac=410., R=D0+vac;
  double S=Pop-I-R;

  gsl_vector_set(dya,0,S);
  gsl_vector_set(dya,1,I);
  gsl_vector_set(dya,2,R);

  int dn=3, dm=300;

  double* Dal15=malloc(sizeof(double)*(dn*dm));
  double* Tal15=malloc(sizeof(double)*(dm));

  dk15=driver(SIR15,da,dya,&Dal15,dm,&Tal15,db,(db-da)/dm,abs,eps);

  gsl_matrix_view D15=gsl_matrix_view_array(Dal15,dk15,dn);
  gsl_vector_view T15=gsl_vector_view_array(Tal15,dk15);

  fprintf(plot,"\n \n #Here begins the SIR model data, for Population=5.5e6, Time=10days and People within contact 15.\n");

  for(int i=0;i<dk15;i++){
    double ti=gsl_vector_get(&T15.vector,i);
    fprintf(plot,"%10g %10g %10g %10g\n", ti, gsl_matrix_get(&D15.matrix,i,0),gsl_matrix_get(&D15.matrix,i,1),gsl_matrix_get(&D15.matrix,i,2));
  }
  
  double* Dal10=malloc(sizeof(double)*(dn*dm));
  double* Tal10=malloc(sizeof(double)*(dm));

  dk10=driver(SIR10,da,dya,&Dal10,dm,&Tal10,db,(db-da)/dm,abs,eps);

  gsl_matrix_view D10=gsl_matrix_view_array(Dal10,dk10,dn);
  gsl_vector_view T10=gsl_vector_view_array(Tal10,dk10);
 
  fprintf(plot,"\n \n #Here begins the SIR model data, for Population=5.5e6, Time=10days and People within contact 10.\n");
  for(int i=0;i<dk10;i++){
    double ti=gsl_vector_get(&T10.vector,i);
    fprintf(plot,"%10g %10g %10g %10g\n", ti, gsl_matrix_get(&D10.matrix,i,0),gsl_matrix_get(&D10.matrix,i,1),gsl_matrix_get(&D10.matrix,i,2));
  }
  
  double* Dal5=malloc(sizeof(double)*(dn*dm));
  double* Tal5=malloc(sizeof(double)*(dm));

  dk5=driver(SIR5,da,dya,&Dal5,dm,&Tal5,db,(db-da)/dm,abs,eps);

  gsl_matrix_view D5=gsl_matrix_view_array(Dal5,dk5,dn);
  gsl_vector_view T5=gsl_vector_view_array(Tal5,dk5);

  fprintf(plot,"\n \n #Here begins the SIR model data, for Population=5.5e6, Time=10days and People within contact 5.\n");
  for(int i=0;i<dk5;i++){
    double ti=gsl_vector_get(&T5.vector,i);
    fprintf(plot,"%10g %10g %10g %10g\n", ti, gsl_matrix_get(&D5.matrix,i,0),gsl_matrix_get(&D5.matrix,i,1),gsl_matrix_get(&D5.matrix,i,2));
  }
  
  gsl_vector_free(ya);
  gsl_vector_free(dya);
  free(Xal);
  free(Yal);
  free(Tal15);
  free(Dal15);
  free(Tal10);
  free(Dal10);
  free(Tal5);
  free(Dal5);
  fclose(plot);
  
  return 0;
}
