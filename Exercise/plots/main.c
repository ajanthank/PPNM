#include<stdio.h>
#include<math.h>
#include<gsl/gsl_sf_erf.h>
#include<gsl/gsl_sf_gamma.h>
#include<complex.h>

double Erf(double);
double Gamma(double);
double LnGamma(double);
double complex cgamma(double complex z);

int main(){

  double xmin=-2,xmax=2;
  for (double x = xmin; x<=xmax; x+=1.0/8){
    printf("%10g %10g %10g %10g\n",x,erf(x),gsl_sf_erf(x),Erf(x));
  }
  FILE* data2=fopen("data2.txt","w");
  xmin=-4+1.0/32;
  xmax=6;
  for (double x = xmin;x<=xmax; x+=1.0/16){
    fprintf(data2,"%10g %10g %10g %10g \n",x,tgamma(x),gsl_sf_gamma(x),Gamma(x));
  }
  FILE* data3=fopen("data3.txt","w");
  xmin=-5.01;
  xmax=10;
  for (double x = xmin;x<=xmax;x+=1.0/900) {
    fprintf(data3,"%10g %10g %10g %10g \n",x,lgamma(x),gsl_sf_lngamma(x),LnGamma(x));
  }
  double eps=1.0/128 , d_re=1.0/64 , d_im=1.0/64;

  FILE* data4=fopen("data4.txt","w");

  for( double re=-4.5+eps;re<=4+eps;re+=d_re){
    for( double im=-4.5;im<=4.5;im+=d_im){
      double complex z= CMPLX(re,im);
      fprintf(data4,"%10g %10g %10g\n",creal(z),cimag(z),cabs(cgamma(z)));
    }
  }
  

  fclose(data2);
  fclose(data3);
  fclose(data4);
  return 0;
}
