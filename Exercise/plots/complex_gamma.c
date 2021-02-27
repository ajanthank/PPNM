#include<stdio.h>
#include<math.h>
#include<complex.h>

double complex cgamma(double complex z){
  double complex Pi = CMPLX(M_PI,0);
  if (creal(z)<0) {
    return Pi/csin(Pi*z)/cgamma(CMPLX(1.0,0.0)-z);
  }
  if (creal(z)<9) {
    return cgamma(z+CMPLX(1.0,0.0))/z;
  }
  double complex a=CMPLX(1,0);
  double complex b=CMPLX(12,0);
  double complex c=CMPLX(10,0);
  double complex d=CMPLX(2,0);
  double complex lncgamma=z*clog(z+a/(b*z-a/z/c))-z+clog(d*Pi/z)/d;
  return cexp(lncgamma);
}
