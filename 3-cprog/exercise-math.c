#include <stdio.h>
#include <math.h>
#include <complex.h>

int math_exercise () {
  printf("############################################################################# \n");
  printf("Here begins the Math-exercises.\n");
  printf("############################################################################# \n");
  char gamma[]= "\u0393";
  printf("This is exercise 1:\n");
  printf("This is %s of 5: %f \n", gamma, tgamma(5));
  printf("This is J1 of 3: %f \n", j1(3));

  double complex z = csqrt(-2.0);
  printf("This is sqrt(-2):                  Re( %g),  Im( %g) \n", creal(z), cimag(z));
  double complex z2 = cexp(I*M_PI);
  printf("This is exp(i*pi):                 Re( %g),  Im( %g) \n", creal(z2), cimag(z2));
  double complex z3 = cexp(I*1);
  printf("This is exp(i):                    Re( %g),  Im( %g) \n", creal(z3), cimag(z3));
  double complex z4 = cpow(I,M_E);
  printf("This is i to the power of e:       Re( %g),  Im( %g) \n", creal(z4), cimag(z4));
  double complex z5 = cpow(I,I);
  printf("This is i to the power of i:       Re( %g),  Im( %g) \n", creal(z5), cimag(z5));

  float x_f=1.f/9;
  double x_d=1./9;
  long double x_ld=1.L/9;

  printf("This is exercise 2: \n");
  printf("Float:              %0.25gf \n",x_f );
  printf("Double:             %0.25lg \n",x_d);
  printf("Long Double:        %0.25Lg \n", x_ld);
  return 0;
}
