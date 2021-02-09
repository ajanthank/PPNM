#include <stdio.h>
#include <math.h>
#include <complex.h>

int math_exercise () {
  printf("############################################################################# \n");
  printf("Here begins the Math-exercises.\n");
  printf("############################################################################# \n");
  
  printf("This is exercise 1:\n");
  printf("This is gamma of 5: %f \n", tgamma(5));
  printf("This is J1 of 3: %f \n", j1(3));
  printf("This is sqrt(-2):                  Re( %f),  Im( %f) \n", creal(sqrt(-2)), cimag(sqrt(-2)));
  printf("This is exp(i*pi):                 Re( %f),  Im( %f) \n", creal(exp(I*M_PI))), cimag(exp(I*M_PI));
  printf("This is exp(i):                    Re( %f),  Im( %f) \n", creal(exp(I))), cimag(exp(I));
  printf("This is i to the power of e:       Re( %f),  Im( %f) \n", creal(pow(I,M_E))), cimag(pow(I,M_E));
  printf("%f\n %f \n%f\n%f\n", creal(I), cimag(I), M_E, M_PI);

  float x_f=1.f/9;
  double x_d=1./9;
  long double x_ld=1.L/9;

  printf("This is exercise 2: \n");
  printf("Float:              %0.25gf \n",x_f );
  printf("Double:             %0.25lg \n",x_d);
  printf("Long Double:        %0.25Lg \n", x_ld);
  return 0;
}
