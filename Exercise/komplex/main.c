#include <stdio.h>
#include "komplex.h"
#define TINY 1e-6

int main(){
  komplex a= {1,-2}, b={3,4};


  komplex r = komplex_add(a,b);
  komplex R = {4,6};
  komplex c = komplex_mul(a,b);
  komplex a_con= komplex_conjugate(a);
  double a_abs= komplex_abs(a);
  komplex div= komplex_div(a,b);
  komplex e= komplex_exp(a);
  komplex sina=komplex_sin(a);
  komplex cosa=komplex_cos(a);
  komplex sqrta=komplex_pow(a);
  //komplex powa=komplex_pow(a);

  komplex_print("a=",a);
  komplex_print("b=",b);
  komplex_print("a+b should   = ", R);
  komplex_print("a+b actually = ", r);
  komplex_print("a*b =", c);
  komplex_print("Complex conjugat of a=",a_con);
  printf("Absolute value of a=%g\n",a_abs);
  komplex_print("a/b=",div);
  komplex_print("exp(a)=",e);
  komplex_print("sin(a)=",sina);
  komplex_print("cos(a)=",cosa);
  komplex_print("sqrt(a)=",sqrta);
  //komplex_print("pow(a,1/2)=",powa);
  return 0;
}
