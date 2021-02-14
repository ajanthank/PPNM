#include <stdio.h>
#include "komplex.h"
#include <math.h>

void komplex_print (char *s, komplex a) {
	printf ("%s (%g,%g)\n", s, a.re, a.im);
}

komplex komplex_new (double x, double y) {
	komplex z = { x, y };
	return z;
}

void komplex_set (komplex* z, double x, double y) {
	(*z).re = x;
	(*z).im = y;
}

komplex komplex_add (komplex a, komplex b) {
	komplex result = { a.re + b.re , a.im + b.im };
	return result;
}

komplex komplex_sub(komplex a, komplex b) {
  komplex result = {a.re - b.re, a.im - a.im};
  return result;
}

komplex komplex_mul (komplex a, komplex b) {
  komplex result = {a.re*b.re-a.im*b.im , a.re*b.im+b.re*a.im};
  return result;
}
komplex komplex_conjugate(komplex z) {
  komplex result = {z.re, -z.im};
  return result;
}
double komplex_abs(komplex z) {
  double result = z.re*z.re-z.im*komplex_conjugate(z).im;
  return result;
}
komplex komplex_div (komplex a, komplex b) {
  b = komplex_conjugate(b);
  komplex result = {(a.re*b.re-a.im*b.im)/komplex_abs(b),(a.re*b.im+b.re*a.im)/komplex_abs(b) };
return result;
}
komplex komplex_exp (komplex z){
  komplex result = {exp(z.re)*cos(z.im),exp(z.re)*sin(z.im)};
  return result;
}
komplex komplex_sin (komplex z) {
  komplex result = {sin(z.re)*cosh(z.im), cos(z.re)*sinh(z.im)};
  return result;
}
komplex komplex_cos (komplex z) {
  komplex result = {cos(z.re)*cosh(z.im), -sin(z.re)*sinh(z.im)};
  return result;
}
komplex komplex_sqrt (komplex z) {
  komplex result = {sqrt((sqrt(pow(z.re,2)+pow(z.im,2))+z.re)/2) , sqrt((sqrt(pow(z.re,2)+pow(z.im,2))-z.re)/2) };
  return result;
}
komplex komplex_pow (komplex z, double n) {
  double r = pow(z.re,2)+pow(z.im,2);
  komplex result = {cos(n*atan(-z.im/z.re))*pow(r,n/2) , sin(n*atan(-z.im/z.re))*pow(r,n/2) };
  return result;
}
