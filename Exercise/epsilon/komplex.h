#ifndef HAVE_KOMPLEX_H /* this is necessary for multiple includes */
#define HAVE_KOMPLEX_H
#include <complex.h>

struct komplex {double re; double im;};
typedef struct komplex komplex;

void komplex_print(char *s, komplex z){
  printf("%s (%g,%g)\n",s,a.re, a.im);
}
void komplex_set(komplex *z, double x, double y){
  (*z).re = x;
  (*z).im = y;
}
