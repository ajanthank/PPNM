#include <stdio.h>

double f(int i){
  double r= 1.0/i;
  return r;
}

int main(){
  int n= 1e7;
  printf("n=%g\n",(double) n);
  double s=0;
  for( int i=1; i<=n; i=i+1)
    s+=f(i);
  printf("s=%g\n",s);
}
