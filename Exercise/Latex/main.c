#include<stdio.h>
#include<math.h>

double ex(double x){
if(x<0)return 1/ex(-x);
if(x>1./8)return pow(ex(x/2),2);
return 1+x*(1+x/2*(1+x/3*(1+x/4*(1+x/5*(1+x/6*(1+x/7*(1+x/8*(1+x/9*(1+x/10)))))))));
}

double e(double x){
if(x<0)return 1/e(-x);
return 1+x*(1+x/2*(1+x/3*(1+x/4*(1+x/5*(1+x/6*(1+x/7*(1+x/8*(1+x/9*(1+x/10)))))))));
}


int main(){
  double dx=1.0/8;
  for (double i=-1.0;i<=10.0;i+=dx){
    printf("%10g %10g %10g %10g\n",i,exp(i),ex(i),e(i));
  }
  return 0;
}
