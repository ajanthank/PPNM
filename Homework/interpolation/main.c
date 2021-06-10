#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_interp.h>

double linterp(int n, double x[], double y[], double z);
double linterp_integ(int n, double x[], double y[], double z);

int main(){
  printf("###################### Problem A ######################\n");
  printf("--------------------- Question 1 ---------------------\n");
  
  int n=20;

  double x[n], y[n],x_min=-4,x_max=4, step=(x_max-x_min)/(n-1);

  FILE* xydata=fopen("XY.txt","w");

  for (int i=0;i<n;i++){
    x[i]=x_min+ (double) i *step;
    y[i]=cos(x[i]);
    fprintf(xydata,"%10g %10g\n", x[i], y[i]);
  }
  double dx= 1.0/320;
  FILE* lindata=fopen("LinData.txt","w");
  for (double i=x_min; i<=x_max;i+=dx){
    fprintf(lindata,"%10g %10g\n",i,linterp(n,x,y,i));
  }

  FILE* lintegdata = fopen("LinIntegData.txt","w");
  for (double i=x_min; i<=x_max;i+=dx){
    fprintf(lintegdata,"%10g %10g\n",i,linterp_integ(n,x,y,i)+sin(x_min));
  }

  double integral=linterp_integ(n,x,y,M_PI);
  printf("Estimate: %10g\n",integral);
  printf("Real: %10g \n",sin(4));
  
  fclose(xydata);
  fclose(lindata);
  fclose(lintegdata);
  return 0;
}
