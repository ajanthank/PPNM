#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_interp.h>

typedef struct{int n; double *x, *y, *b, *c;} qinterp;

double linterp(int n, double x[], double y[], double z);
double linterp_integ(int n, double x[], double y[], double z);
qinterp* qinterp_alloc(int n, double* x, double* y);
double qieval(qinterp* s, double z);
double qiinteg(qinterp* s, double z);
double qideriv(qinterp* s, double z);
void qinterp_free(qinterp* s);

int main(){

  //printf("###################### Problem A ######################\n");
  
  int n=20;

  double x[n], y[n], dy[n], x_min=0,x_max=8;

  FILE* xydata=fopen("XY.txt","w");

  for (int i=0;i<n;i++){
    x[i]=i/2.0;
    y[i]=cos(x[i]);
    dy[i]=-sin(x[i]);
    fprintf(xydata,"%10g %10g\n", x[i], y[i]);
  }
  double dx= 1.0/(320);
  FILE* lindata=fopen("LinData.txt","w");
  for (double i=x_min; i<=x_max;i+=dx){
    fprintf(lindata,"%10g %10g\n",i,linterp(n,x,y,i));
  }

  FILE* lintegdata = fopen("LinIntegData.txt","w");
  for (double i=x_min; i<=x_max;i+=dx){
    fprintf(lintegdata,"%10g %10g\n",i,linterp_integ(n,x,y,i)+sin(x_min));
    }
  
  fclose(lindata);
  fclose(lintegdata);

  printf("###################### Problem B ######################\n");
  printf("#The quadspline:\n");

  qinterp* s=qinterp_alloc(n,x,y);

  int z=0;
  double eps=10;

  while(z<=eps*s->x[n-1]){
    printf("%10g %10g %10g %10g %10g %10g %10g\n",z/eps,qieval(s,z/eps),qiinteg(s,z/eps),qideriv(s,z/eps),cos(z/eps),sin(z/eps),
	   -sin(z/eps));
    z++;
  }

  qinterp_free(s);
  fclose(xydata);
  return 0;
}
