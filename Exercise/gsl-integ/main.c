#include<stdio.h>
#include<math.h>
#include<gsl/gsl_integration.h>

double f (double x, void * params){
  double z = *(double *) params;
  double f = log(z*x) / sqrt(x);
  return f;
}

double erfg (double x, void * param){

  double erfg = (2/sqrt(M_PI))*exp(-pow(x,2));
  return erfg;
}
double myG (double z){
  gsl_function G;
  G.function = &erfg;

  int limit=999;

  gsl_integration_workspace* w = gsl_integration_workspace_alloc (limit);
  double a=0,result,error;

  gsl_integration_qags(&G,a,z,0,1e-8,limit,w,&result,&error);

  gsl_integration_workspace_free (w);
  return result;
}


int main (){
  int limit = 1000;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (limit);

  double result, error;
  double z=1;

  gsl_function F;
  F.function=&f;
  F.params = &z;

  gsl_integration_qags (&F, 0, 1, 0, 1e-8, limit,w, &result, &error);

  printf("Result= %g\n", result);

  FILE * datab=fopen("DataB.txt","w");
  for (double x= -3; x<=3;x+=1.0/100){
    fprintf(datab,"%10g %10g\n",x,myG(x));
  }

  gsl_integration_workspace_free (w);
  fclose(datab);

  return 0;
}
