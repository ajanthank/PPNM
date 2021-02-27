#include<stdio.h>
#include<math.h>
#include<gsl/gsl_integration.h>

double f (double x, void * params){
  double z = *(double *) params;
  double f = log(z*x) / sqrt(x);
  return f;
}

double erfg (double x, void * params){
  double z=*(double *) params;
  double erfg = (2/sqrt(M_PI))*exp(pow(x*z,2));
  return erfg;
}
double myG (double z){
  gsl_function F;
  F.function = &erfg;
  F.params = &z;
  int limit=999;
  //  gsl_set_error_handler_off()
  gsl_integration_workspace* w = gsl_integration_workspace_alloc (limit);
  double a=0,acc=1e-6,eps=1e-6,result,error;

  gsl_integration_qagiu(&F,a,acc,eps,limit,w,&result,&error);

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

  gsl_integration_qags (&F, 0, 1, 0, 1e-8, 1000,w, &result, &error);

  printf("Result= %g\n", result);

  FILE * datab=fopen("DataB.txt","w");
  for (double x= 0; x<=10;x=+1.0){
    fprintf(datab,"%10g %10g",x,myG(x));
  }


  gsl_integration_workspace_free (w);
  return 0;
}