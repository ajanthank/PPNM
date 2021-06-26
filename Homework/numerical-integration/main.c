#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <assert.h>

double RAI24 (double fun(double), double x0, double x1, double abs, double eps, double funval2, double funval3, int nrec, int* eta, double* err){
  assert(nrec<1e6);
  double funval1=fun(x0+1*(x1-x0)/6);
  double funval4=fun(x0+5*(x1-x0)/6);

  double Q=(2*funval1+funval2+funval3+funval4*2)*(x1-x0)/6;
  double q=(funval1+funval2+funval3+funval4)*(x1-x0)/4;

  double tol=abs+eps*fabs(Q);
  double Err=fabs(Q-q);

  if(Err<tol){
    *err += Err*Err;
    if(*eta<nrec){
      *eta=nrec;
    }
    return Q;
  }
  else{
    double Q1=RAI24(fun, x0, (x0+x1)/2, abs/sqrt(2.), eps, funval1, funval2, nrec+1,eta, err);
    double Q2=RAI24(fun, (x0+x1)/2, x1, abs/sqrt(2.), eps, funval3, funval4, nrec+1, eta, err);
    return Q1+Q2;
  }
}

double integ(double fun(double), double x0, double x1, double abs, double eps,int* eta, double* err){
  double funval2=fun(x0+2.*(x1-x0)/6.);
  double funval3=fun(x0+3.*(x1-x0)/6);
  int nrec=0;
  *err=0;
  double val=RAI24(fun, x0,x1,abs,eps,funval2, funval3, nrec,eta,err);
  *err=sqrt(*err);

  return val;
}
//Trapezium
double trpz(double trapfun(double fun(double), double min, double max, double x), double fun(double), double trp2, double trp3, double x0, double x1,
	    double min, double max, double abs, double eps, int nrec,int* eta, double* err){

  assert(nrec<1e6);
  double f1=x0+1./6*(x1-x0);
  double f4=x0+5./6*(x1-x0);
  double trp1=trapfun(fun,min,max,f1);
  double trp4=trapfun(fun,min,max,f4);

  double Q=(x1-x0)/6.*(2*trp1+trp2+trp3+2*trp4);
  double q=(x1-x0)/4.*(trp1+trp2+trp3+trp4);
  double tol=abs+eps*fabs(Q);
  double Err=fabs(Q-q);

  if(Err<tol){
    *err+=Err*Err;
    return Q;
  }
  else{
    double Q1=trpz(trapfun, fun, trp1,trp2,x0,(x0+x1)/2,min, max, abs/sqrt(2), eps, nrec+1,eta, err);
    double Q2=trpz(trapfun, fun, trp1,trp2,(x0+x1)/2,x1,min, max, abs/sqrt(2), eps, nrec+1,eta, err);
    return Q1+Q2;
  }
}

double trapfun(double fun(double), double x0, double x1, double x){
  double trans=((x1-x0)/2.*cos(x)+(x1+x0)/2.);
  return fun(trans)*sin(x)*(x1+x0)/2.;
}

double x0x1inf(double fun(double), double x0, double x1, double x){
  double trans=(cos(x)/(1-cos(x)*cos(x)));
  return fun(trans)*sin(x)*(1+cos(x)*cos(x))/pow((1-cos(x)*cos(x)),2);
}

double x0inf(double fun(double), double x0, double x1, double x){
  double trans=x1-((1-cos(x))/(1+cos(x)));
  return fun(trans)* sin(x) * 2./pow(cos(x)+1,2);
}

double x1inf(double fun(double), double x0, double x1, double x){
  double trans=(x0+(cos(x)+1)/(1-cos(x)));
  return fun(trans)* sin(x) * 2./pow(cos(x)-1,2);
}

// Clenshaw-Curtis
double CC(double fun(double), double x0, double x1, double abs, double eps,int* eta, double* err){
  double Q, funval2, funval3;
  double xc2=2./6.*M_PI, xc3=4./6.*M_PI;
  int nrec=0;
  *err=0.;
  if(isinf(x0)==1 && isinf(x1)==1){
    funval2=x0x1inf(fun,0.,0.,xc2);
    funval3=x0x1inf(fun,0.,0.,xc3);
    Q=trpz(x0x1inf,fun,funval2,funval3, 0.,M_PI, x0, x1, abs, eps,nrec,eta,err);
  }
  else if(isinf(x0)==1 && isinf(x1)==0){
    funval2=x0inf(fun,0.,x1,xc2);
    funval3=x0inf(fun,0.,x1,xc3);
    Q=trpz(x0inf, fun, funval2, funval3, 0., M_PI, x0, x1, abs, eps,nrec,eta,err);
  }
  else if(isinf(x0)==0 && isinf(x1)==1){
    funval2=x1inf(fun,x0,0.,xc2);
    funval3=x1inf(fun,x0,0.,xc3);
    Q=trpz(x1inf, fun, funval2, funval3,0.,M_PI, x0, x1, abs, eps, nrec,eta, err);
  }
  else{
    funval2=trapfun(fun,x0,x1,xc2);
    funval3=trapfun(fun,x0,x1,xc3);
    Q=trpz(trapfun, fun, funval2, funval3, 0., M_PI, x0, x1, abs, eps, nrec,eta, err);
  }
  *err=sqrt(*err);
  return Q;
}

//###################### Function for excersie ########################

double fun1(double x){
   return sqrt(x);
}

double fun2(double x){
  return 4*sqrt(1-pow(x,2));
}

double fun3(double x){
  return 1/sqrt(x);
}

double fun4(double x){
  return log(x)/sqrt(x);
}

double fun5(double x){
  return 4*sqrt(1-pow(x,2));
}

double gslfun5(double x, void* params){
  params=NULL;
  return fun5(x);
}

double fun6(double x){
  return exp(-x*x);
}

double gslfun6(double x, void* params){
  params=NULL;
  return fun6(x);
}

//######################### Print and main function ###########################

void printr(char* text, double integval, double exact, double abs, double eps, double err, int nrec){
  printf("Integration function:\t %s \nExact solution:\t %10g \nNumerical Solution:\t %10g\nError Goal:\t %10g\nActual Error:\t %10g\nComputational Error:\t %10g\nFunction calls:\t %d\n\n\n",text, exact,integval, abs+fabs(exact)*eps, fabs(integval-exact),err, nrec);
}


int main(){
  double x0=0., x1=1.; 
  double abs=1e-3, eps=1e-3, err=0;

  int calls=0;

  printf("-------------------- Exercise A -------------------\n\n");

  double integvalfun1= integ(fun1, x0,x1,abs,eps,&calls, &err);
  printr("∫_0^1 dx √(x)",integvalfun1,2./3.,abs,eps,err,calls);
  
  double integvalfun2= integ(fun2,x0,x1,abs,eps,&calls ,&err);
  printr("∫_0^1 dx 4√(1-x²)",integvalfun2,M_PI,abs,eps,err,calls);

  printf("------------------- Exercise B --------------------\n\n");

  printf("######################################################\n\n");

  double integvalfun3=integ(fun3,x0,x1,abs,eps,&calls,&err);
  printf("The first integrater:\n");
  printr("∫_0^1 dx 1/√(x)",integvalfun3,2.,abs,eps,err,calls);
  printf("The CC-integrater:\n");
  double CCintegvalfun3= CC(fun3,x0,x1,abs,eps,&calls,&err);
  printr("∫_0^1 dx 1/√(x)",CCintegvalfun3,2.,abs,eps,err,calls);

  printf("######################################################\n\n");

  double integvalfun4=integ(fun4,x0,x1,abs,eps,&calls,&err);
  printf("The first integrater:\n");
  printr("∫_0^1 dx 1/√(x)",integvalfun4,-4.,abs,eps,err,calls);
  printf("The CC-integrator:\n");
  double CCintegvalfun4=CC(fun4,x0,x1,abs,eps,&calls,&err);
  printr("∫_0^1 dx 1/√(x)",CCintegvalfun4,-4.,abs,eps,err,calls);

  printf("######################################################\n\n");

  printf("The first integrator:\n");
  double integvalfun5=integ(fun5,x0,x1,abs,eps,&calls,&err);
  printr("∫_0^1 dx 4√(1-x²)",integvalfun5,M_PI,abs,eps,err,calls);
  printf("The CC-integrator:\n");
  double CCintegvalfun5=CC(fun5,x0,x1,abs,eps,&calls,&err);
  printr("∫_0^1 dx 4√(1-x²)",CCintegvalfun5,M_PI,abs,eps,err,calls);

  gsl_function gsl_fun5;
  gsl_fun5.function=&gslfun5;
  gsl_fun5.params=NULL;
  gsl_integration_cquad_workspace* workspace=gsl_integration_cquad_workspace_alloc(1000);

  double gslintegvalfun5;
  double gslabs;
  size_t neval;

  printf("GSL CC integrator:\n");

  gsl_integration_cquad(&gsl_fun5, x0, x1, abs, eps, workspace, &gslintegvalfun5, &gslabs, &neval);
  printr("∫_0^1 dx 4√(1-x²)",gslintegvalfun5,M_PI,abs,eps,gslabs,(int)neval);

  printf("----------------------- Exercise C ------------------------\n");
  printf("Both limits are infinity:\n");
  printf("######################################################\n\n");

  printf("CC integrator:\n");
  double CCintegvalfun61=CC(fun6,-INFINITY, INFINITY,abs,eps,&calls,&err);
  printr("∫_(-∞)^(∞) dx exp(-x*x)",CCintegvalfun61,sqrt(M_PI),abs,eps,err,calls);

  gsl_function gsl_fun6;
  gsl_integration_workspace* workspace2=gsl_integration_workspace_alloc(1000);
  gsl_fun6.function=&gslfun6;
  gsl_fun6.params=NULL;
  double qagiintegvalfun6;
  gslabs=0;
  neval=0;

  printf("Gsl integrator:\n");
  gsl_integration_qagi(&gsl_fun6,abs,eps,1000,workspace2,&qagiintegvalfun6,&gslabs);
  printr("∫_(-∞)^(∞) dx exp(-x*x)",qagiintegvalfun6,sqrt(M_PI),abs,eps,gslabs,(int)neval);
  printf("Number of calls is not available.\n");

  printf("First limit is infinity:\n");
  printf("######################################################\n\n");

  printf("CC integrator:\n");
  double CCintegvalfun62=CC(fun6,0,INFINITY,abs,eps,&calls,&err);
  printr("∫_(-∞)^0 dx exp(-x*x)",CCintegvalfun62,sqrt(M_PI)/2.,abs,eps,err,calls);

  printf("Gsl integrator");
  double qagilintegvalfun6;
  gslabs=0;
  neval=0;
  gsl_integration_qagil(&gsl_fun6,0,abs,eps,1000, workspace2, &qagilintegvalfun6,&gslabs);
  printr("∫_(-∞)^0 dx exp(-x*x)",qagilintegvalfun6,sqrt(M_PI)/2.,abs,eps,gslabs,(int)neval);
  printf("Number of calls is not available.\n");

  printf("Second limit is infinity:\n");
  printf("######################################################\n\n");

  printf("CC integrator:\n");
  double CCintegvalfun63=CC(fun6,0,INFINITY,abs,eps,&calls,&err);
  printr("∫_0^(∞) dx exp(-x*x)",CCintegvalfun63,sqrt(M_PI)/2.,abs,eps,err,calls);

  printf("Gsl integrator");
  double qagiuintegvalfun6;
  gslabs=0;
  neval=0;
  gsl_integration_qagiu(&gsl_fun6,0,abs,eps,1000, workspace2, &qagiuintegvalfun6,&gslabs);
  printr("∫_0^(∞) dx exp(-x*x)",qagiuintegvalfun6,sqrt(M_PI)/2.,abs,eps,gslabs,(int)neval);
  printf("Number of calls is not available.\n");


  gsl_integration_cquad_workspace_free(workspace);
  gsl_integration_workspace_free(workspace2);
  return 0;
}
