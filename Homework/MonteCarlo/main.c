#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#define RND ((double) rand()/RAND_MAX)

double corput(int n, int base){
  double q=0, bk=(double) 1./base;
  while (n>0){
    q+=(n%base)*bk;
    n/=base;
    bk/=base;
  }
  return q;
}

void halton1(int n, int d, double *x){
  int base[]={2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53};
  int maxd=sizeof(base)/sizeof(int);
  assert(d<=maxd);
  for(int i=0;i<d;i++){
    x[i]=corput(n+1,base[i]);
  }
}

void halton2(int n, int d, double *x){
  int base[]={3,5,7,11,13,17,19,23,31,37,41,43,47,53,59,61,67};
  int maxd=sizeof(base)/sizeof(int);
  assert(d<=maxd);
  for(int i=0;i<d;i++){
    x[i]=corput(n+1,base[i]);
  }
}

void randnum(int d, const double *x0, const double *x1, double *x){
  for(int i=0;i<d;i++){
    x[i]=x0[i]+RND*(x1[i]-x0[i]);
  }
}

void randnum_halton_cor(int n, int d, const double *x0, const double *x1, double *x){
  halton1(n,d,x);
  for(int i=0;i<d;i++){
    x[i]=x0[i]+x[i]*(x1[i]-x0[i]);
  }
}

void randnum_halton2_cor(int n, int d, const double *x0, const double *x1, double *x){
  halton2(n,d,x);
  for(int i=0;i<d;i++){
    x[i]=x0[i]+x[i]*(x1[i]-x0[i]);
  }
}

void MC(int d, double *x0, double *x1, double fun(double* x), int points, double* val, double* err){
  double V=1;
  for(int i=0;i<d;i++){
    V*=x1[i]-x0[i];
  }
  double sum=0;
  double sum2=0;
  double funval, x[d];

  for(int i=0;i<points;i++){
    randnum(d, x0,x1,x);
    funval=fun(x);
    sum+=funval;
    sum2+=pow(funval,2);
  }

  double mean=sum/points;
  double variance=sum2/points-pow(mean,2);
  *val=mean*V;
  *err=sqrt(variance/points)*V;
}

void quasiMC(int d, double *x0, double *x1, double fun(double* x), int points, double* val, double* err){
  double V=1;
  for(int i=0;i<d;i++){
    V*=x1[i]-x0[i];
  }

  double sum1=0, sum2=0, funval1, funval2, x1x[d], x2x[d];
  for(int i=0;i<points/2;i++){
    randnum_halton_cor(i,d,x0,x1,x1x);
    randnum_halton2_cor(i,d,x0,x1,x2x);
    funval1=fun(x1x);
    funval2=fun(x2x);
    if(!isinf(funval1) && !isinf(funval2)){
      sum1+=funval1;
      sum2+=funval2;
    }
  }
  double mean=(sum2+sum1)/points;

  *val=mean*V;
  *err=V*fabs(sum1-sum2)/points;
}

double fun1(double* x){
  return 1/(M_PI*M_PI*M_PI*(1-cos(x[0])*cos(x[1])*cos(x[2])));
}

int main(){

  printf("------------------- Exercise A ----------------\n");

  int points=1e5;

  const int d=3;
  double* x0=calloc(d,sizeof(double));
  double* x1=malloc(sizeof(double)*d);
  for(int i=0;i<d;i++){
    x0[i]=0;
    x1[i]=M_PI;
  }
  double val=0, err=0, valq=0, errq=0;

  MC(d,x0,x1,fun1,points,&val,&err);
  quasiMC(d,x0,x1,fun1,points,&valq,&errq);

  printf("Pseudo-random MC:\nResult of MC: %g\nComputed Error:%g\n",val,err);

  printf("--------------------- Exercise B --------------------\n");

  printf("Quasi-random MC:\nResult of MC: %g\nComputed Error:%g\n",valq,errq);
  
  FILE* dat=fopen("Data.txt", "w");
  int Nmax=1e3, steps=1e3;
  for (int i=1; i<Nmax; i++){
    int N=i*steps;
    MC(d, x0, x1, fun1, N, &val,&err);
    quasiMC(d,x0,x1,fun1,N, &valq, &errq);
    fprintf(dat, "%d %10g %10g\n", N, err, errq);
  }

  fclose(dat);
  return 0;
}
