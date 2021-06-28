#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

int binsearch(int n, double* x, double z){/* locates the interval for z by bisection */ 
	assert(n>1 && x[0]<=z && z<=x[n-1]);
	int i=0, j=n-1;
	while(j-i>1){
		int mid=(i+j)/2;
		if(z<x[mid]){
		  j=mid;
		}
		else{
		  i=mid;
		}
	}
	return i;
}

double linterp(int n, double* x, double* y, double z){
  int index=binsearch(n,x,z);
  assert(x[index+1]>x[index]);
  double Delta_x=x[index+1]-x[index];
  double Delta_y=y[index+1]-y[index];
  double p=Delta_y/Delta_x;
  double lin_z=y[index]+p*(z-x[index]);
  
  return lin_z;
}

double linterp_integ(int n, double* x,double* y,double z){
  int i = binsearch(n,x,z);
  double output=0;
  int j=1;
  while(j<=i){
    output+=y[j-1]*(x[j]-x[j-1])+0.5*(y[j]-y[j-1])*(x[j]-x[j-1]);
    j++;
  }
  double yinput=linterp(n,x,y,z);
  output+=y[i]*(z-x[i])+0.5*(yinput-y[i])*(z-x[i]);

  return output;
}

// Quadratic Splin

typedef struct{int n; double *x, *y, *b, *c;} qinterp;

qinterp* qinterp_alloc(int n, double* x, double* y){
  qinterp *s=(qinterp*)malloc(sizeof(qinterp));
  s->b=(double*)malloc((n-1)*sizeof(double));
  s->c=(double*)malloc((n-1)*sizeof(double));
  s->x=(double*)malloc(n*sizeof(double));
  s->y=(double*)malloc(n*sizeof(double));
  s->n=n;

  for(int i=0;i<n;i++){
    s->x[i]=x[i];
    s->y[i]=y[i];
  }

  double dx[n-1], dy[n-1];
  for(int i=0; i<n-1;i++){
    dx[i]=x[i+1]-x[i];
    dy[i]=y[i+1]-y[i];
  }

  s->c[0]=0;
  for(int i=0;i<n-2;i++){
    s->c[i+1]=(dy[i+1]-dy[i]-s->c[i]*dx[i])/dx[i+1];
  }

  s->c[n-2]/=2;
  for(int i=n-3;i>=0;i--){
    s->c[i]=(dy[i+1]-dy[i]-s->c[i+1]*dx[i+1])/dx[i];
  }

  for(int i=0;i<n-1;i++){
    s->b[i]=dy[i]-s->c[i]*dx[i];
  }

  return s;
}

double qieval(qinterp* s, double z){
  int i=binsearch(s->n,s->x,z);
  double dx=z-s->x[i];
  double val=s->y[i]+dx*(s->b[i]+dx*s->c[i]);
  return val;
}

double qiinteg(qinterp* s, double z){
  int i=binsearch(s->n,s->x,z);
  double sum=0, dx=0;
  for(int j=0;j<=i;j++){
    dx=s->x[j]-s->x[j-1];
    sum+=dx*(s->y[j-1]+dx*(s->b[j-1]/2+dx*s->c[j-1]/3));
  }
  dx=z-s->x[i];
  sum+=dx*(s->y[i]+dx*(s->b[i]/2+dx*s->c[i]/3));
  
  return sum;
}

double qideriv(qinterp* s, double z){
  int i=binsearch(s->n,s->x,z);
  double val=s->b[i]+2.*(z-s->x[i])*s->c[i];
  return val*2.;
}

void qinterp_free(qinterp* s){
  free(s->x);
  free(s->y);
  free(s->b);
  free(s->c);
  free(s);
}
