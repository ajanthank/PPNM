#include <math.h>
#include <stdio.h>
#include <assert.h>

int binsearch(int n, double* x, double z){/* locates the interval for z by bisection */ 
	assert(n>1 && x[0]<=z && z<=x[n-1]);
	int i=0, j=n-1;
	while(j-i>1){
		int mid=(i+j)/2;
		if(z>x[mid]) i=mid; else j=mid;
		}
	return i;
}

double linterp(int n, double x[], double y[], double z){
  int index=binsearch(n,x,z);
  assert(x[index+1]>x[index]);
  double Delta_x=x[index+1]-x[index];
  double Delta_y=y[index+1]-y[index];
  double p=Delta_y/Delta_x;
  double lin_z=y[index]+p*(z-x[index]);
  
  return lin_z;
}

double linterp_integ(int n, double x[],double y[],double z){
  int i = binsearch(n,x,z);
  double output=0;
  int j=1;
  while(j<=i){
    output+=y[p-1]*(x[p]-x[p-1])+0.5*(y[p]-y[p-1])*(x[p]-x[p-1]);
    j++;
  }
  double yinput=linterp(n,x,y,z);
  output+=y[i]*(z-x[i])+0.5*(yinput-y[i])*(z-x[i]);

  return output;
}
