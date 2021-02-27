#include <stdio.h>
#include <pthread.h>
#include <math.h>
void* bar(void* arg){

  double* x = (double*)arg;

  for (int i=0;i<1e8;i++){
    *x=cos(*x);
  }
  
  return NULL;
}

int main(){
  double x=0, y=100, z=-100;

  pthread_t threadx, thready;
  pthread_attr_t* attributes =NULL;
  pthread_create( &threadx, attributes, bar, (void*)&x);
  pthread_create( &thready, attributes, bar, (void*)&y);
  //bar((void*)&x);
  //bar((void*)&y);
  bar((void*)&z);
  void* returnval= NULL;
  pthread_join(threadx,returnval);
  pthread_join(thready,returnval);
  printf("x=%g cos(x)=%g \n",x, cos(x));
  printf("y=%g cos(y)=%g \n",y,cos(y));
  printf("z=%G cos(z)=%g \n",z,cos(z));
  return 0;
}
