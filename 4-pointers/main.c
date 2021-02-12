#include <stdio.h>
#include <stdlib.h>

//struct vector {int n, double v[]};
//struct vector my_vector;
//typedef struct vector vector;
typedef struct {int n, double data[]} vector;
vector* vector_allocate(int n){
  vector* v=malloc(sizeof(vector));
  (*v).n=n;
  (*v).data=mallac(n*sizeof(double));
  return v;
}
void vector_free(vector* v){
  free((*v.data));
  free(v);
}

void vector_set(vector* v, int i, double value){
  if (i<0)
    printf()
}

void set0(double x){
  x=0;
}

void set0cor(double *x){
  (*x)=0;
}
int main(){
  double y=1;

  set0(y);
  printf("y after set0 = %g \n",y);

  set0cor(&y);
  printf("y after set0cor = %g \n", y);

  int n = 5;

  double v[5];
  for (int i=0;i<n;i++){
    v[i]=i;
  }

  int i = 0;
  while (i<n) {
    printf("v[%d]=%g \n", i, v[i]);
    i++;
  }
  double w[5]={1,2,3,4,5};
  printf("%g\n", );
  return 0;
 free 
}
