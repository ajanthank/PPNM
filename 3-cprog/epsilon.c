#include <stdio.h>
#include <limits.h>
#include <float.h>

int E1(){
  //Exercise i
  printf("Exercise 1.i:\n\n");
  int i,j,k=1;
  int val1=1;//INT_MAX;
  int val2=1;//INT_MIN;
  while	(i<val1){
    i++;
  }
  printf("My value of i for the \"while\" loop is %i:\n",i);
  for (j=1; j<val1; j=j+1){
  }
  printf("My value of i for the \"for\" loop is %i:\n",j);
  do {
    k=k+1;
  } while (k<val1);
  printf("My value of i for the \"do-while\" loop is %i:\n",k);

  printf("The value of INT_MAX in limits.h is %i:\n", INT_MAX);

  i=j=k=1;

  //Exercise ii
  printf("\n\nExercise 1.ii\n\n");

  while (i>val2){
    i--;
  }
  printf("My value of i for the \"while\" loop is %i:\n",i);
  for (j=1; j>val2; j=j-1){
  }
  printf("My value of i for the \"for\" loop is %i:\n",j);
  do {
    k=k-1;
  } while (k>val2);
  printf("My value of i for the \"do-while\" loop is %i:\n",k);

  printf("The value of INT_MAX in limits.h is %i:\n", INT_MIN);

  //Exercise iii
  printf("\n\nExercise 1.iii\n\n");

  double x = 1;

  while (1+x!=1){
    x/=2;
  }
  x*=2;
  
  printf("%f\n",x);
  printf("%f \n%f \n%Lf \n",FLT_EPSILON, DBL_EPSILON, LDBL_EPSILON);
  return 0;
}

int E2(){
  int max = INT_MAX/3;
  printf("max=%i \n", max);
  int i=1;
  float sum_up_float, sum_down_float=0.0;
  while(i<max){
    sum_up_float = 1.f/i+sum_up_float;
    i++;
  }
  
  printf("The sum_up_float becomes: %f \n",sum_up_float);

  while(i>0){
    sum_down_float = sum_down_float+1.f/i;
    i--;
  }
  printf("The sum_down_float becomes: %f \n", sum_down_float);


  double sum_up_double, sum_down_double=0.0;

  i=1;

  while(i<max){
    sum_up_double=sum_up_double+1.f/i;
    i++;
  }
  while(i>0){
    sum_down_double=sum_down_double+1.f/i;
    i--;
  }
  printf("The sum_up_double is: %f \nThe sum_down_double is: %f \n",sum_up_double, sum_down_double);
  return 0;
}

int equal(double a, double b, double epsilon, double tau);

void name_digit(int i);

int math_exercise();

int main(){
  printf("\n\n\n################## Exercise 1 #################\n");
  E1();
  printf("\n\n\n################## Exercise 2 #################\n");
  E2();
  printf("\n\n\n################## Exercise 3 #################\n");
  printf("The return value of the function is: %i\n", equal(2.3,2.5,0,0));
  printf("\n\n\n################## Exercise 4 #################\n");
  name_digit(12);
  math_exercise();
  return 0;
}
