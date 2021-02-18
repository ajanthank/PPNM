#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main (int argc, char** argv){
  printf("Program running is: %s\n", argv[0]);
  
  if (argc<2) fprintf(stderr, "There were no arguments\n");
  else {
    for (int i=1; i<argc;i++){
      double x = atof(argv[i]);
      printf("x=%4g  sin(x)=%10g  cos(x)=%10g \n",x,sin(x),cos(x));
    }
  }
  return 0;
}
