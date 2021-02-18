#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main(){
  double x;
  int items;
  FILE* my_out_stream=fopen("out_file.txt","w");
  do {
    items=fscanf(stdin,"%lg", &x);
    fprintf(my_out_stream, "x=%4g  sin(x)=%10g  cos(x)=%10g \n",x,sin(x),cos(x));
  }while(items!= EOF);
  fclose(my_out_stream);
  return 0;
}
