#include<stdio.h>
#include<math.h>


int main() {
  FILE *input_file= fopen("input.txt","r");
  FILE *output_file = fopen("output2.txt","w");

  double x;

  int i=0;

  while(i!=EOF){
    i = fscanf(input_file,"%lg",&x);
    fprintf(output_file, "x=%4g  sin(x)=%10g  cos(x)=%10g\n",x,sin(x),cos(x));
  }

  fclose(input_file);
  fclose(output_file);
  return 0;
}
