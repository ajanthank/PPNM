#include<stdio.h>
#include<math.h>

int main(){
  double x;
  int items;
  FILE* my_out_stream=fopen("out.file.txt","w");
  do{
    items=scanf("%lg",&x);
    printf("x=%g sin(x)=%g\n", x, sin(x));
    fprintf()
  } while (items!=EOF);
  return 0;
}
