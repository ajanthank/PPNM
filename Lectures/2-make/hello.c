#include <stdio.h>
#include <stdlib.h>
#include <pwd.h>
#include <math.h>

int main(){

  printf("Hello world and hello %s\n", getenv("USER"));
  
return 0;
}
