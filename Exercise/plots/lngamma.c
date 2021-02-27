#include <math.h>

double LnGamma(double x){
  /// single precision gamma function (Gergo Nemes, from Wikipedia)
  if(x<0)return log(fabs(M_PI/sin(M_PI*x)/exp(LnGamma(1-x))));
  if(x<9)return log(exp(LnGamma(x+1))/x); // move argument up
  double lngamma=x*log(x+1/(12*x-1/x/10))-x+log(2*M_PI/x)/2;
  return lngamma;
}
