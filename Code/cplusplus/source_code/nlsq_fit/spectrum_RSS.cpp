#include "data.h"
#include <cmath>
#include <iostream>
double spectrum_RSS(std::vector<double> arg){
  using namespace Measurement;
  using namespace std;
  double I=arg.at(0);
  double Nw=arg.at(1);
  double lambda_r=arg.at(2);
  double M=arg.at(3);
  double Nd=arg.at(4);
  double h=arg.at(5);
  double PI=3.141592653589793;
  cout.precision(12);

  size_t size=dataX.size();
  vector<double> y(size, 0);
  double RSS=0;

  for(int i=0; i<size; ++i){
    y[i]=I*pow(sin(PI*Nw*(dataX[i]-lambda_r)/lambda_r)/
    	       (PI*Nw*(dataX[i]-lambda_r)/lambda_r)
    	       ,2)*
      ( 1 + M*cos(2*PI*(Nw+Nd)*lambda_r/dataX[i]) )+h;
  }

  for(int i=0; i<size; ++i)
    RSS += pow(y[i]-dataY[i],2);

  return RSS;
}
