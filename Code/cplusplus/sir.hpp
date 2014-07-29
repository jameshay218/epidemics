#ifndef SIR_HPP
#define SIR_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

class SIR{
private:
  double t, step;
  double S,I,R,Pop[3];
  double dPop[3];

  double beta, gamma;
  double tmax;

public:
  SIR();
  SIR(double parameters[7]);
  ~SIR();

  void Diff(double Pop[3]);
  void Runge_Kutta();
  void Solve_Eq(vector<vector<double> >& data);
  void sse_sir(double parameters[3], vector<vector<double> > data);
  double add_arrays(vector<vector<double> > data1, vector<vector<double> > data2);
  void update_params(double parameters[7]);

  void user_params();

};
#endif
