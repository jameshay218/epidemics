#ifndef SIR_HPP
#define SIR_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include "gsl-1.16/multimin/gsl_multimin.h"
#include <iomanip>

using namespace std;

class SIR{
private:
  double t, step;
  double S,I,R,Pop[3];
  double dPop[3];

  double beta, gamma;
  double tmax;
  vector<vector<double> > current_data;
public:
  SIR();
  SIR(double parameters[7], vector<vector<double> > x);
  ~SIR();

  void Diff(double Pop[3]);
  void Runge_Kutta();
  void Solve_Eq(vector<vector<double> >& data);
  double sse_sir(double parameters[3]);
  double add_arrays(vector<vector<double> > data1, vector<vector<double> > data2);
  void update_params(double parameters[7]);
  void user_params();
  void fit_model();
  double operator()(vector<double> param);
  void update_data(vector<vector<double> > x);
  vector<vector<double> > sse_sir2(double parameters[3]);
  vector<double> rand_params();

};
#endif
