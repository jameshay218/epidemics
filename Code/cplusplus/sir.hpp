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
  double t, step, t0;
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
  void Solve_Eq_t0(vector<vector<double> >& _results);
  void Solve_Eq_total(vector<vector<double> >& _results);


  double sse_sir(vector<double> parameters);
  double sse_sir_t0(vector<double> parameters);
  double sse_sir_multi(vector<double> parameters);
  vector<vector<double> > sse_sir3(vector<double> parameters);
  vector<vector<double> > sse_sir2(vector<double> parameters);
  vector<vector<double> > combine_vectors(vector<vector<double> > data1, vector<vector<double> > data2);
  double add_arrays(vector<vector<double> > data1, vector<vector<double> > data2);
  void update_params(double parameters[7]);
  void user_params();
  void fit_model();
  void update_data(vector<vector<double> > x);
  vector<double> rand_params3();
  vector<double> rand_params4();

};
#endif
