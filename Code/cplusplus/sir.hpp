#ifndef SIR_HPP
#define SIR_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

using namespace std;

class SIR{
private:
  double t, step, t0, beta, gamma, tmax;
  double S,I,R,Pop[3];
  double dPop[3];
  vector<vector<double> > current_data;

public:
  SIR();
  SIR(double _tmax, vector<vector<double> > x);
  ~SIR();

  void Diff(double Pop[3]);
  void Runge_Kutta();
  void Solve_Eq(vector<vector<double> >& data);
  void Solve_Eq_t0(vector<vector<double> >& _results);
  void Solve_Eq_total(vector<vector<double> >& _results);
  void Solve_Eq_total2(vector<vector<int> >& _results);


  double sse_sir(vector<double> parameters);
  double sse_sir_t0(vector<double> parameters);
  double sse_sir_multi(vector<double> parameters);
  
  vector<vector<double> > combine_vectors(vector<vector<double> > data1, vector<vector<double> > data2);
  double calculate_SSE(vector<vector<double> > data1, vector<vector<double> > data2);
  vector<vector<int> > combined_model(vector<double> parameters);

  void update_params(double parameters[7]);
  void user_params();
  void fit_model();
  void update_data(vector<vector<double> > x);
  vector<double> rand_params3();
  vector<double> rand_params4();
  bool param_check();

  vector<vector<double> > sse_sir_combined(vector<double> parameters);
  vector<vector<double> > sse_sir_single(vector<double> parameters);
  vector<vector<vector<double> > > sse_sir_components(vector<double> parameters);

  double dpois(vector<vector<double> > model, vector<vector<double> > data);
  double mle_sir(vector<double> parameters); 
};
#endif
