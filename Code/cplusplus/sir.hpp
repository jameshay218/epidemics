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

  // Housekeeping functions
  SIR();
  SIR(double _tmax, vector<vector<double> > x);
  ~SIR();
  void update_params(double parameters[7]);
  vector<double> rand_params4();
  vector<vector<double> > combine_vectors(vector<vector<double> > data1, vector<vector<double> > data2);void update_data(vector<vector<double> > x);
  bool param_check();



  /* NORMAL ODE SOLVER */
  void Diff(double Pop[3]);
  void Runge_Kutta();
  void Solve_Eq_t0(vector<vector<double> >& _results);
  void Solve_Eq_total(vector<vector<double> >& _results);
  


  /* MLE SOLVER */
  double dpois(vector<vector<double> > model, vector<vector<double> > data);
  double mle_sir(vector<double> parameters); 
  double poisson_pmf(const double k, const double lambda);
  


  /* SSE SOLVER */
  double overall_sse(vector<double> parameters);
  double calculate_SSE(vector<vector<double> > data1, vector<vector<double> > data2);
  vector<vector<double> > combined_model(vector<double> parameters);
  vector<vector<double> > ode_solve_combined(vector<double> parameters);
  vector<vector<double> > solve_single(vector<double> parameters);
  vector<vector<vector<double> > > ode_solve_separate(vector<double> parameters);
  
  

  /* Obsolete functions */
  double sse_sir(vector<double> parameters);
  double sse_sir_t0(vector<double> parameters);
  void Solve_Eq(vector<vector<double> >& data);
  void user_params();
  vector<double> rand_params3();

};
#endif
