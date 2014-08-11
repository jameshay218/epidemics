#ifndef EPIDEMIC_HPP
#define EPIDEMIC_HPP

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

enum EpiType{none, sir, irsir, seir, serir, spike};

class Epidemic{
protected:
  double t, step, t0, tmax;
  vector<double> dPop, dPop1, dPop2, dPop3, dPop4, tmpPop, initialPop, populations;
  vector<double> pars;
  vector<vector<double> > current_data, temp_model, total_model;
  int noPops, diffIndex; //resize the pop vectors to fit this
  EpiType type;

public:
  int infectedIndex;
  Epidemic();
  Epidemic(double _tmax, vector<vector<double> > x, EpiType _type);
  virtual ~Epidemic() = 0;

  vector<vector<double> > combine_vectors(vector<vector<double> > data1, vector<vector<double> > data2);
  vector<double> return_parameters();
  EpiType return_type() { return type;};
  void update_data(vector<vector<double> > x);
  virtual bool param_check() = 0;
  void reset_models(int size);

  /* NORMAL ODE SOLVER */
  virtual void Diff(vector<double> Pop) = 0;
  void Runge_Kutta();
  void Solve_Eq_t0(vector<vector<double> >& _results);
  void Solve_Eq_total(vector<vector<double> >& _results);
  

  /* MLE SOLVER */
  double dpois(vector<vector<double> > model, vector<vector<double> > data);
  virtual double mle_sir(vector<double> parameters) = 0; 
  double poisson_pmf(const double k, const double lambda);
  

  /* SSE SOLVER */
  virtual double overall_sse(vector<double> parameters) = 0;
  virtual vector<vector<double> > ode_solve(vector<double> parameters) = 0;
  virtual vector<vector<double> > ode_solve_combined(vector<double> parameters) = 0; 
  virtual vector<vector<vector<double> > > ode_solve_separate(vector<double> parameters) = 0;
  static double calculate_SSE(vector<vector<double> > data1, vector<vector<double> > data2, int index);
};
#endif
