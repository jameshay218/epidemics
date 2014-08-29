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
  int noPops, diffIndex, detectionTime, seedTime, minTime; //resize the pop vectors to fit this
  EpiType type;

public:
  int infectedIndex;
  Epidemic();
  Epidemic(double _tmax, vector<vector<double> > x, EpiType _type, int detection);
  virtual ~Epidemic() = 0;

  vector<double> return_parameters();
  EpiType return_type() { return type;};
  int return_detection_time(){return detectionTime;};
  int return_seed_time(){return seedTime;};
  int return_min_time(){return minTime;};
  void update_data(vector<vector<double> > x);
  void reset_models(int size);

  /* NORMAL ODE SOLVER */
  virtual void Diff(vector<double> Pop) = 0;
  void Runge_Kutta();
  void Solve_Eq_t0(vector<vector<double> >& _results, int index);
  void Solve_Eq_total(vector<vector<double> >& _results, int index);
  

  /* SSE SOLVER */
  virtual vector<vector<double> > ode_solve(vector<double> parameters) = 0;
  virtual vector<vector<double> > ode_solve_combined(vector<double> parameters) = 0; 
};
#endif
