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
  double t, step, t0, tmax, infecteds;
  vector<double> dPop, dPop1, dPop2, dPop3, dPop4, tmpPop, initialPop, populations;
  vector<double> pars;
  vector<vector<double> > current_data, temp_model, total_model;
  int noPops, diffIndex, detectionTime; //resize the pop vectors to fit this
  bool optimI0, optimT0;
  EpiType type;
  
public:
  int infectedIndex;
  Epidemic();
  Epidemic(double _tmax, vector<vector<double> > x, EpiType _type, int detection, double _infected);
  virtual ~Epidemic() = 0;

  vector<double> return_parameters();
  EpiType return_type() { return type;};
  int return_detection_time(){return detectionTime;};
  double return_infecteds(){return infecteds;};
  void update_data(vector<vector<double> > x);
  void reset_models(int size);

  /* NORMAL ODE SOLVER */
  virtual void Diff(vector<double> Pop) = 0;
  void Runge_Kutta();
  void Solve_Eq_t0(vector<vector<double> >& _results, int index);
  void Solve_Eq_total(vector<vector<double> >& _results, int index);
  

  /* SSE SOLVER */
  virtual vector<vector<double> > ode_solve(vector<double> parameters) = 0;
 
};
#endif
