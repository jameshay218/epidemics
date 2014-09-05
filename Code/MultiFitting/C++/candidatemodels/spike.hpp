#ifndef SPIKE_HPP
#define SPIKE_HPP

#include "epidemic.hpp"

using namespace std;

class Spike: public Epidemic{
protected:
  double gamma;
  
public:
  // Housekeeping functions
  Spike(double _tmax, vector<vector<double> > x, EpiType _type, int detection) : Epidemic(_tmax, x, _type, detection){
    // Stores beta and gamma
    gamma = 0.1;
    infectedIndex = 0;

    pars.clear();
    pars.push_back(gamma);
    
    // Initial population sizes
    populations[0] = 500.0;
    pars.push_back(populations[0]);
    
    pars.push_back(1);
    
    minTime = detection - 2;
    seedTime = detection-1;
    detectionTime = detection;

  };
  ~Spike();    

  /* NORMAL ODE SOLVER */
  virtual void Diff(vector<double> Pop);

  /* SSE SOLVER */
  virtual vector<vector<double> > ode_solve(vector<double> parameters);
  virtual vector<vector<double> > ode_solve_combined(vector<double> parameters);
};
#endif
