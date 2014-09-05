#ifndef SPIKE_HPP
#define SPIKE_HPP

#include "epidemic.hpp"

using namespace std;

class Spike: public Epidemic{
protected:
  double gamma;
  
public:
  // Housekeeping functions
  Spike(double _tmax, vector<vector<double> > x, EpiType _type, int detection, double _infecteds, bool includeT0, bool includeI0) : Epidemic(_tmax, x, _type, detection, _infecteds){
    // Stores beta and gamma
    gamma = 0.1;
    infectedIndex = 0;

    pars.clear();
    pars.push_back(gamma);
    
    // Initial population sizes
    populations[0] = 500.0;
    pars.push_back(populations[0]);
 
    optimT0 = includeT0;
    optimI0 = includeI0;

    if(includeT0) pars.push_back(1);    


  };
  ~Spike();    

  /* NORMAL ODE SOLVER */
  virtual void Diff(vector<double> Pop);

  /* SSE SOLVER */
  virtual vector<vector<double> > ode_solve(vector<double> parameters);

};
#endif
