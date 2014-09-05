#ifndef SEIR_HPP
#define SEIR_HPP

#include "epidemic.hpp"

using namespace std;

class SEIR: public Epidemic{
protected:
  double beta, gamma, alpha;
  
public:
  // Housekeeping functions
  SEIR(double _tmax, vector<vector<double> > x, EpiType _type, int detection) : Epidemic(_tmax, x, _type, detection){
    // Stores beta and gamma
    beta = 0.001;
    alpha = 0.01;
    gamma = 0.1;
    infectedIndex = 3;

    pars.clear();
    pars.push_back(beta);
    pars.push_back(alpha);
    pars.push_back(gamma);
    
    // Initial population sizes
    populations[0] = 500.0;
    pars.push_back(populations[0]);
    populations[1] = 0.0;
    populations[2] = 1.0;
    populations[3] = 0.0;
    pars.push_back(1);
    minTime = detection-40;
    seedTime = detection-5;
    
  };
  ~SEIR();    
 
  /* NORMAL ODE SOLVER */
  virtual void Diff(vector<double> Pop);

  
  /* SSE SOLVER */
  virtual vector<vector<double> > ode_solve(vector<double> parameters);
  virtual vector<vector<double> > ode_solve_combined(vector<double> parameters);

  
};
#endif
