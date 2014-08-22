#ifndef IRSIR_HPP
#define IRSIR_HPP

#include "epidemic.hpp"

using namespace std;

class IRSIR: public Epidemic{
protected:
  double beta, gamma;
  
public:
  // Housekeeping functions
  IRSIR(double _tmax, vector<vector<double> > x, EpiType _type, int detection) : Epidemic(_tmax, x, _type, detection){
    // Stores beta and gamma
    beta = 0.001;
    gamma = 0.1;
    infectedIndex = 2;
    pars.clear();
    pars.push_back(beta);
    pars.push_back(gamma);

    // Initial population sizes
    populations[0] = 500.0;
    pars.push_back(populations[0]);
    populations[1] = 1.0;
    populations[2] = 0.0;

    pars.push_back(1);
  };
  ~IRSIR();    


  /* NORMAL ODE SOLVER */
  virtual void Diff(vector<double> Pop);
 

  
  /* SSE SOLVER */
  //double overall_sse(vector<double> parameters);
  virtual vector<vector<double> > ode_solve(vector<double> parameters);
  virtual vector<vector<double> > ode_solve_combined(vector<double> parameters);
  
};
#endif
