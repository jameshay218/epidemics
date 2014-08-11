#ifndef SIR_HPP
#define SIR_HPP

#include "epidemic.hpp"

using namespace std;

class SIR: public Epidemic{
protected:
  double beta, gamma;
  
public:
  // Housekeeping functions
  SIR(double _tmax, vector<vector<double> > x, EpiType _type) : Epidemic(_tmax, x, _type){
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
  ~SIR();    
  bool param_check();



  /* NORMAL ODE SOLVER */
  void Diff(vector<double> Pop);
 

  /* MLE SOLVER */
  double mle_sir(vector<double> parameters); 
  
  
  /* SSE SOLVER */
  double overall_sse(vector<double> parameters);
  vector<vector<double> > ode_solve(vector<double> parameters);
  vector<vector<double> > combined_model(vector<double> parameters);
  vector<vector<double> > ode_solve_combined(vector<double> parameters);
  vector<vector<vector<double> > > ode_solve_separate(vector<double> parameters);
  
};
#endif
