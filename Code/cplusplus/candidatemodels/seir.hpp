#ifndef SEIR_HPP
#define SEIR_HPP

#include "epidemic.hpp"

using namespace std;

class SEIR: public Epidemic{
protected:
  double beta, gamma;
  
public:
  // Housekeeping functions
  SEIR(double _tmax, vector<vector<double> > x, EpiType _type) : Epidemic(_tmax, x, _type){
    // Stores beta and gamma
    beta = 0.001;
    alpha = 0.01;
    gamma = 0.1;
    
    // Initial population sizes
    populations[0] = 500.0;
    populations[1] = 0.0;
    populations[2] = 1.0;
    populations[3] = 0.0;
  };
  ~SEIR();    
  bool param_check();



  /* NORMAL ODE SOLVER */
  void Diff(vector<double> Pop);

  /* MLE SOLVER */
  double mle_sir(vector<double> parameters); 
  
  
  /* SSE SOLVER */
  double overall_sse(vector<double> parameters);
  vector<vector<double> > combined_model(vector<double> parameters);
  vector<vector<double> > ode_solve_combined(vector<double> parameters);
  vector<vector<vector<double> > > ode_solve_separate(vector<double> parameters);
  
};
#endif
