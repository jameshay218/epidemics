using namespace std;

#include "sir.hpp"


SIR::~SIR(){

}

/* ============================== MODEL EQUATIONS ============================== */

void SIR::Diff(vector<double> Pop) {
  // The differential equations
  dPop[0] = - beta*Pop[0]*Pop[1];              // dS/dt
  dPop[1] = beta*Pop[0]*Pop[1] - gamma*Pop[1];   // dI/dt
  dPop[2] = gamma*Pop[1];                    // dR/dt
}



/* ================================== SSE FITTING PROCEDURE ================================= */

/* Calculates a set of ODEs for each set of 4 parameters passed and returns the combined model */
vector<vector<double> > SIR::ode_solve_combined(vector<double> parameters){
  pars = parameters;
  reset_models((tmax/step));
  beta = exp(parameters[0]);
  gamma = exp(parameters[1]);
  populations[0] = exp(parameters[2]);
  t0 = exp(parameters[3]);
  populations[1] = 1.0;
  populations[2] = 0.0;
  Solve_Eq_total(temp_model, 1);
  return(temp_model);
 }

vector<vector<double> > SIR::ode_solve(vector<double> parameters){
  pars = parameters;
  reset_models((current_data.size()/step));
  beta = exp(parameters[0]);
  gamma = exp(parameters[1]);
  populations[0] = exp(parameters[2]);
  t0 = exp(parameters[3]);
  populations[1] = 1.0;
  populations[2] = 0.0;
  Solve_Eq_t0(temp_model, 1);
  return(temp_model);
}

