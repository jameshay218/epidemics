using namespace std;

#include "serir.hpp"

SERIR::~SERIR(){
}

/* ============================== MODEL EQUATIONS ============================== */

void SERIR::Diff(vector<double> Pop) {
  // The differential equations
  dPop[0] = - beta*Pop[0]*Pop[2];              // dS/dt
  dPop[1] = beta*Pop[0]*Pop[2] - alpha*Pop[1] - mu*Pop[1];   // dI/dt
  dPop[2] = alpha*Pop[1] - gamma*Pop[2];      //dE/dt
  dPop[3] = gamma*Pop[2] + mu*Pop[1];                    // dR/dt
}

/* ================================== SSE FITTING PROCEDURE ================================= */


/* Calculates a set of ODEs for each set of 4 parameters passed and returns the combined model */
vector<vector<double> > SERIR::ode_solve_combined(vector<double> parameters){
  pars = parameters;
  reset_models(tmax/step);
  beta = exp(parameters[0]);
  alpha = exp(parameters[1]);
  mu = exp(parameters[2]);
  gamma = exp(parameters[3]);
  populations[0] = exp(parameters[4]);
  t0 = exp(parameters[5]);
    
  populations[1] = 0.0;
  populations[2] = 1.0;
  populations[3] = 0.0;
  Solve_Eq_total(temp_model, 2);
  return(temp_model);
}

vector<vector<double> > SERIR::ode_solve(vector<double> parameters){
  pars = parameters;
  reset_models(current_data.size()/step);
  beta = exp(parameters[0]);
  alpha = exp(parameters[1]);
  mu = exp(parameters[2]);
  gamma = exp(parameters[3]);
  populations[0] = exp(parameters[4]);
  t0 = exp(parameters[5]);
  
  populations[1] = 0.0;
  populations[2] = 1.0;
  populations[3] = 0.0;
  Solve_Eq_t0(temp_model,2);
  return(temp_model);
}
