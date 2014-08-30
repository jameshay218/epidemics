#include "spike.hpp"

using namespace std;


Spike::~Spike(){
}


/* ============================== MODEL EQUATIONS ============================== */

void Spike::Diff(vector<double> Pop) {
  // The differential equations
  dPop[0] = -gamma*Pop[0];                    // dI/dt
}




/* ================================== SSE FITTING PROCEDURE ================================= */

/* Calculates a set of ODEs for each set of 4 parameters passed and returns the combined model */
vector<vector<double> > Spike::ode_solve(vector<double> parameters){
  pars = parameters;
  reset_models((tmax/step));
  gamma = exp(parameters[0]);
  populations[0] = exp(parameters[1]);
  if(optimT0 && parameters.size() > 2) 
    if(optimI0) t0 = exp(parameters[2]);
    else t0 = exp(parameters[2]);
  else t0 = 1;
  Solve_Eq_total(temp_model, 0);
  return(temp_model);
 }

