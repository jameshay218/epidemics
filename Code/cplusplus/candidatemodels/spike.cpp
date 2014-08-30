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

vector<vector<double> > Spike::ode_solve(vector<double> parameters){
  pars = parameters;
  reset_models(current_data.size()/step);
  gamma = parameters[0];
  /*if(optimT0 && parameters.size() > 2) t0 = exp(parameters[2]);
    else t0 = 0;*/
  populations[0] = parameters[1];
  //t0 = parameters[2];
  t0 = detectionTime;
  //t0 = 25;
  Solve_Eq_t0(temp_model,0);
  return(temp_model);
}

/* Calculates a set of ODEs for each set of 4 parameters passed and returns the combined model */
vector<vector<double> > Spike::ode_solve_combined(vector<double> parameters){
  pars = parameters;
  reset_models((tmax/step));
  gamma = parameters[0];
  populations[0] = parameters[1];
  //t0 = parameters[2];
  t0 = detectionTime;
  Solve_Eq_total(temp_model, 0);
  return(temp_model);
 }
