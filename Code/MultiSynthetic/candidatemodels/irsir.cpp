using namespace std;

#include "irsir.hpp"


IRSIR::~IRSIR(){
  cout << "deleting" << endl;
}

/* ============================== MODEL EQUATIONS ============================== */

void IRSIR::Diff(vector<double> Pop) {
  // The differential equations
  dPop[0] = - beta*Pop[0]*Pop[1];              // dS/dt
  dPop[1] = beta*Pop[0]*Pop[1] - gamma*Pop[1]*Pop[2];   // dI/dt
  dPop[2] = gamma*Pop[1]*Pop[2];                    // dR/dt
}



/* ================================== SSE FITTING PROCEDURE ================================= */

/* Calculates a set of ODEs for each set of 4 parameters passed and returns the combined model */
vector<vector<double> > IRSIR::ode_solve_combined(vector<double> parameters){
  pars = parameters;
  reset_models((tmax/step));
  beta = parameters[0];
  gamma = parameters[1];
  populations[0] = parameters[2];
  if(active) t0 = parameters[3];
  else t0 = detectionTime;
  populations[1] = 1.0;
  populations[2] = 0.0;
  Solve_Eq_total(temp_model, 1);
  return(temp_model);
 }

vector<vector<double> > IRSIR::ode_solve(vector<double> parameters){
  pars = parameters;
  reset_models((current_data.size()/step));
  beta = parameters[0];
  gamma = parameters[1];
  populations[0] = parameters[2];
  if(active) t0 = parameters[3];
  else t0 = detectionTime;
  populations[1] = 1.0;
  populations[2] = 0.0;
  Solve_Eq_t0(temp_model, 1);
  return(temp_model);
}

