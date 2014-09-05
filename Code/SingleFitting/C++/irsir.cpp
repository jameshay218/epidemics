using namespace std;

#include "irsir.hpp"


IRSIR::~IRSIR(){
}

/* ============================== MODEL EQUATIONS ============================== */

void IRSIR::Diff(vector<double> Pop) {
  // The differential equations
  dPop[0] = - beta*Pop[0]*Pop[1];              // dS/dt
  dPop[1] = beta*Pop[0]*Pop[1] - gamma*Pop[1]*Pop[2];   // dI/dt
  dPop[2] = gamma*Pop[1]*Pop[2];                    // dR/dt
}



/* ================================== SSE FITTING PROCEDURE ================================= */

vector<vector<double> > IRSIR::ode_solve(vector<double> parameters){
  pars = parameters;
  reset_models((current_data.size()/step));
  beta = exp(parameters[0]);
  gamma = exp(parameters[1]);
  populations[0] = exp(parameters[2]);

  if(optimT0 && parameters.size() > 2) 
    if(optimI0) t0 = exp(parameters[4]);
    else t0 = exp(parameters[3]);
  else t0 = 0;
  if(optimI0 && parameters.size() > 2) populations[1] = exp(parameters[3]);
  else populations[1] = 1.0;
  
  populations[2] = 1.0;
  Solve_Eq_t0(temp_model, 1);
  return(temp_model);
}

