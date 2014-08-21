using namespace std;

#include "sir.hpp"


SIR::~SIR(){
  cout << "deleting" << endl;
}

/* ============================== MODEL EQUATIONS ============================== */

void SIR::Diff(vector<double> Pop) {
  // The differential equations
  dPop[0] = - beta*Pop[0]*Pop[1];              // dS/dt
  dPop[1] = beta*Pop[0]*Pop[1] - gamma*Pop[1];   // dI/dt
  dPop[2] = gamma*Pop[1];                    // dR/dt
}


/* ============================== MLE ================================== */

/* Returns the negative log likelihood of a given set of parameters given
   the current set of data */
double SIR::mle_sir(vector<double> parameters) {
  pars = parameters;
  reset_models(current_data.size());
  // For each set of 4 parameters in the passed vector, solve a set of ODEs and save this result
  beta = exp(parameters[0]);
  gamma = exp(parameters[1]);
  populations[0] = exp(parameters[2]);
  t0 = exp(parameters[3]);
  populations[1] = 1.0;
  populations[2] = 0.0;
  Solve_Eq_t0(temp_model, 1);
  // Calculate the overall SSE between the combined model and current data
  double mle = dpois(temp_model,current_data);
  return(mle);
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


bool SIR::param_check(){
  if(beta < 0.00001 || beta > gamma){
    return false;
  }
  if(gamma < 0.001 || gamma > 0.5){
    return false;
  }
  if(populations[0] < 50 || populations[0] > 5000){
    return false;
  }
  if(t0 < 0 || t0 > 100){
    return false;
  }
  return true;
}

