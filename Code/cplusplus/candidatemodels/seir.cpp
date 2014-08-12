using namespace std;

#include "seir.hpp"

SEIR::~SEIR(){
}

/* ============================== MODEL EQUATIONS ============================== */

void SEIR::Diff(vector<double> Pop) {
  // The differential equations
  dPop[0] = - beta*Pop[0]*Pop[2];              // dS/dt
  dPop[1] = beta*Pop[0]*Pop[2] - alpha*Pop[1];   // dI/dt
  dPop[2] = alpha*Pop[1] - gamma*Pop[2];      //dE/dt
  dPop[3] = gamma*Pop[2];                    // dR/dt
}


/* ============================== MLE ================================== */

/* Returns the negative log likelihood of a given set of parameters given
   the current set of data */
double SEIR::mle_sir(vector<double> parameters) {
  pars = parameters;
  reset_models(current_data.size());
  // For each set of 4 parameters in the passed vector, solve a set of ODEs and save this result
  beta = exp(parameters[0]);
  alpha = exp(parameters[1]);
  gamma = exp(parameters[2]);
  populations[0] = exp(parameters[3]);
  t0 = exp(parameters[4]);
  
  populations[1] = 0.0;
  populations[2] = 1.0;
  populations[3] = 0.0;
  
  Solve_Eq_t0(temp_model, 2);
  total_model = combine_vectors(total_model, temp_model);

  // Calculate the overall SSE between the combined model and current data
  double mle = dpois(total_model,current_data);
  return(mle);
}




/* ================================== SSE FITTING PROCEDURE ================================= */


/* Calculates a set of ODEs for each set of 4 parameters passed and returns the combined model */
vector<vector<double> > SEIR::ode_solve_combined(vector<double> parameters){
  pars = parameters;
  reset_models(tmax);
  beta = exp(parameters[0]);
  alpha = exp(parameters[1]);
  gamma = exp(parameters[2]);
  populations[0] = exp(parameters[3]);
  t0 = exp(parameters[4]);
    
  populations[1] = 0.0;
  populations[2] = 1.0;
  populations[3] = 0.0;
  cout << beta << ' ' << alpha << ' ' << gamma << ' ' << populations[0] << ' ' << t0 << endl;      
  Solve_Eq_total(temp_model, 2);
  total_model = combine_vectors(total_model, temp_model);

  return(total_model);
}

vector<vector<double> > SEIR::ode_solve(vector<double> parameters){
  pars = parameters;
  reset_models(current_data.size());
  
  beta = exp(parameters[0]);
  alpha = exp(parameters[1]);
  gamma = exp(parameters[2]);
  populations[0] = exp(parameters[3]);
  t0 = exp(parameters[4]);
    
  cout << beta << ' ' << alpha << ' ' << gamma << ' ' << populations[0] << ' ' << t0 << endl;
  populations[1] = 0.0;
  populations[2] = 1.0;
  populations[3] = 0.0;
  Solve_Eq_t0(temp_model,2);
  return(temp_model);
}

bool SEIR::param_check(){
  if(beta < 0.0005 || beta > 0.05){
    return false;
  }
  if(alpha < 0.0005 || alpha > 0.5){
    return false;
  }
  if(gamma < 0.001 || gamma > 0.5){
    return false;
  }
  if(populations[0] < 100 || populations[0] > 2000){
    return false;
  }
  if(t0 < 1 || t0 > 100){
    return false;
  }
  return true;
}

