using namespace std;

#include "seir.hpp"

SEIR::~SEIR(){
}

/* ============================== MODEL EQUATIONS ============================== */

void SIR::Diff(vector<double> Pop) {
  // The differential equations
  dPop[0] = - beta*Pop[0]*Pop[2];              // dS/dt
  dPop[1] = beta*Pop[0]*Pop[2] - alpha*Pop[1];   // dI/dt
  dPop[2] = alpha*Pop[1] - gamma*Pop[2]
  dPop[3] = gamma*Pop[2];                    // dR/dt
}


/* ============================== MLE ================================== */

/* Returns the negative log likelihood of a given set of parameters given
   the current set of data */
double SIR::mle_sir(vector<double> parameters) {
  pars = parameters;
  reset_models(current_data.size());
  // For each set of 4 parameters in the passed vector, solve a set of ODEs and save this result
  for(unsigned int k = 0; k<(parameters.size()/4);++k){
    beta = exp(parameters[4*k]);
    alpha = exp(parameters[1+(4*k)]);
    gamma = exp(parameters[2+(4*k)]);
    populations[0] = exp(parameters[3+(4*k)]);
    t0 = exp(parameters[4+(4*k)]);
    
    populations[1] = 0.0;
    populations[2] = 1.0;
    populations[3] = 0.0;
	       
    Solve_Eq_t0(temp_model);
    total_model = combine_vectors(total_model, temp_model);
  }
  // Calculate the overall SSE between the combined model and current data
  double mle = dpois(total_model,current_data);
  return(mle);
}




/* ================================== SSE FITTING PROCEDURE ================================= */

/* Calculates the overall SSE from the current data set given a vector of parameters */
double SIR::overall_sse(vector<double> parameters){  
  pars = parameters;
  reset_models(current_data.size());
  
  // For each set of 4 parameters in the passed vector, solve a set of ODEs and save this result
  for(unsigned int k = 0; k<(parameters.size()/4);++k){
    beta = exp(parameters[4*k]);
    alpha = exp(parameters[1+(4*k)]);
    gamma = exp(parameters[2+(4*k)]);
    populations[0] = exp(parameters[3+(4*k)]);
    t0 = exp(parameters[4+(4*k)]);
    
    populations[1] = 0.0;
    populations[2] = 1.0;
    populations[3] = 0.0;
	       
    if(!param_check()){
      return(99999999999.9);
    }
    Solve_Eq_t0(temp_model);
    total_model = combine_vectors(total_model, temp_model);
  }
  // Calculate the overall SSE between the combined model and current data
  double sse = calculate_SSE(current_data,total_model, 2);
  return(sse);
}

/* Calculates a set of ODEs for each set of 4 parameters passed and returns the combined model */
vector<vector<double> > SIR::ode_solve_combined(vector<double> parameters){
  pars = parameters;
  reset_models(tmax);
  
  for(unsigned int k = 0; k<(parameters.size()/4);++k){
 beta = exp(parameters[4*k]);
    alpha = exp(parameters[1+(4*k)]);
    gamma = exp(parameters[2+(4*k)]);
    populations[0] = exp(parameters[3+(4*k)]);
    t0 = exp(parameters[4+(4*k)]);
    
    populations[1] = 0.0;
    populations[2] = 1.0;
    populations[3] = 0.0;
	       
    Solve_Eq_total(temp_model);
    total_model = combine_vectors(total_model, temp_model);
  }
  return(total_model);
}

/* As above, calculates a set of ODEs for each set of 4 parameters passed, but returns a vector of
   each sub epidemic */
vector<vector<vector<double> > > SIR::ode_solve_separate(vector<double> parameters){
  pars = parameters;
  reset_models(tmax);
  vector<vector<vector<double> > > models;  

  for(unsigned int k = 0; k<(parameters.size()/4);++k){
    beta = exp(parameters[4*k]);
    alpha = exp(parameters[1+(4*k)]);
    gamma = exp(parameters[2+(4*k)]);
    populations[0] = exp(parameters[3+(4*k)]);
    t0 = exp(parameters[4+(4*k)]);
    
    populations[1] = 0.0;
    populations[2] = 1.0;
    populations[3] = 0.0;
	       

    Solve_Eq_total(temp_model);
    models.push_back(temp_model);
  }
  return(models);
}


bool SIR::param_check(){
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

