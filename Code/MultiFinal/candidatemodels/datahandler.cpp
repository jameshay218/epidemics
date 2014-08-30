using namespace std;

#include "datahandler.hpp"
#include "epidemic.hpp"
#include "sir.hpp"
#include "seir.hpp"
#include "spike.hpp"
#include "serir.hpp"
#include "irsir.hpp"
#include "simplex.hpp"
#include <string>

const float EPSILON = 0.1;
const double upperS0 = 5000;
const double lowerS0 = 100;
const double seedT0 =  200;
const double lowerT0 = -200;



/* ================================ FINAL LEAST SQUARES ============================= */

void Handler::realtime_fit_multi(double targetRSq){
  clock_t t1, t2,t3; // To record the total run time
  EpiType newEpidemic; // Keep track of the type of a new epidemic to be added
  Epidemic *toRemove, *bestFit; // Keep track of a newly created epidemic to allow deletion
  vector<EpiType> candidateModels,selectedTypes;
  vector<double> finalParams, finalParamsTemp; // Keep track of parameters from optimisation
  vector<vector<double> > combinedResults, residuals, combinedResultsTemp,saveResults; // 
  vector<vector<vector<double> > > componentResults, componentResultsTemp;
  double RSquare, tempRSquare, SSE, tempSSE;
  int iterations, tempIterations;
  double bestAICc, tempAICc, AICc,RMSE;

  //  candidateModels.push_back(spike);
  candidateModels.push_back(sir);
  //candidateModels.push_back(seir);
 //candidateModels.push_back(serir);
 //candidateModels.push_back(irsir);

  // Start of model fitting process
  t1=clock();
  srand(time(NULL));
  
  // Start off with a baseline model of the mean of the first 4 points
  for(unsigned int j = 0; j <= 4; ++j){
    this->temp_data.push_back(this->current_data[j]);
  }
  this->baseModel = this->current_model = base_model(this->temp_data);

  epidemics.push_back(new_epidemic(sir,10.0));
  epidemicSizes.push_back(4);

   // For each time point in the current data set, carry out the fitting procedure
  for(unsigned int i = 30;i<this->current_data.size();++i){
    cout << endl << "-----------------" << endl;
    cout << "Iteration number " << i << endl;
    
    // Clear temp vector and temp epidemics
    this->temp_data.clear();
    this->tempEpidemics.clear();
    this->tempSizes.clear();

    toRemove = NULL;

    // Store all data up to the current index and update epidemics with data
    for(unsigned int j = 0; j < i; ++j){
      this->temp_data.push_back(this->current_data[j]);
    } 
    for(unsigned int f = 0;f<this->epidemics.size();++f){
      this->epidemics[f]->update_data(this->temp_data);
    }
    this->empty_model = create_empty_data_vector(temp_data.size());

    // Create an optimised model fit for the given data (temp) and calculate RSquare from the SSE
    SSE = optimiseEpidemics(finalParams, combinedResults, componentResults, iterations);
    RSquare = 1 - SSE/SStot(temp_data, 1);
    cout << "Optimised fit, RSquare: " << RSquare << endl;
    AICc = aicc(SSE, temp_data.size(), finalParams.size());
    // If fit is sufficiently good with k epidemics, try fitting with k-1 epidemics. If this fit
    // is sufficient, use k-1 epidemics.
    if(RSquare > targetRSq && this->epidemics.size() > 0){
      cout << "Considering removal of an epidemic" << endl;
      // Store current set of epidemics
      this->tempEpidemics = this->epidemics;
      this->tempSizes = this->epidemicSizes;

      // Try removing each epidemic so far
      toRemove = NULL;
      for(unsigned int j = 0;j<this->epidemics.size();++j){
	// Create vector of epidemic types that have been removed
	vector<EpiType> removedEpidemics;

	// Try removing the j'th epidemic. If that type has already been removed, increase j.
	// if epidemics[j]->return_type() is not in removedEpidemics, then that. Otherwise, j++ 
	// and try again. This ensures that only one epidemic of each type is removed at a time.
	if(!check_already_tested(removedEpidemics, this->tempEpidemics[j]->return_type())){
	  this->epidemics = fewer_epidemics(j);
	  this->epidemicSizes = fewer_sizes(j);
	  	
	  // Calculate the fit without this epidemic
	  tempSSE = optimiseEpidemics(finalParamsTemp, combinedResultsTemp, componentResultsTemp, tempIterations);
	  tempRSquare = 1 - tempSSE/SStot(temp_data,1);
	  tempAICc = aicc(tempSSE, temp_data.size(), finalParamsTemp.size());
	  // If sufficient fit from removing an epidemic, permanently remove epidemic from set
	  // Otherwise, restore epidemics to full number. Also track which removal was best
	  if(tempAICc < AICc){
	    if(toRemove != NULL) delete toRemove;
	    toRemove = this->tempEpidemics[j];
	    cout << "Fit sufficient. Removing an epidemic" << endl;
	    SSE = tempSSE;
	    AICc = tempAICc;
	    RSquare = tempRSquare;
	    finalParams = finalParamsTemp;
	    combinedResults = combinedResultsTemp;
	    componentResults = componentResultsTemp;
	    iterations = tempIterations;
	  }
	  else{
	    this->epidemics = this->tempEpidemics;
	    this->epidemicSizes = this->tempSizes;
	  }
	}
      }
      if(toRemove != NULL) remove_epidemic(toRemove);
    }

    // Reset pointers in case epidemic was deleted (check for memory leak)
    reset_epidemics();
    
    // Check if a new epidemic has started. If so, and RSquare from current K epidemics 
    // insufficient, try fitting another. If this is a significantly better fit, permanently
    // add an extra epidemic.
    residuals = get_residuals(this->temp_data, this->current_model, 1);
    newEpidemic = check_epidemic(residuals);

    // If so, try adding epidemic
    if(newEpidemic != none && RSquare < targetRSq){
      cout << "Epidemic detected!" << endl;
      bestFit = NULL;
      bestAICc = 99999999999999.9;
      // Temporarily store the current epidemic state and add a new epidemic to the current set
      this->tempEpidemics = this->epidemics;
      this->tempSizes = this->epidemicSizes;

      for(unsigned int f = 0; f < candidateModels.size(); ++f){
	toRemove = NULL;
	newEpidemic = candidateModels[f];
	cout << "Considering addition of ";
	print_epidemic_type(newEpidemic);
	cout << endl;
	toRemove = new_epidemic(newEpidemic, i);
	this->epidemics.push_back(toRemove);
	this->epidemicSizes.push_back(toRemove->return_parameters().size());
	
	// Produce optimised fit with added epidemic
	tempSSE = optimiseEpidemics(finalParamsTemp, combinedResultsTemp, componentResultsTemp, tempIterations);
	tempRSquare = 1 - tempSSE/SStot(this->temp_data,1);
	tempAICc = aicc(tempSSE, temp_data.size(), toRemove->return_parameters().size());
	cout << "Fit with additional epidemic: " << tempRSquare << endl;
	cout << "SSE with additional epidemic: " << tempSSE << endl;
	cout << "AICc: " << tempAICc << endl;
	cout << "Iterations taken: " << tempIterations << endl;
	// If fit improved, keep track of epidemic. If best fit, add to epidemics at the end
	if(tempRSquare > RSquare && tempAICc < bestAICc){
	  cout << "AICc improved with additional epidemic. Add epidemic" << endl;
	  if(bestFit != NULL) delete bestFit;
	  bestFit = new_epidemic(newEpidemic,i);
	  SSE = tempSSE;
	  bestAICc = tempAICc;
	  RSquare = tempRSquare;
	  finalParams = finalParamsTemp;
	  combinedResults = combinedResultsTemp;
	  componentResults = componentResultsTemp;
	  iterations = tempIterations;
	}
	remove_epidemic(toRemove);
	reset_epidemics();
	this->epidemicSizes = this->tempSizes;
      }
      if(bestFit != NULL){
	  cout << "Adding epidemic type: ";
	  print_epidemic_type(bestFit->return_type());
	  selectedTypes.push_back(bestFit->return_type());
	  cout << endl;
	  epidemics.push_back(bestFit);
	  this->epidemicSizes.push_back(bestFit->return_parameters().size());
	  bestFit = NULL;
	}
      }
    // Reset pointers in case epidemic was deleted;
    reset_epidemics();
       
    // Transform parameters back to normal space
    finalParams = convert_all_params_back(finalParams);
      
    plotGraphMulti(componentResults, combinedResults, this->temp_data, i, finalParams, RSquare, 2);    
    
    // Print out results for this iteration
    cout << "Final SSE was: " << SSE << endl;
    cout << "Final RSquare: " << RSquare << endl;
    cout << "Final parameters: ";
    printcon(finalParams) ;
    finalParams.push_back(RSquare);
    finalParams.push_back(RMSE);
    saveResults.push_back(finalParams);
    RMSE = sqrt((SSE/temp_data.size()));
    cout << "RMSE: " << RMSE << endl;
    cout << endl << "-----------------" << endl;
    this->current_model = combinedResults;
  }
  t2=clock();
  t3 = (t2-t1)/CLOCKS_PER_SEC;
  finalParams.clear();
  finalParams.push_back(double(t3));
  saveResults.push_back(finalParams);
  print_vector(saveResults);

  ofstream myfile; 
  myfile.open ((saveLocation + "/results.csv")); 
  
  for(unsigned int i = 0; i < saveResults.size();++i){
    for(unsigned int j = 0;j < saveResults[i].size();++j){
      myfile << saveResults[i][j] << ",";
    }
    myfile << "\n";
  }
  for(unsigned int i = 0; i < selectedTypes.size();++i){
    myfile << selectedTypes[i] << ",";
  }
  myfile << "\n";
  myfile.close();
  cout << "Time taken: " << (t2-t1)/CLOCKS_PER_SEC << endl << endl;
}


double Handler::optimiseEpidemics(vector<double> &parameters, vector<vector<double> > &results, vector<vector<vector<double> > > &allResults, int& itr){
  double SSE = 99999999999999999.9;
  double tempSSE;
  vector<double> tempPar, seedParams;
  Simplex simplex;
  int iterations = 10001;
  // If no epidemics have yet been detected, use the mean of the current data as the 
  // current model and return the corresponding SSE.
  if(epidemics.size() == 0){
    results = base_model(temp_data);
    SSE = calculate_SSE(results, temp_data);
    parameters.clear();
    allResults.clear();
    //allResults.push_back(results);
    allResults.push_back(baseModel);
    return(SSE);
  }
  srand(time(NULL));
  parameters = generate_seed_parameters();
  // If there are epidemics to be fitted, perform 40 random fits and keep best fitting model
  for(int index=0;index<10;index++){
    
    
    // Create a list of random seed parameters
    seedParams.clear();      
    seedParams = generate_seed_parameters();
    //printcon(convert_all_params_back(seedParams));
    
    while(fitEpidemics(seedParams) != fitEpidemics(seedParams)){
      //cout << "Changing seed" << endl;
      seedParams.clear();            
      seedParams = generate_seed_parameters();
    }
    
    // Get the optimised parameters from nelder mead algorithm
    tempPar = simplex.neldermead(&Handler::fitEpidemics, *this,  seedParams, iterations);
    // Store the SSE value for this
    tempSSE = fitEpidemics(tempPar);
	
    // If this SSE value is better than the previous, store it and the
    // corresponding parameters
    if(tempSSE < SSE){
      itr = iterations;
      SSE = tempSSE;
      parameters=tempPar;
    }
    cout << "." << flush;
  }
  /* TEMPORARY STUFF */
  /*tempPar.clear();
  tempPar.push_back(log(0.001));
  tempPar.push_back(log(0.1));
  tempPar.push_back(log(500));
  tempPar.push_back(log(10));
  itr = 100;
  SSE = fitEpidemics(tempPar);
  parameters = tempPar;
  */

  // Get the combined values from these parameters, as well as a vector of 
  // each sub-epidemic
  results = ode_solve(parameters);
  //  print_vector(results);
  allResults.clear();
  allResults = ode_solve_separate(parameters);
  
  cout << endl;
  return(SSE);

}



/* ========================================= MODEL FITTING PROCEDURES ==============================*/
// Goes through each epidemic in the current list and does the following: 1. Store the correct number
// of parameters from the provided vector; 2. Solve the current epidemic model with those parameters
// 3. Add the contribution of that epidemic to the overall model. 4. Calculate the SSE between the
// current model and current data.  
double Handler::fitEpidemics(vector<double> params){
  int z = 0;
  current_model = empty_model;
  current_model = combine_vectors(current_model,this->baseModel);
  params = convert_all_params_back(params);
  for(unsigned int i = 0;i<epidemics.size();++i){
    tempParams.clear();
    temp_model.clear();
    for(int x = 0;x<this->epidemicSizes[i];++x){
      tempParams.push_back(params[z]);
      z++;
    }
    temp_model = epidemics[i]->ode_solve(tempParams);
    current_model = combine_vectors(current_model, temp_model);
  }
  //cout << calculate_SSE(current_model,temp_data) << endl;
  return(calculate_SSE(current_model, temp_data));
}

double Handler::fitEpidemicsMLE(vector<double> params){
  int z = 0;
  current_model = empty_model;
  current_model = combine_vectors(current_model,this->baseModel);
  params = convert_all_params_back(params);
  for(unsigned int i = 0;i<epidemics.size();++i){
    tempParams.clear();
    temp_model.clear();
    for(unsigned int x = 0;x<epidemics[i]->return_parameters().size();++x){
      tempParams.push_back(params[z]);
      z++;
    }
    //tempParams = convert_parameters_back(tempParams,epidemics[i]->return_parameters().size());
    temp_model = epidemics[i]->ode_solve(tempParams);
    current_model = combine_vectors(current_model, temp_model);
  }
  return(dpois(current_model, temp_data));
}
  
// As above, but returns the overall model.
vector<vector<double> > Handler::ode_solve(vector<double> params){
  vector<vector<double> > overallResults = create_empty_data_vector(current_data.size());
  vector<double> parr = convert_all_params_back(params);
  overallResults = combine_vectors(overallResults, baseModel);
  int z = 0;
  for(unsigned int i = 0;i<epidemics.size();++i){
    tempParams.clear();
    temp_model.clear();
    for(int x = 0;x<epidemicSizes[i];++x){
      tempParams.push_back(parr[z]);
      z++;
    }
    //tempParams = convert_parameters_back(tempParams,epidemics[i]->return_parameters().size());
    temp_model = epidemics[i]->ode_solve_combined(tempParams);
    
    overallResults = combine_vectors(overallResults,temp_model);

  }
  //print_vector(overallResults);
  return(overallResults);
}

// As above, but returns a vector of the component models.
vector<vector<vector<double> > > Handler::ode_solve_separate(vector<double> params){
  vector<vector<vector<double> > > separateResults;
  vector<double> parr = convert_all_params_back(params);
  int x = 0;
  separateResults.push_back(this->baseModel);
  //print_vector(baseModel);
  for(unsigned int i = 0;i<epidemics.size();++i){
    tempParams.clear();
    temp_model.clear();
    for(unsigned int z = 0;z<epidemics[i]->return_parameters().size();++z){
      tempParams.push_back(parr[x]);
      x++;
    }
    //tempParams = convert_parameters_back(tempParams,epidemics[i]->return_parameters().size());
    temp_model = epidemics[i]->ode_solve_combined(tempParams);
    //print_vector(temp_model);
    separateResults.push_back(temp_model);
  }
  return(separateResults);
}

  

/* Takes the residuals as a 2xN vector. Returns true if the latest X residuals are a certain number
   of standard deviations away from the previous total-X residuals */
EpiType Handler::check_epidemic(vector<vector<double> > residuals){
  EpiType detected = none;
  double mean, sd,sirlimit;
  vector<vector<double> > previousResiduals;
  previousResiduals.clear();
 
  if(residuals.size() < 3) return none;

  // Take all but the latest 3 residuals
  for(unsigned int j = 0; j < residuals.size() - 3;++j){
    previousResiduals.push_back(residuals[j]);
  }

  // Calculate the SD and mean of these residuals
  mean = calculate_mean(previousResiduals, 1);
  sd = calculate_sd(previousResiduals, 1);

  sirlimit = mean + 3*sd;
  
 
  if((residuals[residuals.size()-1][1] > sirlimit)) detected=sir;
  return detected;
}

  
/* ========================= Parameter Housekeeping ===============================*/


vector<double> Handler::convert_all_params_back(vector<double> pars){
    
  vector<double> temporary;
  vector<double> temp2;
  int z = 0;
  temp2.clear();
  for(unsigned int i = 0;i<this->epidemicSizes.size();++i){
    temporary.clear();
    for(int f = 0;f<this->epidemicSizes[i];++f){
      temporary.push_back(pars[z]);
      z++;
    }
    temporary = convert_parameters_back(temporary,this->epidemicSizes[i]);
    temp2 = concatenate_vectors(temp2, temporary);
  }
  return(temp2);
}





vector<double> Handler::convert_parameters_back(vector<double> pars, int number){
  vector<double> temporary;
  
  if(number == 3){
    temporary.push_back(exp(pars[0]));

    if(useLogistic) temporary.push_back(logistic(pars[1],lowerS0,upperS0));
    else{ temporary.push_back(exp(pars[1])); temporary.push_back((exp(pars[2])+5.0));}

  }
  else if(number == 4){
    temporary.push_back(exp(pars[0]));
    temporary.push_back(exp(pars[1]));
  
    if(useLogistic) temporary.push_back(logistic(pars[2],lowerS0,upperS0));
    else {temporary.push_back(exp(pars[2]));
      temporary.push_back((exp(pars[3])+5.0));}
  }
  else if(number == 5){
    temporary.push_back(exp(pars[0]));
    temporary.push_back(exp(pars[1]));
    temporary.push_back(exp(pars[2]));

    if(useLogistic) temporary.push_back(logistic(pars[3],lowerS0,upperS0));
    else {temporary.push_back(exp(pars[3]));
    temporary.push_back((exp(pars[4])+5.0));}
    
  }
  else{
    temporary.push_back(exp(pars[0]));
    temporary.push_back(exp(pars[1]));
    temporary.push_back(exp(pars[2]));
    temporary.push_back((exp(pars[3])));

    if(useLogistic) temporary.push_back(logistic(pars[4],lowerS0,upperS0));
    else {temporary.push_back(exp(pars[4]));
      temporary.push_back(((exp(pars[5])+5.0)));}
  }
  //if(useLogistic) temporary.push_back(logistic((time-seedT0),(time-lowerT0),time));
  if(useLogistic) temporary.push_back(logistic(pars[5],lowerT0,seedT0));
 
  //  printcon(temporary);
  return(temporary);
}

vector<double> Handler::convert_parameters_forward(vector<double> pars, int number, int time){
  vector<double> temporary;
 
  if(number == 3){
    temporary.push_back(log(pars[0]));
    if(useLogistic) temporary.push_back(logit(pars[1],lowerS0,upperS0));
    else temporary.push_back(log(pars[1]));

  }
  else if(number == 4){
    temporary.push_back(log(pars[0]));
    temporary.push_back(log(pars[1]));
  
    if(useLogistic) temporary.push_back(logit(pars[2],lowerS0,upperS0));
    else temporary.push_back(log(pars[2]));
  }
  else if(number == 5){
    temporary.push_back(log(pars[0]));
    temporary.push_back(log(pars[1]));
    temporary.push_back(log(pars[2]));
  
    if(useLogistic) temporary.push_back(logit(pars[3],lowerS0,upperS0));
    else temporary.push_back(log(pars[3]));

  }
  else{
    temporary.push_back(log(pars[0]));
    temporary.push_back(log(pars[1]));
    temporary.push_back(log(pars[2]));
    temporary.push_back(log(pars[3]));
    if(useLogistic) temporary.push_back(logit(pars[4],lowerS0,upperS0));
    else temporary.push_back(log(pars[4]));
  }
  if(useLogistic) temporary.push_back(logit(time,lowerT0,seedT0)); 
  else temporary.push_back(log(time-5.0));
  //  if(useLogistic) temporary.push_back(logit((time-seedT0),(time-lowerT0),time));
  cout << "Time: " << time << endl;
  cout << "Time2: " << logistic(time,lowerT0,seedT0) << endl;
  printcon(temporary);
  printcon(convert_parameters_back(temporary,number));
   
  return(temporary);
}



vector<double> Handler::generate_seed_parameters(){
  vector<double> params;
  vector<double> tempor;
  for(unsigned int i = 0;i<epidemics.size();++i){
    tempor = rand_params(epidemics[i]->return_type());
    tempor = convert_parameters_forward(tempor, epidemics[i]->return_parameters().size(), epidemics[i]->return_detection_time());
    params = concatenate_vectors(params, tempor);
    
  }
  return(params);
}

// Generates a vector of random parameters
// THIS IS WHERE WE GENERATE SEED VALUES
vector<double> Handler::rand_params(EpiType _type){
  vector<double> params;
  double beta,gamma,s0,alpha;
  setprecision(9);
 
  s0 = rand()%2000+500;
  //t0 = rand()%80+1;
  
  //beta = 0.0005;
  //alpha = 0.001;
  //gamma = 0.05;
  //s0 = 1000;
  //  t0 = 15;
  switch(_type){
  case sir:
    beta = (rand()%500+50)/100000.0;
    gamma = (rand()%200+beta)/1000.0;
    params.push_back((beta));
    params.push_back((gamma));
    params.push_back((s0));
    break;
    //params.push_back(log(t0));
    //return(params);
  case seir:
    beta = (rand()%100+1)/10000.0;
    gamma = (rand()%200+1)/1000.0;
    alpha = (rand()%100+1)/10000.0;
    params.push_back((beta));
    params.push_back((alpha));
    params.push_back((gamma));
    params.push_back((s0));
    break;
    //params.push_back(log(t0));
    //return(params);
  case spike:
    gamma = (rand()%200+1)/1000.0;
    params.push_back((gamma));
    params.push_back((s0));
    break;
    //params.push_back(log(t0));
  default:
    params.push_back((gamma));
    params.push_back((s0));
    //params = convert_parameters_forward(params,numb);
   
  }
  return(params);
}


vector<double> Handler::generate_seed_parameters_old(){
  vector<double> params;
  for(unsigned int i = 0;i<epidemics.size();++i){
    params = concatenate_vectors(params,rand_params_old(epidemics[i]->return_type()));
    params.push_back(epidemics[i]->return_detection_time());
  }
  return(params);
}

vector<double> Handler::rand_params_old(EpiType _type){
  vector<double> params;
  double beta,gamma,s0,alpha;
  setprecision(9);
 
  s0 = rand()%2000+500;
 
  switch(_type){
  case sir:
    beta = (rand()%100+1)/10000.0;
    gamma = (rand()%200+beta)/1000.0;
    params.push_back((beta));
    params.push_back((gamma));
    params.push_back((s0));
    break;
  case seir:
    beta = (rand()%100+1)/10000.0;
    gamma = (rand()%200+1)/1000.0;
    alpha = (rand()%100+1)/10000.0;
    params.push_back((beta));
    params.push_back((alpha));
    params.push_back((gamma));
    params.push_back((s0));
    break;
  case spike:
    gamma = (rand()%200+1)/1000.0;
    params.push_back((gamma));
    params.push_back((s0));
    break;
  default:
    params.push_back((gamma));
    params.push_back((s0));
  }
  return(params);
}






/* ================================= MATHEMATICAL FUNCTIONS =====================================*/





/* Calculates and returns the sum of squared errors between to sets of data. Checks 
   the first column to see if on the same time point, and finds squared difference if
   so. If one dataset is at an earlier time point, use the entire squared value as
   the squared residual.
*/
double Handler::calculate_SSE(vector<vector<double> > data1, vector<vector<double> > data2){
  unsigned int i = 0;
  unsigned int j = 0;
  double sse = 0;
  int increment;
  vector<vector<double> > big, small;

  if(data1.size() < data2.size()){
    small = data1;
    big = data2;
    increment = (data2.size()/data1.size());
  }
  else if (data2.size() < data1.size()){
    small = data2;
    big = data1;
    increment = (data1.size()/data2.size());
  }
  else{
    small = data1;
    big = data2;
    increment = 1;
  }

  while(j < small.size()){
    //If indices give same time point, find sse
    if(fabs(small[i][0]- big[j][0])<EPSILON){
      sse += pow((small[i][1] - big[j][1]),2.0);
      i++;
      j+=increment;
    }
    //If first dataset is still before second, full difference
    else if(small[i][0] < big[j][0]){
      sse += pow(small[i][1],2.0);
      i++;
    }
    //...
    else{
      sse += pow(big[j][1],2.0);
      j+=increment;
    }
  }
  return(sse);
}



/* Returns the negative log likelihood of a model given a set of data */
double Handler::dpois(vector<vector<double> > model, vector<vector<double> > data){
  unsigned int i = 0;
  unsigned int j = 0;
  double logLikelihood = 0;
 
  while(j < model.size() && i < data.size()){
    //If indices give same time point, find sse
    if(fabs(model[i][0]- data[j][0])<EPSILON){
      if(model[i][0] == 0 || data[j][0] == 0) logLikelihood += log(1);
      else logLikelihood += pow((model[i][1] - data[j][1]),2.0);
      i++;
      j++;
    }
    //If first dataset is still before second, full difference
    else if((model[i][0] - data[j][0]) > EPSILON){
      j++;
    }
    //...
    else{
      i++;
    }
  }
  return(logLikelihood);
}




/* ================================= Math Functions ==================================== */

double Handler::poisson_pmf(const double k, const double lambda) {
  return exp(k * log(lambda) - lgamma(k + 1.0) - lambda);
}



// Calculates the standard deviation of a given column of a given set of data
double Handler::calculate_sd(vector<vector<double> > data, int column){
  double mean = calculate_mean(data, column);
  double temp = 0;
  for(unsigned int i = 0; i < data.size(); ++i){
    temp += (mean-data[i][column])*(mean-data[i][column]);
  }
  return(sqrt(temp/(double)data.size()));
}
  


// Returns the mean of a given column of a given set of data
double Handler::calculate_mean(vector<vector<double> > data, int column){
  double total = 0;
  int N = data.size();
  for(int i=0;i<N;++i){
    total += data[i][column];
  }
  return(total/N);
}
  

// Returns the SStot of a given column of a given set of data
double Handler::SStot(vector<vector<double> > data, int column){
  double sum = 0;
  double mean = calculate_mean(data, column);
  for(unsigned int i =0; i<data.size();++i){
    sum += pow((data[i][column]-mean),2.0);
  }
  return(sum);
}



double Handler::logistic(double x, double xmin, double xmax){
  return(((xmax-xmin)/(1+exp(-x)))+xmin);
}

double Handler::logit(double y, double xmin, double xmax){
  return(-log(((xmax-xmin)/(y-xmin))-1));
}
    


double Handler::aicc(double sse, int n, int k){
  if(n == 0 || (n-k-1) == 0) return 99999999999.9;
  double aic = n*log(sse/n) - 2*k;
  return(aic + (2*k*(k+1))/(n-k-1));
}




/* ================================= House Keeping Functions ============================ */

Handler::Handler(){
  useMLE = false;
  optimT0 = true;
  optimI0 = false;
  singleEpi = false;
  plot = true;
  save = false;
  useLogistic = false;
}

Handler::~Handler(){
  cout << "Delete Data Handler" << endl;
  for(unsigned int i = 0;i<epidemics.size();++i){
    if(epidemics[i] != NULL){
      delete epidemics[i];
    }
  }
  for(unsigned int i = 0;i<tempEpidemics.size();++i){
    if(tempEpidemics[i] != NULL){
      delete tempEpidemics[i];
    }
  }
}



void Handler::update_options(bool mle, bool useT0, bool useI0, bool _singleEpi, bool savePlot, bool saveResults,string location, bool useLogist){
  useMLE = mle;
  optimT0 = useT0;
  optimI0 = useI0;
  singleEpi = _singleEpi;
  plot = savePlot;
  save = saveResults;
  saveLocation = location;
  useLogistic = useLogist;
}



// Returns a base model, where each data point is simply the mean of the entire data set
vector<vector<double> > Handler::base_model(vector<vector<double> > data){
  vector<vector<double> > toReturn = current_data;
  double tempMean;
  int numberCols = data[0].size();
  
  // For each column, set all values to the mean of the data
  for(int i=1;i<numberCols;++i){
    tempMean = calculate_mean(data, i);
    for(unsigned int j=0;j<current_data.size();++j){
      toReturn[j][i] = tempMean;
    }
  }
  return(toReturn);
}

// Returns a two dimensional vector of the residuals between two sets of data
vector<vector<double> > Handler::get_residuals(vector<vector<double> > data1, vector<vector<double> > data2, int index){
  unsigned int i = 0;
  unsigned int j = 0;
  vector<vector<double> > results;
  vector<double> row;
  //cout << "Data1:"<<endl;
  //print_vector(data1);
  //  cout << "Data2:" << endl;
  //print_vector(data2);
  while(i < data1.size() && j < data2.size()){
    //If indices give same time point, use difference
    if(data1[i][0] == data2[j][0]){
      row.clear();
      row.push_back(data1[i][0]);
      row.push_back((data1[i][index] - data2[j][index]));
      results.push_back(row);
      i++;
      j++;
    }
    //If first dataset is still before second, full difference
    else if(data1[i][0] < data2[j][0]){
      row.clear();
      row.push_back(data1[i][0]);
      row.push_back(data1[i][index]);
      results.push_back(row);
      i++;
    }
    //...
    else{
      row.clear();
      row.push_back(data2[j][0]);
      row.push_back(data2[j][index]);
      results.push_back(row);
      j++;
    }
  }
  //cout << "Residuals: " << endl;
  //print_vector(results);
  return(results);
}




vector<vector<double> > Handler::combine_vectors(vector<vector<double> > model, vector<vector<double> > data){
  unsigned int i = 0;
  unsigned int j = 0;
  unsigned int k = 0;
  vector<vector<double> > final = model;
  
  while(i < model.size()){
    //cout << model[i][0] << " " << data[j][0] << endl;
    if(fabs(model[i][0] - data[j][0]) < EPSILON){
      final[k][0] = model[i][0];
      final[k][1] = model[i][1] + data[j][1];
      k++;
      i++;
      j++;
    }
    else if(model[i][0] < data[j][0]){
      final[k][0] = model[i][0];
      final[k][1] = model[i][1];
      k++;
      i++;
    }
    else{
      j++;
    }
  }
  return(final);
}

// Prints the epidemic type in standard output
void Handler::print_epidemic_type(EpiType epi){
  switch(epi){
  case none:
    cout << "none";
    break;
  case sir:
    cout << "SIR";
    break;
  case seir:
    cout << "SEIR";
    break;
  case serir:
    cout << "SERIR";
    break;
  case irsir:
    cout << "irSIR";
    break;
  case spike:
    cout << "EXP";
    break;
  default:
    break;
  }
}



/* Quick function to add vector b to the end of vector a and to return the result */
vector<double> Handler::concatenate_vectors(vector<double> a, vector<double> b){
  for(unsigned int i =0; i<b.size();++i){
    a.push_back(b[i]);
  }
  return(a);
}



bool Handler::check_already_tested(vector<EpiType> tested, EpiType toCheck){
  for(unsigned int p = 0; p <tested.size();++p){
    if(toCheck == tested[p]){
      return true;
    }
  }
  return false;

}  

// Returns a copy of the current epidemics without the j'th epidemic
vector<Epidemic*> Handler::fewer_epidemics(int j){
  vector<Epidemic*> epi;
  for(unsigned int i = 0;i<tempEpidemics.size();++i){
    if((int)i != j){
      epi.push_back(tempEpidemics[i]);
    }
  }
  return(epi);
}

vector<int> Handler::fewer_sizes(int j){
  vector<int> sizes1;
  for(unsigned int i = 0;i<tempSizes.size();++i){
    if((int)i != j){
      sizes1.push_back(tempSizes[i]);
    }
  }
  return(sizes1);
}

void Handler::reset_epidemics(){
  vector<Epidemic*> epi;
  for(unsigned int i = 0;i<epidemics.size();++i){
    if(epidemics[i] != NULL){
      epi.push_back(epidemics[i]);
    }
   }
  epidemics = epi;
}

void Handler::remove_epidemic(Epidemic* remove){
  for(unsigned int i = 0;i < epidemics.size();++i){
    if(epidemics[i] == remove){
      epidemics[i] = NULL;
      delete remove;
    }
  }
}

// Adds an epidemic of the provided type to the current vector of epidemics
Epidemic* Handler::new_epidemic(EpiType _newEpidemic, int time){
  Epidemic *additionalEpi;
  switch(_newEpidemic){
  case sir:
    additionalEpi = new SIR(current_data.size(),current_data,_newEpidemic, time);
    break;
  case seir:
    additionalEpi = new SEIR(current_data.size(),current_data,_newEpidemic, time);
    break;
  case spike:
    additionalEpi = new Spike(current_data.size(),current_data,_newEpidemic, time);
    break;
  default:
    additionalEpi = new SIR(current_data.size(),current_data,_newEpidemic, time);
    break;
  }
  return(additionalEpi);
}



/* =============================== DATA HANDLING HELPER FUNCTIONS ================================*/
/* Imports data from a .csv file and saves this to the current data property of the object.
   Returns the number of rows of the imported file as a double, or -1 if the file can not be
   opened.
*/
double Handler::import_data(const char* file){
  ifstream in_stream;
  string temp;
  int rows = 0;
  current_data.clear();
  
  cout << "Importing file... " << file << endl;
  
  // Attempt to open the file
  in_stream.open(file);
  if (in_stream.fail()) {
    cout << "Error opening file" << endl;
    return -1;
  }

  // Get line of data and split using commas. Push this to the 
  while(!in_stream.eof()){
    if(!getline(in_stream,temp)) break;
    istringstream ss(temp);
    vector<double> record;
    while(ss){
      double temp1;
      if(!getline(ss,temp,',')) break;
      stringstream convert(temp);
      if(!(convert >> temp1)) temp1 = 0.0;
      record.push_back(temp1);
    }
    current_data.push_back(record);
    rows++;
  }
  empty_model = create_empty_data_vector(rows);
  cout << "File imported! Number of rows: " << current_data.size() << endl << endl;
  return double(rows);
}

vector<vector<double> > Handler::create_empty_data_vector(int _rows){
  vector<vector<double> > empty;
  empty.resize(_rows);
  for(int i = 0;i<_rows;++i){
    empty[i].resize(2);
    empty[i][0] = i;
    empty[i][1] = 0.0;
  }
  return empty;
}
    

/* Function to print out a two dimensional vector of data */
void Handler::print_vector(vector< vector<double> > my_data){
  for(vector< vector<double> >::const_iterator i = my_data.begin(); i !=my_data.end(); ++i){
    for(vector<double>::const_iterator j = i->begin(); j!=i->end();++j){
      cout << *j << ' ';
    }
    cout << endl;  
  }
}
/* Function to print out a two dimensional vector of data */
void Handler::print_vector(vector< vector<int> > my_data){
  for(vector< vector<int> >::const_iterator i = my_data.begin(); i !=my_data.end(); ++i){
    for(vector<int>::const_iterator j = i->begin(); j!=i->end();++j){
      cout << *j << ' ';
    }
    cout << endl;  
  }
}

vector<Epidemic*> Handler::copy_epidemics(vector<Epidemic*> epi){
  vector<Epidemic*> tempEpi;
  return(tempEpi);
}



/* ================================ GRAPH PLOTTING FUNCTIONS ===================================*/


/* Plots a gnuplot graph with the following components:
             1. Dashed lines for each component sub epidemic
	     2. Solid line of total combined epidemic
	     3. Line with explicit points of all data considered so far
	     4. Title with the current data point index
	     5. Key with the parameters for each sub epidemic 
	     6. XLabel of the overall RSquare fit value
	     7. Vertical lines for epidemic detection times
   
   Note that data must be based as a two dimensional vector, where the first column is the 
   x values and the 3rd column (usually the infecteds) are the y values.
*/
void Handler::plotGraphMulti(vector<vector<vector<double> > > finalResults, vector<vector<double> > totalResults, vector<vector<double> > data, int index, vector<double> parameters, double _RSquare, int column){

  Gnuplot gp;   // Need instance of the Gnuplot class to pipe commands to gnuplot
  string name = "graphs10/output"; // The save location and general name of the graph to be saved
  string _index = to_string((index-1));
  string label, xlab, ylab,xlab2;
  string colours[6] = {"blue","red","orange","cyan","violet","yellow"};
  int j = 0;
  int xmax = current_data[current_data.size()-1][0];

  ylab = "Number of infected individuals";
  xlab = "Time";
  xlab2 = "Run number: " + _index + ". " + "RSquare: " + to_string(_RSquare);  // Xlabel is the RSquare value
  name = name + _index + ".png";  // Graph name

  gp << "set terminal png size 1000,800 enhanced font \"Helvetica, 10\"\n";    // Here edits type, size and formatting of graph
  gp << "set xlabel '" << xlab << "'\n";
  gp << "set x2label '" << xlab2 << "'\n";
  gp << "set ylab '" << ylab << "'\n";
  gp << "set xrange[0:" << xmax << "]\n";
  gp << "set output '" << name << "'\n";  // Set output to the specified file name
  gp << "set termoption dash\n"; // Allow dashes
  gp << "set style rect fc lt -1 fs solid 0.15 noborder \n";
  
  gp << "set obj rect from 0, graph 0 to " + _index + ", graph 1\n";
  // Firstly, plot the actual data (3rd column)
  gp << "plot '-' using 1:2 with linespoints lw 1.5 lt 19 title 'Future Data'";
  gp << ", '-' using 1:2 with linespoints lw 1.5 lt 19 linecolor rgb \"brown\" title 'Current Data'";
  
  
  // For each sub epidemic, add a dashed line to the graph and a corresponding label to the key
  
  if(finalResults.size() > 0){
    label = "Baseline Level";
    gp << ", '-' using 1:2 with lines lw 1.5 lt 0 lc " << 0 << " title '" << label << "'";
  }
  
  if(epidemics.size() > 0 && parameters.size() != 0 && finalResults.size() > 1){
    for(int f = 0;f<(int)epidemics.size();++f){
    // Create the label depending on the type of epidemic
    label = "Sub-Epidemic " + to_string(f+1);
    label += create_label(epidemics[f], parameters, j);
    gp << ", '-' using 1:2 with lines lw 1.5 lt 0 linecolor rgb \"" << colours[f] << "\" title '" << label << "'";
    }
  }
  // Secondly, plot the overall model values
  gp << ", '-' using 1:2 with lines lw 1.5 linecolor rgb \"green\" title 'Total'";
  

  gp << "\n";
  
// Send the data to Gnuplot
 
gp.send1d(current_data);
gp.send1d(data);
//gp.send1d(finalResults[0]);
for(unsigned int j = 0;j<finalResults.size();++j){
  gp.send1d(finalResults[j]);
    }
  gp.send1d(totalResults);
}

string Handler::create_label(Epidemic* epi, vector<double> par1, int& i){
  string toReturn;
  switch(epi->return_type()){
  case sir:
    toReturn = " (SIR): beta = " + to_string(par1[i]) + "; gamma = " + to_string(par1[i+1]) + "; S0 = " + to_string(par1[i+2]) + "; t0 = " + to_string(par1[i+3]);
    i = i+4;
    break;
  case seir:
    toReturn = " (SEIR): beta = " + to_string(par1[i]) + "; alpha = " + to_string(par1[i+1]) +"; gamma = " + to_string(par1[i+2]) + "; S0 = " + to_string(par1[i+3]) + "; t0 = " + to_string(par1[i+4]);
    i = i+5;
    break;
  case spike:
    toReturn = " (EXP): gamma = " + to_string(par1[i]) + "; I0 = " + to_string(par1[i+1]) + "; t0 = " + to_string(par1[i+2]);
    i = i+3;
    break;
  case serir:
    toReturn = " (SERIR): beta = " + to_string(par1[i]) + "; alpha = " + to_string(par1[i+1]) + "; mu = " + to_string(par1[i+2]) + "; gamma = " + to_string(par1[i+3]) + "; S0 = " + to_string(par1[i+4]) + "; t0 = " + to_string(par1[i+5]);
    i = i+6;
    break;
  case irsir:
    toReturn = " (irSIR): beta = " + to_string(par1[i]) + "; gamma = " + to_string(par1[i+1]) + "; S0 = " + to_string(par1[i+2]) + "; t0 = " + to_string(par1[i+3]);
    i = i+4;
    break;
  default:
    break;
  
  }
  return(toReturn);
}


/* Older function to plot graph of only data and the combined results */
void Handler::plotGraph(vector<vector<double> > finalResults, vector<vector<double> > data, int index){
  Gnuplot gp;
  string name = "graphs7/output";
  string _index = to_string(index);
  name = name + _index + ".jpeg";
 
  gp << "set terminal jpeg size 1000,800 enhanced font \"Helvetica, 20\"\n";
  gp << "set output '" << name << "'\n";
  gp << "plot '-' using 1:2 with lines title 'I', '-' using 1:2 with linespoints title 'Data'\n";
  gp.send1d(finalResults);
  gp.send1d(data);
}