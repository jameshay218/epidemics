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

Handler::Handler(){
  useMLE = false;
  optimT0 = false;
  optimI0 = false;
  singleEpi = false;
  plot = true;
  save = false;
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


void Handler::update_options(bool mle, bool useT0, bool useI0, bool _singleEpi, bool savePlot, bool saveResults, string location){
  useMLE = mle;
  optimT0 = useT0;
  optimI0 = useI0;
  singleEpi = _singleEpi;
  plot = savePlot;
  save = saveResults;
  saveLocation = location;
}




/* ================================ FINAL LEAST SQUARES ============================= */

void Handler::realtime_fit_single(double targetRSq, EpiType _epi){
  clock_t t1, t2, t3;
  Epidemic *tempEpi,*bestFit;
  vector<EpiType> candidateModels;
  vector<double> finalParams,tempPar;
  vector<EpiType> selectedTypes;
  vector<vector<double> > combinedResults, tempCombinedResults, residuals, saveResults; // 
  vector<vector<vector<double> > > tempAll, tempAll1;
  double RSquare,tempRSquare, RMSE;
  double SSE = 99999999999.9;
  double tempSSE = 99999999999.9;
  double currentBestAICc, AICc, tempAICc = 999999999999.9;
  int iterations = 0;
  int tempIterations = 0;
  t1=clock();
  srand(time(NULL));

  //  candidateModels.push_back(spike);
  candidateModels.push_back(sir);
  candidateModels.push_back(serir);
  candidateModels.push_back(irsir);
  candidateModels.push_back(seir);
 
  // Start off with a baseline model of the mean of the first 4 points
  finalParams = generate_seed_parameters();
  for(unsigned int j = 0; j <= 4; ++j){
    this->temp_data.push_back(this->current_data[j]);
  }
  this->baseModel = this->current_model = base_model(this->temp_data);


  if(_epi != none) epidemics.push_back(new_epidemic(_epi,1,1.0));
  cout << "Considering: ";
  print_epidemic_type(_epi);

  finalParams = tempPar = generate_seed_parameters();
  for(unsigned int i = 7;i<this->current_data.size();++i){
    cout << endl << "-----------------" << endl;
    cout << "Iteration number " << i << endl;
    
    // Clear temp vector and temp epidemics
    this->temp_data.clear();

    // Store all data up to the current index and update epidemics with data
    for(unsigned int j = 0; j < i; ++j){
      this->temp_data.push_back(this->current_data[j]);
    } 
    if(_epi == none && (epidemics.size() == 0 || i%10 == 0)){
      if(epidemics.size() != 0) remove_epidemic(epidemics[0]); reset_epidemics();
      bestFit = NULL;
      currentBestAICc = 99999999999999999999.9;
      cout << "Attempting to find best fitting model" << endl;
      cout << "Starting with: " << epidemics.size() << " epidemics" << endl;
      for(unsigned int k = 0; k<candidateModels.size();++k){
	cout << "Considering: ";
	//currentBestSSE = 999999999999.9;
	print_epidemic_type(candidateModels[k]);
	cout << endl;
	tempEpi = NULL;
	tempEpi = new_epidemic(candidateModels[k],1,1.0);
	epidemics.push_back(tempEpi);
	tempSSE=optimise_single(tempPar, tempCombinedResults,tempAll1, tempIterations);
	tempRSquare = 1 - tempSSE/SStot(temp_data, 1);
	cout << "Parameters: " << tempEpi->return_parameters().size() << endl;
	cout << "Temp SSE: " << tempSSE << endl;
	cout << "Temp RSquare: " << tempRSquare << endl;
	tempAICc = aicc(tempSSE, temp_data.size(), tempEpi->return_parameters().size());
	cout << "Model AICc: " << tempAICc << endl;
	//if((tempRSquare > targetRSq && tempSSE < currentBestSSE) || (tempRSquare > targetRSq && bestFit != NULL && (tempEpi->return_parameters().size() < bestFit->return_parameters().size())) || (bestFit != NULL && tempEpi->return_parameters().size() == bestFit->return_parameters().size() && tempSSE < currentBestSSE && tempIterations < iterations)){
	if(tempRSquare > targetRSq && tempAICc < currentBestAICc){
	  cout << "Model improvement" << endl;
	  if(bestFit != NULL) delete bestFit;
	  bestFit = new_epidemic(candidateModels[k],1,1.0);
	  SSE = tempSSE; 
	  currentBestAICc = tempAICc;
	  RSquare = tempRSquare;
	  finalParams = tempPar;
	  combinedResults = tempCombinedResults;
	  tempAll = tempAll1;
	  iterations = tempIterations;
	}
	remove_epidemic(tempEpi);
	reset_epidemics();
      }
      cout << "Fitting complete" << endl;
      if(bestFit != NULL){
	cout << "Adding epidemic type: ";
	print_epidemic_type(bestFit->return_type());
	selectedTypes.push_back(bestFit->return_type());
	cout << endl;
	epidemics.push_back(bestFit);
	bestFit = NULL;
      }
      reset_epidemics();
    }
    else SSE = optimise_single(finalParams, combinedResults, tempAll, iterations);
    if(epidemics.size() == 0){ cout << "No decent fit. Skipping to next data point." << endl; continue;}
    cout << "Number of epidemics: " << epidemics.size() << endl;
    RSquare = 1 - SSE/SStot(temp_data, 1);
    AICc = aicc(SSE, temp_data.size(), epidemics[0]->return_parameters().size());
    // Transform parameters back to normal space
    for(unsigned int j=0;j<finalParams.size();++j){
      finalParams[j] = exp(finalParams[j]);
    }
    // Plot graph
    if(plot) plotGraphMulti(tempAll, combinedResults, this->temp_data, i, finalParams, RSquare, 2);    
    
    // Save results
    if(save) cout << "Save results" << endl;
    
    // Print out results for this iteration
    cout << "Final SSE was: " << SSE << endl;
    cout << "Final RSquare: " << RSquare << endl;
    cout << "Final AICc: " << AICc << endl;
    cout << "Final parameters: ";
    //finalParams.push_back(double(t1));
    printcon(finalParams) ;
    finalParams.push_back(RSquare);
    finalParams.push_back(RMSE);
    saveResults.push_back(finalParams);
    RMSE = sqrt((SSE/temp_data.size()));
    cout << "RMSE: " << RMSE << endl;
    cout << endl << "-----------------" << endl;
  }
  t2 = clock();
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

  cout << "Fitting procedure complete" << endl;
}


double Handler::optimise_single(vector<double> &parameters, vector<vector<double> > &results, vector<vector<vector<double> > > &allResults, int& itr){
  double SSE = 99999999999999.9;
  double tempSSE;
  vector<double> tempPar, seedParams;
  Simplex simplex;

  parameters.clear();
  allResults.clear();
  parameters = generate_seed_parameters();
  for(int index = 0;index<20;index++){
    seedParams.clear();
    seedParams = generate_seed_parameters();
    //printcon(seedParams);
    if(this->useMLE == false) tempPar = simplex.neldermead(&Handler::fitEpidemicsMLE, *this,  seedParams, itr);
    else tempPar = simplex.neldermead(&Handler::fitEpidemics, *this,  seedParams, itr);
    // Store the SSE value for this
    tempSSE = fitEpidemics(tempPar);
    if(tempSSE < SSE){
      SSE = tempSSE;
      parameters=tempPar;
    }
    cout << "." << flush;
  }
  results = ode_solve(parameters);
  allResults = ode_solve_separate(parameters);
  //allResults.clear();
  cout << endl;
  return SSE;
}


void Handler::realtime_fit_multi(double targetRSq){
  clock_t t1, t2; // To record the total run time
  EpiType newEpidemic; // Keep track of the type of a new epidemic to be added
  Epidemic *toRemove, *bestFit; // Keep track of a newly created epidemic to allow deletion
  vector<EpiType> candidateModels;
  vector<double> finalParams, finalParamsTemp; // Keep track of parameters from optimisation
  vector<vector<double> > combinedResults, residuals, combinedResultsTemp; // 
  vector<vector<vector<double> > > componentResults, componentResultsTemp;
  double RSquare, tempRSquare, SSE, tempSSE, currentBestSSE;
  int iterations, tempIterations;
  double bestAICc, tempAICc;
  srand(clock());

  //candidateModels.push_back(spike);
  candidateModels.push_back(sir);
  candidateModels.push_back(seir);

  // Start of model fitting process
  t1=clock();
  epidemics.push_back(new_epidemic(sir, 1, 1.0));
  epidemics.push_back(new_epidemic(sir, 4, 4.0));
  //add_epidemic(sir,20);
  // Start off with a baseline model of the mean of the first 4 points
  for(unsigned int j = 0; j <= 4; ++j){
    this->temp_data.push_back(this->current_data[j]);
  }
  this->baseModel = this->current_model = base_model(this->temp_data);

  // For each time point in the current data set, carry out the fitting procedure
  for(unsigned int i = 100;i<this->current_data.size();++i){
    cout << endl << "-----------------" << endl;
    cout << "Iteration number " << i << endl;
    
    // Clear temp vector and temp epidemics
    this->temp_data.clear();
    this->tempEpidemics.clear();
    toRemove = NULL;

    // Store all data up to the current index and update epidemics with data
    for(unsigned int j = 0; j < i; ++j){
      this->temp_data.push_back(this->current_data[j]);
    } 
    for(unsigned int f = 0;f<this->epidemics.size();++f){
      this->epidemics[f]->update_data(this->temp_data);
    }

    // Create an optimised model fit for the given data (temp) and calculate RSquare from the SSE
    SSE = optimiseEpidemics(finalParams, combinedResults, componentResults, iterations);
    RSquare = 1 - SSE/SStot(temp_data, 1);
    cout << "Optimised fit, RSquare: " << RSquare << endl;
    // If fit is sufficiently good with k epidemics, try fitting with k-1 epidemics. If this fit
    // is sufficient, use k-1 epidemics.
    if(RSquare > targetRSq && this->epidemics.size() > 0){
      cout << "Considering removal of an epidemic" << endl;
      // Store current set of epidemics
      this->tempEpidemics = this->epidemics;
      
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
	  currentBestSSE = 999999999.9;
	
	  // Calculate the fit without this epidemic
	  tempSSE = optimiseEpidemics(finalParamsTemp, combinedResultsTemp, componentResultsTemp, tempIterations);
	  tempRSquare = 1 - tempSSE/SStot(temp_data,1);
	
	  // If sufficient fit from removing an epidemic, permanently remove epidemic from set
	  // Otherwise, restore epidemics to full number. Also track which removal was best
	  if(tempRSquare > targetRSq && tempSSE < currentBestSSE){
	    toRemove = this->tempEpidemics[j];
	    cout << "Fit sufficient. Removing an epidemic" << endl;
	    SSE = tempSSE;
	    currentBestSSE = tempSSE;
	    RSquare = tempRSquare;
	    finalParams = finalParamsTemp;
	    combinedResults = combinedResultsTemp;
	    componentResults = componentResultsTemp;
	    iterations = tempIterations;
	  }
	  else{
	    this->epidemics = this->tempEpidemics;
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
      cout << "Infecteds unaccounted for: " << residuals[i-1][1] << endl;
      bestFit = NULL;
      bestAICc = 999999999999.9;
	//currentBestSSE = 999999999.9;
      // Temporarily store the current epidemic state and add a new epidemic to the current set
      this->tempEpidemics = this->epidemics;
      for(unsigned int f = 0; f < candidateModels.size(); ++f){
	
	toRemove = NULL;
	newEpidemic = candidateModels[f];
	cout << "Considering addition of ";
	print_epidemic_type(newEpidemic);
	cout << endl;
	
	toRemove = new_epidemic(newEpidemic, i, residuals[i-1][1]);
	this->epidemics.push_back(toRemove);
	cout << "Epidemics size: " << epidemics.size() << endl;
	// Produce optimised fit with added epidemic
	tempSSE = optimiseEpidemics(finalParamsTemp, combinedResultsTemp, componentResultsTemp, tempIterations);
	tempAICc = aicc(tempSSE, temp_data.size(), toRemove->return_parameters().size());
	cout << "AICc: " << tempAICc << endl;
	tempRSquare = 1 - tempSSE/SStot(this->temp_data,1);
	
	cout << "Fit with additional epidemic: " << tempRSquare << endl;
	cout << "Iterations taken: " << tempIterations << endl;
	// If fit improved, keep track of epidemic. If best fit, add to epidemics at the end
	//	if((tempRSquare > RSquare && tempSSE < currentBestSSE) || (tempRSquare > targetRSq && bestFit != NULL && (toRemove->return_parameters().size() < bestFit->return_parameters().size())) || (bestFit != NULL && toRemove->return_parameters().size() == bestFit->return_parameters().size() && tempSSE < currentBestSSE && tempIterations < iterations)){
	if(tempRSquare > RSquare && tempAICc < bestAICc){
	cout << "Fit improved with additional epidemic. Add epidemic" << endl;
	  if(bestFit != NULL) delete bestFit;
	  bestFit = new_epidemic(newEpidemic,i, residuals[i-1][1]);
	  cout << "Better AICc value" << endl;
	  SSE = tempSSE;
	  //currentBestSSE = tempSSE;
	  bestAICc = tempAICc;
	  RSquare = tempRSquare;
	  finalParams = finalParamsTemp;
	  combinedResults = combinedResultsTemp;
	  componentResults = componentResultsTemp;
	  iterations = tempIterations;
	}
	remove_epidemic(toRemove);
	reset_epidemics();
      }
      if(bestFit != NULL){
	  cout << "Adding epidemic type: ";
	  print_epidemic_type(bestFit->return_type());
	  cout << endl;
	  epidemics.push_back(bestFit);
	  bestFit = NULL;
	}
      }
      // Reset pointers in case epidemic was deleted;
    reset_epidemics();
    // Transform parameters back to normal space
    for(unsigned int j=0;j<finalParams.size();++j){
      finalParams[j] = exp(finalParams[j]);
    }
    // Plot graph
    if(plot) plotGraphMulti(componentResults, combinedResults, this->temp_data, i, finalParams, RSquare, 2);    
    
    // Save results
    if(save) cout << "Save results" << endl;
    
    // Print out results for this iteration
    cout << "Final SSE was: " << SSE << endl;
    cout << "Final RSquare: " << RSquare << endl;
    cout << "Final parameters: ";
    printcon(finalParams) ;
    cout << endl << "-----------------" << endl;
    this->current_model = combinedResults;
  }
  t2=clock();
  cout << "Time taken: " << (t2-t1)/CLOCKS_PER_SEC << endl << endl;
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
Epidemic* Handler::new_epidemic(EpiType _newEpidemic, int time, double infected){
  Epidemic *additionalEpi;
  switch(_newEpidemic){
  case sir:
    additionalEpi = new SIR(current_data.size(),current_data,_newEpidemic, time, infected, optimT0, optimI0);
    break;
  case seir:
    additionalEpi = new SEIR(current_data.size(),current_data,_newEpidemic, time, infected, optimT0, optimI0);
    break;
  case spike:
    additionalEpi = new Spike(current_data.size(),current_data,_newEpidemic, time, infected, optimT0, optimI0);
    break;
  case serir:
    additionalEpi = new SERIR(current_data.size(),current_data,_newEpidemic, time, infected, optimT0, optimI0);
    break;
  case irsir:
    additionalEpi = new IRSIR(current_data.size(),current_data,_newEpidemic, time, infected, optimT0, optimI0);
    break;
  default:
    additionalEpi = new SIR(current_data.size(),current_data,_newEpidemic, time, infected, optimT0, optimI0);
    break;
  }
  return(additionalEpi);
}

double Handler::optimiseEpidemics(vector<double> &parameters, vector<vector<double> > &results, vector<vector<vector<double> > > &allResults, int& itr){
  double SSE = 99999999999.9;
  double tempSSE;
  vector<double> tempPar, seedParams;
  vector<vector<vector<double> > > tempAll;
  Simplex simplex;
  int iterations = 10001;
  // If no epidemics have yet been detected, use the mean of the current data as the 
  // current model and return the corresponding SSE.
  if(epidemics.size() == 0){
    results = this->baseModel;
    SSE = calculate_SSE(results, temp_data);
    parameters.clear();
    allResults.clear();
    allResults.push_back(results);
    return(SSE);
  }
  parameters = generate_seed_parameters();
  // If there are epidemics to be fitted, perform 40 random fits and keep best fitting model
  for(int index=0;index<40;index++){    
            
    // Clear the temporary seed parameters and seed rand
   

    // Create a list of random seed parameters
    seedParams.clear();      
    seedParams = generate_seed_parameters();
    
    // Get the optimised parameters from nelder mead algorithm
    if(this->useMLE == false) tempPar = simplex.neldermead(&Handler::fitEpidemicsMLE, *this,  seedParams, iterations);
    else tempPar = simplex.neldermead(&Handler::fitEpidemics, *this,  seedParams, iterations);
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
 
  // Get the combined values from these parameters, as well as a vector of 
  // each sub-epidemic
  results = ode_solve(parameters);
  
  tempAll = ode_solve_separate(parameters);
  
  allResults.clear();
  allResults.push_back(this->baseModel);
  
  for(unsigned int x = 0;x<tempAll.size();++x){
    allResults.push_back(tempAll[x]);
  }
  
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
  current_model = create_empty_data_vector(temp_data.size());
  //current_model = baseModel;
  for(unsigned int i = 0;i<epidemics.size();++i){
    tempParams.clear();
    temp_model.clear();
    for(unsigned int x = 0;x<epidemics[i]->return_parameters().size();++x){
      tempParams.push_back(params[z]);
      z++;
    }
    
    //if(!params_check(tempParams, epidemics[i]->return_type())) return(999999999999.9);
    temp_model = epidemics[i]->ode_solve(tempParams);  
    current_model = combine_vectors(current_model, temp_model);
  }
   return(calculate_SSE(current_model, temp_data));
 }

 double Handler::fitEpidemicsMLE(vector<double> params){
   int z = 0;
   current_model = create_empty_data_vector(temp_data.size());
   for(unsigned int i = 0;i<epidemics.size();++i){
     tempParams.clear();
     temp_model.clear();
     for(unsigned int x = 0;x<epidemics[i]->return_parameters().size();++x){
       tempParams.push_back(params[z]);
       z++;
     }
     temp_model = epidemics[i]->ode_solve(tempParams);
     current_model = combine_vectors(current_model, temp_model);
   }
   return(dpois(current_model, temp_data));
 }

 // As above, but returns the overall model.
 vector<vector<double> > Handler::ode_solve(vector<double> params){
   vector<vector<double> > overallResults = create_empty_data_vector(temp_data.size());
   //vector<vector<double> > overallResults = baseModel;
   int z = 0;
   for(unsigned int i = 0;i<epidemics.size();++i){
     tempParams.clear();
     temp_model.clear();
     for(unsigned int x = 0;x<epidemics[i]->return_parameters().size();++x){
       tempParams.push_back(params[z]);
       z++;
     }

     temp_model = epidemics[i]->ode_solve(tempParams);
     overallResults = combine_vectors(overallResults,temp_model);
   }
   
   return(overallResults);
 }

 // As above, but returns a vector of the component models.
 vector<vector<vector<double> > > Handler::ode_solve_separate(vector<double> params){
   vector<vector<vector<double> > > separateResults;
   for(unsigned int i = 0;i<epidemics.size();++i){
     tempParams.clear();
     temp_model.clear();
     for(unsigned int z = 0;z<epidemics[i]->return_parameters().size();++z){
       tempParams.push_back(epidemics[i]->return_parameters()[z]);
     }
     temp_model = epidemics[i]->ode_solve(tempParams);
     separateResults.push_back(temp_model);
   }
   return(separateResults);
 }



 /* Takes the residuals as a 2xN vector. Returns true if the latest X residuals are a certain number
    of standard deviations away from the previous total-X residuals */
 EpiType Handler::check_epidemic(vector<vector<double> > residuals){
   EpiType detected = none;
   double mean, sd, SIRresidual, EXPresidual, sirlimit,explimit;
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

   sirlimit = mean + 2*sd;
   explimit = mean + 6*sd;

   SIRresidual = residuals[residuals.size()-1][1];
   EXPresidual = residuals[residuals.size()-1][1];

   for(unsigned int i = residuals.size()-1; i > previousResiduals.size();--i){
     if((residuals[i][1] - EXPresidual) < 0) EXPresidual = residuals[i][1];
     //if(residuals[i][1] < SIRresidual) SIRresidual = residuals[i][1];
   }

   double diffEXP = EXPresidual - explimit;
   double diffSIR = SIRresidual - sirlimit;
   if(diffEXP > 0){ detected = spike;}
   else if(diffSIR > 0){detected = sir;}

   return detected;
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


 /* ================================= MATHEMATICAL FUNCTIONS =====================================*/
double Handler::aicc(double sse, int n, int k){
  if(n == 0 || (n-k-1) == 0) return 99999999999.9;
  double aic = n*log(sse/n) - 2*k;
  return(aic + (2*k*(k+1))/(n-k-1));
}



 bool Handler::params_check(vector<double> pars, EpiType epi){
   switch(epi){
   case spike:
     if(pars[0] < 0.00001 || pars[0] > 1) return false;
     if(pars[1] < 0 || pars[1] > 50000) return false;
     if(optimT0 && (pars[0] < 0 || pars[0] > 300)) return false;
     break;
   case sir:
     if(pars[0] < 0.00001 || pars[0] > 1) return false;
     if(pars[1] > 1 || pars[1] <= pars[0]) return false;
     if(pars[2] < 0 || pars[2] > 50000 || (optimI0 && pars[2] <= pars[3])) return false;
     if(optimI0 && (pars[3] < 0 || pars[3] > 50000)) return false;
     if(optimT0){
       if(!optimI0 && (pars[4] < 0 || pars[4] > 300)) return false;
       else if(pars[3] < 0 || pars[3] > 300) return false;
     }
     break;
   case seir:
     if(pars[0] < 0.00001 || pars[0] > 1) return false;
     if(pars[1] < 0.00001 || pars[1] > 1) return false;
     if(pars[2] > 1 || pars[2] <= pars[0]) return false;
     if(pars[3] < 0 || pars[3] > 50000 || (optimI0 && pars[3] <= pars[4])) return false;
     if(optimI0 && (pars[4] < 0 || pars[4] > 50000)) return false;
     if(optimT0){
       if(!optimI0 && (pars[5] < 0 || pars[5] > 300)) return false;
       else if(pars[4] < 0 || pars[4] > 300) return false;
     }
     break;
   case serir:
     if(pars[0] < 0.00001 || pars[0] > 1) return false;
     if(pars[1] < 0.00001 || pars[1] > 1) return false;
    if(pars[2] < 0.00001 || pars[2] > 1) return false;
    if(pars[3] > 1 || pars[3] <= pars[0]) return false;
    if(pars[4] < 0 || pars[4] > 50000 || (optimI0 && pars[4] <= pars[5])) return false;
    if(optimI0 && (pars[5] < 0 || pars[5] > 50000)) return false;
    if(optimT0){
      if(!optimI0 && (pars[6] < 0 || pars[6] > 300)) return false;
      else if(pars[5] < 0 || pars[5] > 300) return false;
    }
    break;
  case irsir:
    if(pars[0] < 0.00001 || pars[0] > 1) return false;
    if(pars[1] > 1 || pars[1] <= pars[0]) return false;
    if(pars[2] < 0 || pars[2] > 50000 || (optimI0 && pars[2] <= pars[3])) return false;
    if(optimI0 && (pars[3] < 0 || pars[3] > 50000)) return false;
    if(optimT0){
      if(!optimI0 && (pars[4] < 0 || pars[4] > 300)) return false;
      else if(pars[3] < 0 || pars[3] > 300) return false;
    }
    break;
  default:
    break;
  }
  return true;
}
  

vector<double> Handler::generate_seed_parameters(){
  vector<double> params;
  for(unsigned int i = 0;i<epidemics.size();++i){
    params = concatenate_vectors(params, rand_params(epidemics[i]->return_type()));
    
    //if(epidemics[i]->return_type() == spike || optimI0) {
    //params.push_back(log(epidemics[i]->return_infecteds()));
    //}
    if(optimT0) params.push_back(log(epidemics[i]->return_detection_time()));
  }
  return(params);
}


// Generates a vector of random parameters
// THIS IS WHERE WE GENERATE SEED VALUES
vector<double> Handler::rand_params(EpiType _type){
  vector<double> params;
  double beta,gamma,s0,i0,alpha,mu;
  
  setprecision(9);
 
  s0 = rand()%100000+1000;
  i0 = rand()%20+1;

  switch(_type){
  case sir:
    beta = (rand()%100+1)/10000.0;
    gamma = (rand()%500+beta)/1000.0;
    params.push_back(log(beta));
    params.push_back(log(gamma));
    params.push_back(log(s0));
    break;
  case seir:
    beta = (rand()%500+10)/10000.0;
    gamma = (rand()%200+10)/10000.0;
    alpha = (rand()%500+10)/10000.0;
    params.push_back(log(beta));
    params.push_back(log(alpha));
    params.push_back(log(gamma));
    params.push_back(log(s0));
    break;
    return(params);
  case serir:
    beta = (rand()%100+1)/10000.0;
    gamma = (rand()%200+10)/10000.0;
    alpha = (rand()%100+1)/10000.0;
    mu = (rand()%200+1)/10000.0;
    params.push_back(log(beta));
    params.push_back(log(alpha));
    params.push_back(log(mu));
    params.push_back(log(gamma));
    params.push_back(log(s0));
    break;
  case irsir:
    beta = (rand()%100+1)/100000.0;
    gamma = (rand()%100+1)/100000.0;
    params.push_back(log(beta));
    params.push_back(log(gamma));
    params.push_back(log(s0));
    break;
  case spike:
    gamma = (rand()%200+1)/1000.0;
    params.push_back(log(gamma));
    params.push_back(log(s0));
    break;
  default:
    break;
  }
  if(optimI0 && _type != spike){
    params.push_back(log(i0));
  }
  return(params);
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
  string name = saveLocation + "/output"; // The save location and general name of the graph to be saved
  string _index = to_string((index-1));
  string label, xlab;
  string colours[5] = {"orange","yellow","cyan","darkblue","violet"};
  int j = 0;

  xlab = "RSquare: " + to_string(_RSquare) + ". Run number: " + _index;  // Xlabel is the RSquare value
  name = name + _index + ".png";  // Graph name

  gp << "set terminal png size 1000,800 enhanced font \"Helvetica, 11\"\n";    // Here edits type, size and formatting of graph
  gp << "set xlabel '" << "Time" << "'\n";
  gp << "set x2label '" << xlab << "'\n";
  gp << "set ylabel '" << "Number of Infected Individuals" << endl;
  gp << "set output '" << name << "'\n";  // Set output to the specified file name
  gp << "set termoption dash\n"; // Allow dashes
  gp << "set style rect fc lt -1 fs solid 0.15 noborder \n";
  
  gp << "set obj rect from 0, graph 0 to " + _index + ", graph 1\n";
  // Firstly, plot the actual data (3rd column)
  gp << "plot '-' using 1:2 with linespoints lw 1.5 lt 19 title 'Future Data'";
  gp << ", '-' using 1:2 with linespoints lw 1.5 lt 19 linecolor rgb \"brown\" title 'Current Data'";
  
  
  // For each sub epidemic, add a dashed line to the graph and a corresponding label to the key
  if(finalResults.size() > 1){
    label = "Baseline Level";
    gp << ", '-' using 1:2 with lines lw 1.5 lt 0 lc " << 0 << " title '" << label << "'";
  }
 

  for(int f = 0;f<(int)epidemics.size();++f){
    // Create the label depending on the type of epidemic
    label = "Sub-Epidemic " + to_string(f+1);
    label += create_label(epidemics[f], parameters, j);
    gp << ", '-' using 1:2 with lines lw 1.5 lt 0 linecolor rgb \"" << colours[f] << "\" title '" << label << "'";
  }

  // Secondly, plot the overall model values
  gp << ", '-' using 1:2 with lines lw 1.5 linecolor rgb \"green\" title 'Total'";
  

  gp << "\n";
  
  // Send the data to Gnuplot
 
  gp.send1d(current_data);
  gp.send1d(data);
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
  string name = "graphs/output";
  string _index = to_string(index);
  name = name + _index + ".jpeg";
 
  gp << "set terminal jpeg size 1000,800 enhanced font \"Helvetica, 20\"\n";
  gp << "set output '" << name << "'\n";
  gp << "plot '-' using 1:2 with lines title 'I', '-' using 1:2 with linespoints title 'Data'\n";
  gp.send1d(finalResults);
  gp.send1d(data);
}



/* Quick function to add vector b to the end of vector a and to return the result */
vector<double> Handler::concatenate_vectors(vector<double> a, vector<double> b){
  for(unsigned int i =0; i<b.size();++i){
    a.push_back(b[i]);
  }
  return(a);
}





/* Calculates and returns the sum of squared errors between to sets of data. Checks 
   the first column to see if on the same time point, and finds squared difference if
   so. If one dataset is at an earlier time point, use the entire squared value as
   the squared residual.
*/
double Handler::calculate_SSE(vector<vector<double> > data1, vector<vector<double> > data2){
  unsigned int i = 0;
  unsigned int j = 0;
  double sse = 0;
  vector<vector<double> > big, small;

  if(data1.size() < data2.size()){
    small = data1;
    big = data2;
  }
  else{
    small = data2;
    big = data1;
  }

  while(j < small.size()){
    //If indices give same time point, find sse
    if(fabs(small[i][0]- big[j][0])<EPSILON){
      sse += pow((small[i][1] - big[j][1]),2.0);
      i++;
      j++;
    }
    //If first dataset is still before second, full difference
    else if(small[i][0] < big[j][0]){
      i++;
    }
    //...
    else{
      j++;
    }
  }
  return(sse);
}



/* Returns the negative log likelihood of a model given a set of data */
double Handler::dpois(vector<vector<double> > model, vector<vector<double> > data){
  unsigned int i = 0;
  unsigned int j = 0;
  double logLikelihood = 0;
  vector<vector<double> > big, small;

  if(model.size() < data.size()){
    small = model;
    big = data;
  }
  else{
    small = data;
    big = model;
  }

  while(j < small.size()){
    //If indices give same time point, find sse
    if(fabs(small[i][0]- big[j][0])<EPSILON){
      if(small[i][0] == 0 || big[j][0] == 0) logLikelihood += log(1);
      else logLikelihood += pow((small[i][1] - big[j][1]),2.0);
      i++;
      j++;
    }
    //If first dataset is still before second, full difference
    else if(small[i][0] < big[j][0]){
      i++;
    }
    //...
    else{
      j++;
    }
  }
  return(logLikelihood);
}


double Handler::poisson_pmf(const double k, const double lambda) {
  return exp(k * log(lambda) - lgamma(k + 1.0) - lambda);
}






void Handler::overall_test(double targetRsq){
  Epidemic* toRemove1;
  EpiType newEpidemic;
  vector<double> finalParams, finalParamsTemp;
  vector<vector<double> > baseModel, combinedResults, residuals, combinedResultsTemp, temp;
  vector<vector<vector<double> > > componentResults, componentResultsTemp;
  double SSE;
  int iter;
  newEpidemic = sir;
  
  toRemove1 = new_epidemic(newEpidemic, 0, 1);
  epidemics.push_back(toRemove1);
  
  for(unsigned int f = 0;f<epidemics.size();++f){
    epidemics[f]->update_data(current_data);
  }
  
  cout << "Data size: " << current_data.size() << endl;
 
  temp_data = current_data;
  finalParams.push_back(log(0.001));  
  //finalParams.push_back(log(0.01));  
  //finalParams.push_back(log(0.01));  
  finalParams.push_back(log(0.1));  
  finalParams.push_back(log(500));
  finalParams.push_back(log(1));
  // temp = epidemics[0]->ode_solve(finalParams);
  temp=ode_solve(finalParams);
  printcon(finalParams);
  finalParams[3] = log(1);
  printcon(finalParams);
  baseModel = ode_solve(finalParams);

  cout << "Temp size: " << temp.size() << endl;
  //  print_vector(temp);
  //  print_vector(temp_data);
  //return;
  cout << dpois(temp,baseModel) << endl;
  plotGraph(temp,temp_data,1);
  return;
  SSE = optimiseEpidemics(finalParams, combinedResults, componentResults, iter);
  
  cout << SSE << endl;
  for(unsigned int i = 0;i<finalParams.size();++i){
    finalParams[i] = exp(finalParams[i]);
  }
  printcon(finalParams);
  //print_vector(combinedResults);
  //print_vector(temp_data);
  plotGraphMulti(componentResults, combinedResults, temp_data, 1, finalParams, 1, 2);    
 
  return;
  finalParams.push_back(log(0.001));  
  finalParams.push_back(log(0.1));  
  //finalParams.push_back(log(0.05));  
  finalParams.push_back(log(500));
  finalParams.push_back(log(10));
  
  //temp = epidemics[0]->ode_solve(finalParams);
  temp=ode_solve(finalParams);
  cout << "Temp size: " << temp.size() << endl;
  //print_vector(temp);
  /*for(int i = 0; i < 100;++i){
    cout << temp[i][0] << " " << temp[i][1] << endl;
    }*/
  //print_vector(temp);
  cout << calculate_SSE(temp,current_data) << endl;
  plotGraph(temp,current_data,1);
  
  for(int i =0;i < 100;++i){
    baseModel.push_back(temp[i]);
  }
  //print_vector(baseModel);
  //print_vector(current_data);
  return;
  temp = ode_solve(finalParams);
  print_vector(temp);
  current_data = temp;
  SSE = fitEpidemics(finalParams);
  cout << "SSE: " << SSE << endl;
  plotGraph(temp, current_data, 1);
}

