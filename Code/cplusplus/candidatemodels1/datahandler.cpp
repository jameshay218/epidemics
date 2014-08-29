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
  optimT0 = true;
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

vector<double> Handler::convert_parameters_back(vector<double> pars, int number, int maxTime, int minTime){
  vector<double> temporary;
  
  if(number == 3){
    temporary.push_back(logistic(pars[0],0.01,0.2));
    temporary.push_back(logistic(pars[1],100,15000));
    temporary.push_back(logistic(pars[2],minTime,maxTime));
  }
  else if(number == 4){
    temporary.push_back(logistic(pars[0],0.0001,0.02));
    temporary.push_back(logistic(pars[1],0.01,0.2));
    temporary.push_back(logistic(pars[2],100,15000));
    temporary.push_back(logistic(pars[3],minTime,maxTime));
  }
  else if(number == 5){
    temporary.push_back(logistic(pars[0],0.0001,0.02));
    temporary.push_back(logistic(pars[1],0.0001,0.1));
    temporary.push_back(logistic(pars[2],0.01,0.2));
    temporary.push_back(logistic(pars[3],100,15000));
    temporary.push_back(logistic(pars[4],minTime,maxTime));
  }
  else{
    temporary.push_back(logistic(pars[0],0.0001,0.02));
    temporary.push_back(logistic(pars[1],0.0001,0.1));
    temporary.push_back(logistic(pars[2],0.0001,0.1));
    temporary.push_back(logistic(pars[3],0.01,0.2));
    temporary.push_back(logistic(pars[4],100,15000));
    temporary.push_back(logistic(pars[5],minTime,maxTime));
  }
  return(temporary);
}

vector<double> Handler::convert_parameters_forward(vector<double> pars, int number, int maxTime, int minTime){
  vector<double> temporary;
 
  if(number == 3){
    temporary.push_back(logit(pars[0],0.01,0.2));
    temporary.push_back(logit(pars[1],100,15000));
    temporary.push_back(logit(pars[2],minTime,maxTime));
  }
  else if(number == 4){
    temporary.push_back(logit(pars[0],0.0001,0.02));
    temporary.push_back(logit(pars[1],0.01,0.2));
    temporary.push_back(logit(pars[2],100,15000));
    temporary.push_back(logit(pars[3],minTime,maxTime));
  }
  else if(number == 5){
    temporary.push_back(logit(pars[0],0.0001,0.02));
    temporary.push_back(logit(pars[1],0.0001,0.1));
    temporary.push_back(logit(pars[2],0.01,0.2));
    temporary.push_back(logit(pars[3],100,15000));
    temporary.push_back(logit(pars[4],minTime,maxTime));
  }
  else{
    temporary.push_back(logit(pars[0],0.0001,0.02));
    temporary.push_back(logit(pars[1],0.0001,0.1));
    temporary.push_back(logit(pars[2],0.0001,0.1));
    temporary.push_back(logit(pars[3],0.01,0.2));
    temporary.push_back(logit(pars[4],100,15000));
    temporary.push_back(logit(pars[5],minTime,maxTime));
  }
  return(temporary);
}



void Handler::update_options(bool mle, bool useT0, bool useI0, bool _singleEpi, bool savePlot, bool saveResults){
  useMLE = mle;
  optimT0 = useT0;
  optimI0 = useI0;
  singleEpi = _singleEpi;
  plot = savePlot;
  save = saveResults;
}



/* ================================ FINAL LEAST SQUARES ============================= */

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
  double bestAICc, tempAICc, AICc;

  epidemics.push_back(new_epidemic(sir,15));
  epidemicSizes.push_back(4);
  epidemics.push_back(new_epidemic(spike,40));
  epidemicSizes.push_back(3);
  epidemics.push_back(new_epidemic(sir,95));
  epidemicSizes.push_back(4);



  for(unsigned int j = 0; j <= 4; ++j){
    this->temp_data.push_back(this->current_data[j]);
  }
  this->baseModel = this->current_model = base_model(this->temp_data);
  cout << "Here" << endl;
  this->empty_model = create_empty_data_vector(current_data.size());
  epidemics[0]->update_data(current_data);
  //print_vector(baseModel);
  
  plotGraph(current_data,current_data,1);

  finalParams.push_back(0.001);
  finalParams.push_back(0.1);
  finalParams.push_back(500);
  finalParams.push_back(5);
  finalParamsTemp = convert_parameters_forward(finalParams, 4, 15, -5);

  finalParams.clear();
  finalParams.push_back(0.1);
  finalParams.push_back(800);
  finalParams.push_back(30);

  finalParams = convert_parameters_forward(finalParams,3,40,20);
  finalParamsTemp = concatenate_vectors(finalParamsTemp, finalParams);
  
  finalParams.clear();
  finalParams.push_back(0.0008);
  finalParams.push_back(0.1);
  finalParams.push_back(1000);
  finalParams.push_back(85);
  finalParamsTemp = concatenate_vectors(finalParamsTemp, convert_parameters_forward(finalParams, 4, 95, 75));
  printcon(finalParamsTemp);
  finalParams = convert_all_params_back(finalParamsTemp);
  printcon(finalParams);

  residuals = ode_solve(finalParamsTemp);
  plotGraph(residuals,current_data,1);

  cout << calculate_SSE(residuals,current_data) << endl;
  temp_data = current_data;
  cout << "SSE2: " << fitEpidemics(finalParamsTemp) << endl;

  RSquare = 1- fitEpidemics(finalParamsTemp)/SStot(temp_data,1);
  cout << "RSquare before: " << RSquare << endl;
  //  print_vector(get_residuals(residuals,current_data,1));
  vector<double> tempPar;
  Simplex simplex;
  tempPar = simplex.neldermead(&Handler::fitEpidemics, *this,  finalParamsTemp, iterations);
  printcon(tempPar);
  SSE = fitEpidemics(tempPar);
  RSquare = 1- SSE/SStot(temp_data,1);
  cout << "RSquare: " << RSquare << endl;
  residuals = ode_solve(tempPar);
  plotGraph(temp_data,residuals,2);


  return;
 
  printcon(finalParamsTemp);
  finalParams = convert_parameters_back(finalParamsTemp,4,10,-10);
 
  //print_vector(residuals);
  plotGraph(residuals,current_data,1);
  temp_data = current_data;
 
  combinedResults = get_residuals(temp_data,residuals,1);
  //print_vector(combinedResults);
 
  
  

  return;

  finalParams.clear();
  finalParams.push_back(0.0001);
  finalParams.push_back(0.2);
  finalParams.push_back(500);
  finalParams.push_back(10);
  residuals = epidemics[0]->ode_solve_combined(finalParams);
  epidemics[0]->update_data(current_data);
  //print_vector(residuals);
  plotGraph(residuals,residuals,2);


  finalParams.clear();
  finalParams.push_back(0.02);
  finalParams.push_back(0.2);
  finalParams.push_back(500);
  finalParams.push_back(10);
  residuals = epidemics[0]->ode_solve_combined(finalParams);
  epidemics[0]->update_data(current_data);
  //print_vector(residuals);
  plotGraph(residuals,residuals,3);


  finalParams.clear();
  finalParams.push_back(0.02);
  finalParams.push_back(0.01);
  finalParams.push_back(500);
  finalParams.push_back(10);
  residuals = epidemics[0]->ode_solve_combined(finalParams);
  epidemics[0]->update_data(current_data);
  //print_vector(residuals);
  plotGraph(residuals,residuals,4);


  finalParams.clear();
  finalParams.push_back(0.000001);
  finalParams.push_back(0.01);
  finalParams.push_back(500);
  finalParams.push_back(10);
  residuals = epidemics[0]->ode_solve_combined(finalParams);
  epidemics[0]->update_data(current_data);
  //print_vector(residuals);
  plotGraph(residuals,residuals,5);


  return;






















  candidateModels.push_back(spike);
  candidateModels.push_back(sir);
  //candidateModels.push_back(seir);

  // Start of model fitting process
  t1=clock();
  //epidemics.push_back(new_epidemic(sir, 1));
  this->epidemics.push_back(new_epidemic(sir, 10));
  this->epidemicSizes.push_back(4);
  this->epidemics.push_back(new_epidemic(sir, 30));
  this->epidemicSizes.push_back(4);
  /*this->epidemics.push_back(new_epidemic(sir, 50));
  this->epidemicSizes.push_back(4);
  */
  // Start off with a baseline model of the mean of the first 4 points
  for(unsigned int j = 0; j <= 4; ++j){
    this->temp_data.push_back(this->current_data[j]);
  }
  this->baseModel = this->current_model = base_model(this->temp_data);

   // For each time point in the current data set, carry out the fitting procedure
  for(unsigned int i = 20;i<this->current_data.size();++i){
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
	  
	  currentBestSSE = 999999999.9;
	
	  // Calculate the fit without this epidemic
	  tempSSE = optimiseEpidemics(finalParamsTemp, combinedResultsTemp, componentResultsTemp, tempIterations);
	  tempRSquare = 1 - tempSSE/SStot(temp_data,1);
	  tempAICc = aicc(tempSSE, temp_data.size(), finalParamsTemp.size());
	  // If sufficient fit from removing an epidemic, permanently remove epidemic from set
	  // Otherwise, restore epidemics to full number. Also track which removal was best
	  if(tempAICc < AICc){
	    toRemove = this->tempEpidemics[j];
	    cout << "Fit sufficient. Removing an epidemic" << endl;
	    SSE = tempSSE;
	    AICc = tempAICc;
	    currentBestSSE = tempSSE;
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
	currentBestSSE = 999999999.9;
	newEpidemic = candidateModels[f];
	cout << "Considering addition of ";
	print_epidemic_type(newEpidemic);
	cout << endl;
	toRemove = new_epidemic(newEpidemic, i);
	this->epidemics.push_back(toRemove);
	this->epidemicSizes.push_back(toRemove->return_parameters().size());
	cout << "Epidemic sizes: ";
	printcon(epidemicSizes);
	
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
	  currentBestSSE = tempSSE;
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
  
    // Plot graph
    //cout << "Number of epidemics: " << componentResults.size() << endl;
    //cout << "Model size: " << combinedResults.size() << endl;
    //print_vector(combinedResults);
    //print_vector(componentResults[0]);
    //print_vector(componentResults[1]);
    
    plotGraphMulti(componentResults, combinedResults, this->temp_data, i, finalParams, RSquare, 2);    
    
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

double Handler::optimiseEpidemics(vector<double> &parameters, vector<vector<double> > &results, vector<vector<vector<double> > > &allResults, int& itr){
  double SSE = 9999999999999999999999.9;
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
    allResults.push_back(results);
    allResults.push_back(results);
    return(SSE);
  }
  srand(time(NULL));
  // If there are epidemics to be fitted, perform 40 random fits and keep best fitting model
  for(int index=0;index<5;index++){
    
            
    // Clear the temporary seed parameters and seed rand
   
    

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
 
  // Get the combined values from these parameters, as well as a vector of 
  // each sub-epidemic
  results = ode_solve(parameters);
  //  print_vector(results);
  allResults.clear();
  allResults = ode_solve_separate(parameters);
  
  cout << endl;
  return(SSE);

}

vector<double> Handler::convert_all_params_back(vector<double> pars){
    
  vector<double> temporary;
  vector<double> temp2;
  int z = 0;
  temp2.clear();
  for(int i = 0;i<this->epidemicSizes.size();++i){
    temporary.clear();
    for(int f = 0;f<this->epidemicSizes[i];++f){
      temporary.push_back(pars[z]);
      z++;
    }
    temporary = convert_parameters_back(temporary,this->epidemicSizes[i], this->epidemics[i]->return_detection_time(), this->epidemics[i]->return_min_time());
    temp2 = concatenate_vectors(temp2, temporary);
  }
  return(temp2);
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
    for(unsigned int x = 0;x<this->epidemicSizes[i];++x){
      tempParams.push_back(params[z]);
      z++;
    }
    temp_model = epidemics[i]->ode_solve(tempParams);
    current_model = combine_vectors(current_model, temp_model);
  }
  //return(current_model);
  cout << calculate_SSE(current_model,temp_data) << endl;
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
    for(unsigned int x = 0;x<epidemics[i]->return_parameters().size();++x){
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

  sirlimit = mean + 3*sd;
  //  explimit = mean + 6*sd;
 
  SIRresidual = residuals[residuals.size()-1][1];
  //EXPresidual = residuals[residuals.size()-1][1];
  if((residuals[residuals.size()-1][1] > sirlimit)) detected=sir;
    //    cout << "1 detect" << endl;
    /*
    if(residuals[residuals.size()-2][1] > sirlimit){
      cout << "2 detect" << endl;
      if((residuals[residuals.size()-3][1] > sirlimit)){
	cout << "3 detect" << endl;
	detected = sir;
      }
    }
    }*/

  /*

  for(unsigned int i = residuals.size()-1; i > previousResiduals.size();--i){
    if((residuals[i][1] - EXPresidual) < 0) EXPresidual = residuals[i][1];
    if(residuals[i][1] < SIRresidual) SIRresidual = residuals[i][1];
  }
 
  double diffEXP = EXPresidual - explimit;
  double diffSIR = SIRresidual - sirlimit;
  if(diffEXP > 0){ detected = spike;}
  else if(diffSIR > 0){detected = sir;}
  */
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


vector<double> Handler::generate_seed_parameters(){
  vector<double> params;
  vector<double> tempor;
  for(unsigned int i = 0;i<epidemics.size();++i){
    tempor = rand_params(epidemics[i]->return_type());
    tempor.push_back(epidemics[i]->return_seed_time());
    tempor = convert_parameters_forward(tempor, epidemics[i]->return_parameters().size(), epidemics[i]->return_detection_time(),epidemics[i]->return_min_time());
    params = concatenate_vectors(params, tempor);

  }
  return(params);
}


// Generates a vector of random parameters
// THIS IS WHERE WE GENERATE SEED VALUES
vector<double> Handler::rand_params(EpiType _type){
  vector<double> params;
  double beta,gamma,s0,alpha;
  int numb = 0;
  setprecision(9);
  beta = (rand()%100+1)/10000.0;
  alpha = (rand()%100+1)/10000.0;
  gamma = (rand()%200+beta)/1000.0;
  s0 = rand()%2000+1500;
  //t0 = rand()%80+1;
  
  //beta = 0.0005;
  //alpha = 0.001;
  //gamma = 0.05;
  //s0 = 1000;
  //  t0 = 15;
  switch(_type){
  case sir:
    params.push_back((beta));
    params.push_back((gamma));
    params.push_back((s0));
    numb = 4;
    break;
    //params.push_back(log(t0));
    //return(params);
  case seir:
    params.push_back((beta));
    params.push_back((alpha));
    params.push_back((gamma));
    params.push_back((s0));
    numb = 5;
    break;
    //params.push_back(log(t0));
    //return(params);
  case spike:
    params.push_back((gamma));
    params.push_back((s0));
    numb = 3;
    break;
    //params.push_back(log(t0));
  default:
    params.push_back((gamma));
    params.push_back((s0));
    numb = 3;
    //params = convert_parameters_forward(params,numb);
   
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
  string name = "graphs7/output"; // The save location and general name of the graph to be saved
  string _index = to_string((index-1));
  string label, xlab, ylab,xlab2;
  string colours[6] = {"orange","yellow","cyan","red","blue","violet"};
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
  gp << "plot '-' using 1:2 with linespoints lw 1 lt 19 title 'Future Data'";
  gp << ", '-' using 1:2 with linespoints lw 1 lt 19 linecolor rgb \"brown\" title 'Current Data'";
  
  
  // For each sub epidemic, add a dashed line to the graph and a corresponding label to the key
  if(finalResults.size() > 1){
    label = "Baseline Level";
    gp << ", '-' using 1:2 with lines lw 1 lt 0 lc " << 0 << " title '" << label << "'";
  }
  
  if(epidemics.size() > 0){
    for(int f = 0;f<(int)epidemics.size();++f){
    // Create the label depending on the type of epidemic
    label = "Sub-Epidemic " + to_string(f+1);
    label += create_label(epidemics[f], parameters, j);
    gp << ", '-' using 1:2 with lines lw 1 lt 0 linecolor rgb \"" << colours[f] << "\" title '" << label << "'";
    }
  }
  // Secondly, plot the overall model values
  gp << ", '-' using 1:2 with lines lw 1 linecolor rgb \"green\" title 'Total'";
  

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


double Handler::poisson_pmf(const double k, const double lambda) {
  return exp(k * log(lambda) - lgamma(k + 1.0) - lambda);
}






void Handler::overall_test(double targetRsq){
  Epidemic* toRemove1, *toRemove2;
  EpiType newEpidemic;
  vector<double> finalParams, finalParamsTemp;
  vector<vector<double> > baseModel, combinedResults, residuals, combinedResultsTemp, temp;
  vector<vector<vector<double> > > componentResults, componentResultsTemp;
  double SSE;
  int iter;
  newEpidemic = sir;
  
  toRemove1 = new_epidemic(newEpidemic, 0);
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
  for(int i = 0;i<finalParams.size();++i){
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






/* ================================ FINAL LEAST SQUARES ============================= */

void Handler::realtime_fit_single(double targetRSq, EpiType _epi){
  clock_t t1, t2;
  Epidemic *tempEpi,*bestFit;
  vector<EpiType> candidateModels;
  vector<double> finalParams,tempPar;
  vector<vector<double> > combinedResults, tempCombinedResults, residuals; // 
  vector<vector<vector<double> > > tempAll, tempAll1;
  double RSquare,tempRSquare;
  double SSE = 99999999999.9;
  double tempSSE = 99999999999.9;
  double currentBestAICc, tempAICc = 999999999999.9;
  int iterations = 0;
  int tempIterations = 0;
   
  srand(clock());

  candidateModels.push_back(spike);
  candidateModels.push_back(sir);
  candidateModels.push_back(serir);
  candidateModels.push_back(irsir);
  candidateModels.push_back(seir);
 
  // Start off with a baseline model of the mean of the first 4 points
  for(unsigned int j = 0; j <= 4; ++j){
    this->temp_data.push_back(this->current_data[j]);
  }
  this->baseModel = this->current_model = base_model(this->temp_data);


  if(_epi != none) epidemics.push_back(new_epidemic(_epi,10));
 
  finalParams = tempPar = generate_seed_parameters();
  for(unsigned int i = 38;i<this->current_data.size();++i){
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
      cout << "Attempting to find best fitting model" << endl;
      for(unsigned int k = 0; k<candidateModels.size();++k){
	cout << "Considering: ";
	//currentBestSSE = 999999999999.9;
	print_epidemic_type(candidateModels[k]);
	cout << endl;
	tempEpi = NULL;
	tempEpi = new_epidemic(candidateModels[k], 1.0);
	epidemics.push_back(tempEpi);
	tempSSE=optimise_single(tempPar, tempCombinedResults,tempAll1, tempIterations);
	tempRSquare = 1 - tempSSE/SStot(temp_data, 1);
	tempAICc = aicc(tempSSE, temp_data.size(), tempEpi->return_parameters().size());
	cout << "Model AICc: " << tempAICc << endl;
	//if((tempRSquare > targetRSq && tempSSE < currentBestSSE) || (tempRSquare > targetRSq && bestFit != NULL && (tempEpi->return_parameters().size() < bestFit->return_parameters().size())) || (bestFit != NULL && tempEpi->return_parameters().size() == bestFit->return_parameters().size() && tempSSE < currentBestSSE && tempIterations < iterations)){
	if(tempRSquare > targetRSq && tempAICc < currentBestAICc){
	  cout << "Model improvement" << endl;
	  bestFit = new_epidemic(candidateModels[k],1.0);
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
      if(bestFit != NULL){
	cout << "Adding epidemic type: ";
	print_epidemic_type(bestFit->return_type());
	cout << endl;
	epidemics.push_back(bestFit);
	bestFit = NULL;
      }
      reset_epidemics();
    }
    else SSE = optimise_single(finalParams, combinedResults, tempAll, iterations);
    RSquare = 1 - SSE/SStot(temp_data, 1);
    // Transform parameters back to normal space
    for(unsigned int j=0;j<finalParams.size();++j){
      finalParams[j] = exp(finalParams[j]);
    }
    // Plot graph
    //if(plot) plotGraphMulti(tempAll, combinedResults, this->temp_data, i, finalParams, RSquare, 2);    
    
    // Save results
    if(save) cout << "Save results" << endl;
    
    // Print out results for this iteration
    cout << "Final SSE was: " << SSE << endl;
    cout << "Final RSquare: " << RSquare << endl;
    cout << "Final parameters: ";
    printcon(finalParams) ;
    cout << endl << "-----------------" << endl;
    return;
  }
  cout << "Fitting procedure complete" << endl;
}


double Handler::optimise_single(vector<double> &parameters, vector<vector<double> > &results, vector<vector<vector<double> > > &allResults, int& itr){
  double SSE = 9999999999.9;
  double tempSSE;
  vector<double> tempPar, seedParams;
  Simplex simplex;

  parameters.clear();
  allResults.clear();

  for(int index = 0;index<5;index++){
    seedParams.clear();
    seedParams = generate_seed_parameters();
    cout << "Pre optimisation SSE: " << fitEpidemics(seedParams) << endl;
    if(this->useMLE == true) tempPar = simplex.neldermead(&Handler::fitEpidemicsMLE, *this,  seedParams, itr);
    else tempPar = simplex.neldermead(&Handler::fitEpidemics, *this,  seedParams, itr);
    // Store the SSE value for this
    tempSSE = fitEpidemics(tempPar);
    cout << "Iterations: " << itr << endl;
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