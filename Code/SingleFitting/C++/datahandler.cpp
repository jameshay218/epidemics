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

// Amount of allowable difference between data points. NOTE THAT THIS MUST MATCH THE RUNGE-KUTTA STEP SIZE IN EPIDEMIC.CPP
const float EPSILON = 0.1;

// Default constructor
Handler::Handler(){
  useMLE = false;
  optimT0 = false;
  optimI0 = false;
  singleEpi = false;
  plot = true;
  save = false;
}

// Default deconstructor. Ensures that all epidemics are deleted correctly
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


// Function to allow user to update optimisation options.
void Handler::update_options(bool mle, bool useT0, bool useI0, bool _singleEpi, bool savePlot, bool saveResults, string location){
  useMLE = mle; // Use MLE based fitting
  optimT0 = useT0; // Include T0 in optimisation
  optimI0 = useI0; // Include I0 in optimisation
  singleEpi = _singleEpi; // Single epidemic model
  plot = savePlot; // Save graphs to folder
  save = saveResults; // Save results to folder
  saveLocation = location; // Save location
}

// Ensures that the values provided by a given set of model data within a sensible range.
// If not, replaces the nonsense value with a copy of the previous value.
void Handler::sense_check(vector<vector<double> > &model){
  if(model[0][1] > 5000 || model[0][1] <0 || (model[0][1] != model[0][1])) model[0][1] = 0;
  for(unsigned int j = 1;j<model.size();++j){
    if(model[j][1] > 5000 || model[j][1] <0 || (model[j][1] != model[j][1])){
      model[j][1] = model[j-1][1];
    }
  }
}



/* ================================ FINAL LEAST SQUARES ============================= */
/* ------------------------ Main single epidemic fitting framework ------------------ */
/* Takes a target R Square value which acts as the threshold for epidemic addition/
   removal. Also takes an epidemic type to be optimised. If this type is `none', then
   the unknown model fitting procedure is exected.

*/


void Handler::realtime_fit_single(double targetRSq, EpiType _epi){
  clock_t t1, t2, t3;           // Record time taken
  Epidemic *tempEpi,*bestFit;   // Temporary Epidemic pointers to track addition and removal of epidemics
  vector<EpiType> candidateModels;  // A list of candidate model types to be considered
  vector<double> finalParams,tempPar;   // Stores the current best and temporary overall model parameters
  vector<EpiType> selectedTypes;        // Tracks the epidemic types that have already been tested
  vector<vector<double> > combinedResults, tempCombinedResults, residuals, saveResults; // Stores the overall epidemic model, residuals, and results of the fitting procedure to be saved
  vector<vector<vector<double> > > tempAll, tempAll1;  // Stores each sub epidemic model separately
  double RSquare,tempRSquare, RMSE; // Statistics
  double SSE = 99999999999.9;  
  double tempSSE = 99999999999.9;
  double currentBestAICc, AICc, tempAICc = 999999999999.9;
  int iterations = 0; // Records the number of iterations taken by the Nelder Mead algorithm
  int tempIterations = 0; // Temporary iterations store

  // Start timer and seed random
  t1=clock();
  srand(time(NULL));

  if(useMLE) cout << "Using MLE" << endl;


  /* Add the candidate models that we wish to consider in the fitting procedure */
  // candidateModels.push_back(spike);
  candidateModels.push_back(sir);
  //candidateModels.push_back(serir);
  candidateModels.push_back(irsir);
  candidateModels.push_back(seir);
 
  // Start off with a baseline model of the mean of the first 4 points
  finalParams = generate_seed_parameters();
  for(unsigned int j = 0; j <= 4; ++j){
    this->temp_data.push_back(this->current_data[j]);
  }
  this->baseModel = this->current_model = base_model(this->temp_data);

  // If an epidemic type was specified, add this to the list of epidemics to be optimised
  if(_epi != none) epidemics.push_back(new_epidemic(_epi,1,1.0));
  cout << "Considering: ";
  print_epidemic_type(_epi);

  finalParams = tempPar = generate_seed_parameters();

  /* Begin the iterative fitting procedure, starting at a given index and proceeding until the end of the data set */
  for(unsigned int i = 30;i<this->current_data.size();++i){
    cout << endl << "-----------------" << endl;
    cout << "Iteration number " << i << endl;
    
    // Clear temp vector and temp epidemics
    this->temp_data.clear();

    // Store all data up to the current index and update epidemics with data
    for(unsigned int j = 0; j < i; ++j){
      this->temp_data.push_back(this->current_data[j]);
    } 

    /* If no epidemic type was specified, then begin the procedure to select the best epidemic model
       based on AICc score. Repeat every tenth data point */
    if(_epi == none && (epidemics.size() == 0 || i%10 == 0)){
      // Clear epidemic vector
      if(epidemics.size() != 0) remove_epidemic(epidemics[0]); reset_epidemics();
      bestFit = NULL;
      currentBestAICc = 99999999999999999999.9;
      cout << "Attempting to find best fitting model" << endl;
      cout << "Starting with: " << epidemics.size() << " epidemics" << endl;
      
      // Go through the list of candidate models
      for(unsigned int k = 0; k<candidateModels.size();++k){
	cout << "Considering: ";
	print_epidemic_type(candidateModels[k]);
	cout << endl;
	tempEpi = NULL;
	
	// Create a new epidemic to be considered
	tempEpi = new_epidemic(candidateModels[k],1,1.0);
	epidemics.push_back(tempEpi);

	//
	/* CALL TO OPTIMISATION PROCEDURE */
	tempSSE=optimise_single(tempPar, tempCombinedResults,tempAll1, tempIterations);
	//
	//

	if(useMLE) {
	  cout << "Negative log likelihood: " << tempSSE << endl; 
	  tempSSE = fitEpidemics(tempPar);
	}

	// Calculate R Square and AICc
	tempRSquare = 1 - tempSSE/SStot(temp_data, 1);
	cout << "Parameters: " << tempEpi->return_parameters().size() << endl;
	cout << "Temp SSE: " << tempSSE << endl;
	cout << "Temp RSquare: " << tempRSquare << endl;
	tempAICc = aicc(tempSSE, temp_data.size(), tempEpi->return_parameters().size());
	cout << "Model AICc: " << tempAICc << endl;

	// If R Square is sufficient and AICc is better than the current best, then record this epidemic
	// to be added
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
	// Memory clean up of temporary epidemics
	remove_epidemic(tempEpi);
	reset_epidemics();
      }
      cout << "Fitting complete" << endl;

      // If a best fit epidemic was found, add this to the list of epidemics
      if(bestFit != NULL){
	cout << "Adding epidemic type: ";
	print_epidemic_type(bestFit->return_type());
	selectedTypes.push_back(bestFit->return_type());
	cout << endl;
	epidemics.push_back(bestFit);
	bestFit = NULL;
      }
      // Clean up pointers
      reset_epidemics();
    }
    else SSE = optimise_single(finalParams, combinedResults, tempAll, iterations); // If type was specified, call normal optimisation procedure

    if(useMLE) {cout << "Negative log likelihood: " << tempSSE << endl; tempSSE = fitEpidemics(tempPar); }

    // If no good model fit was found, then skip this data point
    if(epidemics.size() == 0){ cout << "No decent fit. Skipping to next data point." << endl; continue;}

    // Calculate and save final results
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

  // At end of procedure, record total run time and save all results to a .csv file
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

/* Optimisation procedure. Takes a reference to parameters, combined results, component results and number of iterations taken. The function updates these
   arguments with the results of the Nelder-Mead optimisation procedure. */
double Handler::optimise_single(vector<double> &parameters, vector<vector<double> > &results, vector<vector<vector<double> > > &allResults, int& itr){
  double SSE = 99999999999999.9;
  double tempSSE;
  vector<double> tempPar, seedParams;
  Simplex simplex; // Create an instance of the Nelder-Mead simplex class.

  parameters.clear();
  allResults.clear();
  parameters = generate_seed_parameters();

  // Find the best of twenty independent optimisation runs
  for(int index = 0;index<20;index++){
    seedParams.clear();
    seedParams = generate_seed_parameters();

    if(useMLE == true){ tempPar = simplex.neldermead(&Handler::fitEpidemicsMLE, *this,  seedParams, itr);} // If using MLE, call MLE based objective function NOT WORKING
    else  tempPar = simplex.neldermead(&Handler::fitEpidemics, *this,  seedParams, itr); // Otherwise, call the normal SSE based objective function

    // Store the SSE value for this
    if(useMLE == true) tempSSE = fitEpidemicsMLE(tempPar);
    else tempSSE = fitEpidemics(tempPar);

    // If recent SSE was improvement on previous best, store these parameters
    if(tempSSE < SSE){
      SSE = tempSSE;
      parameters=tempPar;
    }
    cout << "." << flush;
  }
  // Store the results of solving these parameters
  results = ode_solve(parameters);
  allResults = ode_solve_separate(parameters);

  cout << endl;
  return SSE;
}

// Checks if the given Epidemic type is within the provided list of EpiTypes
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

// Resets the pointers of the current epidemic vector to remove any NULL pointers
void Handler::reset_epidemics(){
  vector<Epidemic*> epi;
  for(unsigned int i = 0;i<epidemics.size();++i){
    if(epidemics[i] != NULL){
      epi.push_back(epidemics[i]);
    }
   }
  epidemics = epi;
}

// Finds the given epidemic pointer in the list of epidemics and removes it
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



/* ========================================= MODEL FITTING PROCEDURES ==============================*/
// Goes through each epidemic in the current list and does the following: 1. Store the correct number
// of parameters from the provided vector; 2. Solve the current epidemic model with those parameters
// 3. Add the contribution of that epidemic to the overall model. 4. Calculate the SSE between the
// current model and current data.  
double Handler::fitEpidemics(vector<double> params){
  int z = 0;
  current_model = create_empty_data_vector(temp_data.size());
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

// As above, but minimises the dpois value
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

 // As above, but returns the overall model rather than the SSE.
 vector<vector<double> > Handler::ode_solve(vector<double> params){
   vector<vector<double> > overallResults = create_empty_data_vector(temp_data.size());
   int z = 0;
   for(unsigned int i = 0;i<epidemics.size();++i){
     tempParams.clear();
     temp_model.clear();
     for(unsigned int x = 0;x<epidemics[i]->return_parameters().size();++x){
       tempParams.push_back(params[z]);
       z++;
     }
     temp_model = epidemics[i]->ode_solve(tempParams);
     sense_check(temp_model); // Checks that the solved model is within a reasonable range
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
     sense_check(temp_model);
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

   // Get line of data and split using commas. Push this to the data vector
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


// Creates an empty vector of data of the given size, where the time goes from 0 to _rows.
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

// Calculates the AICc score for a given SSE, N and number of parameters, K
double Handler::aicc(double sse, int n, int k){
  if(n == 0 || (n-k-1) == 0) return 99999999999.9;
  double aic = n*log(sse/n) + 2*k;
  return(aic + (2*k*(k+1))/(n-k-1));
}

// Old function to check if the provided parameters are outside of a sensible range. Unique checks for 
// each epidemic type.
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
  
// Generates a vector of log transformed seed parameters for each ongoing epidemic
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
 
  s0 = rand()%5000+500;
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
    beta = (rand()%400+40)/10000.0;
    gamma = (rand()%200+50)/10000.0;
    alpha = (rand()%500+50)/10000.0;
    params.push_back(log(beta));
    params.push_back(log(alpha));
    params.push_back(log(gamma));
    params.push_back(log(s0));
    break;
    return(params);
  case serir:
    beta = (rand()%100+50)/10000.0;
    gamma = (rand()%200+20)/10000.0;
    alpha = (rand()%100+50)/10000.0;
    mu = (rand()%200+20)/10000.0;
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


// Takes two X 2D vectors. Returns a new two dimensional vector where the second column of both provided vectors
// are summed at the same time point. Note that the first vector, model, should be the shorter of the two vectors.
// The idea behind this function is to sum the infected individuals from two data vectors where the time steps are not
// necessarily the same
vector<vector<double> > Handler::combine_vectors(vector<vector<double> > model, vector<vector<double> > data){
  unsigned int i = 0;
  unsigned int j = 0;
  unsigned int k = 0;
  vector<vector<double> > final = model;
  
  while(i < model.size()){
    // If time points the same, then combine and store
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
      else logLikelihood += poisson_pmf(small[i][1], big[j][0]);
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

// Poisson pmf function
double Handler::poisson_pmf(double k, double lambda) {
  return (k * log(lambda) - lgamma(k + 1.0) - lambda);
}





