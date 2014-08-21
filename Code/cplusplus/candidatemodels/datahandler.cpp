using namespace std;

#include "datahandler.hpp"
#include "epidemic.hpp"
#include "sir.hpp"
#include "seir.hpp"
#include "spike.hpp"
#include "simplex.hpp"
#include <string>

const float EPSILON = 0.1;

Handler::Handler(){
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





/* ================================ FINAL LEAST SQUARES ============================= */

void Handler::realtime_fit_multi(double targetRSq){
  clock_t t1, t2; // To record the total run time
  EpiType newEpidemic;
  Epidemic* toRemove;
  vector<double> finalParams, finalParamsTemp;
  vector<vector<double> > baseModel, combinedResults, residuals, combinedResultsTemp;
  vector<vector<vector<double> > > componentResults, componentResultsTemp;
  double RSquare, RSquareK1, SSE, tempSSE, currentBestSSE;

  // Start of model fitting process
  t1=clock();
  
  // Test starting with 2 epidemics
  //newEpidemic = sir;

  //add_epidemic(newEpidemic, 10);
  //newEpidemic = spike;
  //add_epidemic(newEpidemic, 40);

  // Start off with a baseline model of the mean of the first 4 points
  for(unsigned int j = 0; j <= 4; ++j){
    temp_data.push_back(current_data[j]);
  }
  current_model = base_model(temp_data);
  baseModel = current_model;

  // For each time point in the current data set, carry out the fitting procedure
  for(unsigned int i = 5;i<current_data.size();++i){
    cout << endl << "-----------------" << endl;
    cout << "Iteration number " << i << endl;
    
    // Clear temp vector and temp epidemics
    temp_data.clear();
    tempEpidemics.clear();
    toRemove = NULL;

    // Store all data up to the current index and update epidemics
    for(unsigned int j = 0; j < i; ++j){
      temp_data.push_back(current_data[j]);
    } 
    for(unsigned int f = 0;f<epidemics.size();++f){
      epidemics[f]->update_data(temp_data);
    }
    cout << "About to optimise" << endl;
    // Create an optimised model fit for the given data (temp) and calculate RSquare from the SSE
    SSE = optimiseEpidemics(finalParams, combinedResults, componentResults);
    RSquare = 1 - SSE/SStot(temp_data, 1);
    
    // If fit is sufficiently good with k epidemics, try fitting with k-1 epidemics. If this fit
    // is sufficient, use k-1 epidemics.
    if(RSquare > targetRSq && epidemics.size() > 0){
      cout << "Considering removal of an epidemic" << endl;
      // Store current set of epidemics
      tempEpidemics = epidemics;
      
      // Try removing each epidemic so far
      for(unsigned int j = 0;j<epidemics.size();++j){
	// Create vector of epidemic types that have been removed
	vector<EpiType> removedEpidemics;
	
	// Try removing the j'th epidemic. If that type has already been removed, increase j.
	// if epidemics[j]->return_type() is not in removedEpidemics, then that. Otherwise, j++ and try again.
	if(!check_already_tested(removedEpidemics, tempEpidemics[j]->return_type())){
	  epidemics = fewer_epidemics(j);
	  currentBestSSE = 999999999.9;
	
	  // Calculate the fit without this epidemic
	  tempSSE = optimiseEpidemics(finalParamsTemp, combinedResultsTemp, componentResultsTemp);
	  RSquareK1 = 1 - tempSSE/SStot(temp_data,1);
	
	  // If sufficient fit from removing an epidemic, permanently remove epidemic from set
	  // Otherwise, restore epidemics to full number. Also track which removal was best
	  if(RSquareK1 > targetRSq && tempSSE < currentBestSSE){
	    toRemove = this->tempEpidemics[j];
	    cout << "Fit improved. Removing an epidemic" << endl;
	    SSE = tempSSE;
	    currentBestSSE = tempSSE;
	    RSquare = RSquareK1;
	    finalParams = finalParamsTemp;
	    combinedResults = combinedResultsTemp;
	    componentResults = componentResultsTemp;
	  }
	  else{
	    cout << "No need to remove" << endl;
	    toRemove = NULL;
	    epidemics = tempEpidemics;
	  }
	}
      }
      if(toRemove != NULL) remove_epidemic(toRemove);
    }

    // Reset pointers in case epidemic was deleted
    reset_epidemics();
    cout << "Epidemics reset " << endl;
    // Check if a new epidemic has started. If so, and RSquare from current K epidemics 
    // insufficient, try fitting another. If this is a significantly better fit, permanently
    // add an extra epidemic.
    residuals = get_residuals(temp_data, current_model, 1);
    newEpidemic = check_epidemic(residuals);
    if(newEpidemic != none && RSquare < targetRSq){
      cout << "Epidemic detected!" << endl;
      cout << "Epidemic type: " << newEpidemic << endl;
      // Temporarily store the current epidemic state and add a new epidemic to the current set
      tempEpidemics = epidemics;
      toRemove = add_epidemic(newEpidemic, i);
      cout << "Number of epidemics: " << epidemics.size() << endl;
      tempSSE = optimiseEpidemics(finalParamsTemp, combinedResultsTemp, componentResultsTemp);
      RSquareK1 = 1 - tempSSE/SStot(temp_data,1);
      cout << "Fit with additional epidemic: " << RSquareK1 << endl;
      if(RSquareK1 > RSquare){
	cout << "Fit improved with additional epidemic. Adding epidemic" << endl;
	SSE = tempSSE;
	currentBestSSE = tempSSE;
	RSquare = RSquareK1;
	finalParams = finalParamsTemp;
	combinedResults = combinedResultsTemp;
	componentResults = componentResultsTemp;
      }
      else{
	//cout << "Deleting!" << endl;
	//cout << "Address1: " << toRemove << endl;
	remove_epidemic(toRemove);
      }
    }
    // Reset pointers in case epidemic was deleted;
    reset_epidemics();
       
    // Transform parameters back to normal space
    for(unsigned int j=0;j<finalParams.size();++j){
      finalParams[j] = exp(finalParams[j]);
    }
    // Plot graph
    plotGraphMulti(componentResults, combinedResults, temp_data, i, finalParams, RSquare, 2);    
    
    // Print out results for this iteration
    cout << "Final SSE was: " << SSE << endl;
    cout << "Final RSquare: " << RSquare << endl;
    cout << "Final parameters: ";
    printcon(finalParams) ;
    cout << endl << "-----------------" << endl;
    current_model = combinedResults;
  }
  t2=clock();
  cout << "Time taken: " << (t2-t1)/CLOCKS_PER_SEC << endl << endl;
}


bool Handler::check_already_tested(vector<EpiType> tested, EpiType toCheck){
  for(unsigned int p = 0; p <tested.size();++p){
    if(toCheck == tested[p]){
      cout << "Already tried removing epidemic type. Skip." << endl;
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
  cout << "Reseting epidemics" << endl;
  vector<Epidemic*> epi;
  cout << "Number of epidemics: " << epidemics.size() << endl;
  for(unsigned int i = 0;i<epidemics.size();++i){
    if(epidemics[i] != NULL){
      cout << "Pushing back epidemic" << endl;
      cout << "Address: " << epidemics[i] << endl;
      cout << epidemics[i]->return_type() << endl;
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
Epidemic* Handler::add_epidemic(EpiType _newEpidemic, int time){
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
  epidemics.push_back(additionalEpi);
  return(additionalEpi);
}

double Handler::optimiseEpidemics(vector<double> &parameters, vector<vector<double> > &results, vector<vector<vector<double> > > &allResults){
  double SSE = 99999999999.9;
  double tempSSE;
  vector<double> tempPar, seedParams;
  vector<vector<vector<double> > > tempAll;
  Simplex simplex;
  
  // If no epidemics have yet been detected, use the mean of the current data as the 
  // current model and return the corresponding SSE.
  if(epidemics.size() == 0){
    results = base_model(temp_data);
    SSE = calculate_SSE(results, temp_data);
    parameters.clear();
    allResults.clear();
    allResults.push_back(results);
    return(SSE);
  }

  // If there are epidemics to be fitted, perform 40 random fits and keep best fitting model
  for(int index=0;index<10;index++){
    
            
    // Clear the temporary seed parameters and seed rand
   
    srand(clock());

    // Create a list of random seed parameters
    seedParams.clear();      
    seedParams = generate_seed_parameters();
    
    // Get the optimised parameters from nelder mead algorithm
    tempPar = simplex.neldermead(&Handler::fitEpidemics, *this,  seedParams);
    // Store the SSE value for this
    tempSSE = fitEpidemics(tempPar);
	
    // If this SSE value is better than the previous, store it and the
    // corresponding parameters
    if(tempSSE < SSE){
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
  allResults.push_back(base_model(temp_data));
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
  return(calculate_SSE(current_model, temp_data));
}
  
// As above, but returns the overall model.
vector<vector<double> > Handler::ode_solve(vector<double> params){
  vector<vector<double> > overallResults = empty_model;
  int z = 0;
  for(unsigned int i = 0;i<epidemics.size();++i){
    tempParams.clear();
    temp_model.clear();
    for(unsigned int x = 0;x<epidemics[i]->return_parameters().size();++x){
      tempParams.push_back(params[z]);
      z++;
    }
    temp_model = epidemics[i]->ode_solve_combined(tempParams);
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
    temp_model = epidemics[i]->ode_solve_combined(tempParams);
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

  cout << "Latest residual value: " << residuals[residuals.size()-1][1] << endl;

  // Calculate the SD and mean of these residuals
  mean = calculate_mean(previousResiduals, 1);
  sd = calculate_sd(previousResiduals, 1);

  cout << "Current residual mean: " << mean << endl;
  cout << "Current residual SD: " << sd << endl;

  sirlimit = mean + 2*sd;
  explimit = mean + 6*sd;
  
  cout << "SIR limit: " << sirlimit << endl;
  cout << "EXP limit: " << explimit << endl;
  
  SIRresidual = residuals[residuals.size()-1][1];
  EXPresidual = residuals[residuals.size()-1][1];

  cout << "SIR residual: " << SIRresidual << endl;
  for(unsigned int i = residuals.size()-1; i > previousResiduals.size();--i){
    if(residuals[i][1] < EXPresidual) EXPresidual = residuals[i][1];
    //if(residuals[i][1] < SIRresidual) SIRresidual = residuals[i][1];
  }

  if(EXPresidual > explimit) detected = spike;
  else if(SIRresidual > sirlimit) cout << "Detected SIR" << endl; detected = sir;
  
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




vector<double> Handler::generate_seed_parameters(){
  vector<double> params;
  for(unsigned int i = 0;i<epidemics.size();++i){
    params = concatenate_vectors(params, rand_params(epidemics[i]->return_type()));
    params.push_back(log(epidemics[i]->return_detection_time()));
  }
  return(params);
}


// Generates a vector of random parameters
// THIS IS WHERE WE GENERATE SEED VALUES
vector<double> Handler::rand_params(EpiType _type){
  vector<double> params;
  double beta,gamma,s0,alpha;
  
  setprecision(9);
  beta = (rand()%100+1)/10000.0;
  alpha = (rand()%100+1)/10000.0;
  gamma = (rand()%200+beta)/1000.0;
  s0 = rand()%1000+100;
  //t0 = rand()%80+1;
  
  /*beta = 0.001;
  alpha = 0.001;
  gamma = 0.1;
  s0 = 500;
  t0 = 15;*/
  switch(_type){
  case sir:
    params.push_back(log(beta));
    params.push_back(log(gamma));
    params.push_back(log(s0));
    //params.push_back(log(t0));
    return(params);
  case seir:
    params.push_back(log(beta));
    params.push_back(log(alpha));
    params.push_back(log(gamma));
    params.push_back(log(s0));
    //params.push_back(log(t0));
    return(params);
  case spike:
    params.push_back(log(gamma));
    params.push_back(log(s0));
    //params.push_back(log(t0));
  default:
    return(params);
  }
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
  string name = "graphs/output"; // The save location and general name of the graph to be saved
  string _index = to_string(index);
  string label, xlab;
  int j = 0;

  xlab = "RSquare: " + to_string(_RSquare);  // Xlabel is the RSquare value
  name = name + _index + ".jpeg";  // Graph name

  gp << "set terminal jpeg size 1000,800 enhanced font \"Helvetica, 10\"\n";    // Here edits type, size and formatting of graph
  gp << "set xlabel '" << xlab << "'\n";
  gp << "set output '" << name << "'\n";  // Set output to the specified file name
  gp << "set termoption dash\n"; // Allow dashes

  // Firstly, plot the actual data (3rd column)
  gp << "plot '-' using 1:2 with linespoints lt 19 title 'Data'";
  
  // Secondly, plot the overall model values
  gp << ", '-' using 1:2 with lines lw 1 title 'Total'";
  
  // For each sub epidemic, add a dashed line to the graph and a corresponding label to the key
  label = "Baseline Level";
  gp << ", '-' using 1:2 with lines lt 0 lc " << 0 << " title '" << label << "'";
 

  for(int f = 0;f<(int)epidemics.size();++f){
    // Create the label depending on the type of epidemic
    label = "Sub-Epidemic " + to_string(f+1);
    label += create_label(epidemics[f], parameters, j);
    gp << ", '-' using 1:2 with lines lt 0 lc " << f << " title '" << label << "'";
  }


  gp << "\n";
  
  // Send the data to Gnuplot
  gp.send1d(data);
  gp.send1d(totalResults);
  for(unsigned int j = 0;j<finalResults.size();++j){
    gp.send1d(finalResults[j]);
  }
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






void Handler::likelihood_test(vector<double> &params){
  vector<double> parameters, tempParams;
  vector<vector<double> > data, model;
  //print_vector(current_data);
  SIR mySIR(double(current_data.size()), current_data, sir, 1);
  Simplex simplex;
  parameters.push_back(log(0.001));
  parameters.push_back(log(0.1));
  parameters.push_back(log(500));
  parameters.push_back(log(10));
  
  parameters.push_back(log(0.001));
  parameters.push_back(log(0.05));
  parameters.push_back(log(650));
  parameters.push_back(log(50.0));
  
  parameters.push_back(log(0.002));
  parameters.push_back(log(0.075));
  parameters.push_back(log(450));
  parameters.push_back(log(75.0));
  
 
  data = mySIR.ode_solve_combined(parameters);
 
  mySIR.update_data(data);

  parameters.clear();
  parameters.push_back(log(0.002));
  parameters.push_back(log(0.11));
  parameters.push_back(log(500));
  parameters.push_back(log(10));

  parameters.push_back(log(0.0011));
  parameters.push_back(log(0.055));
  parameters.push_back(log(670));
  parameters.push_back(log(50));
  parameters.push_back(log(0.002));
  parameters.push_back(log(0.075));
  parameters.push_back(log(450));
  parameters.push_back(log(75.0));
  
  cout << "Test: " << mySIR.mle_sir(parameters) << endl;
  
  tempParams = simplex.neldermead(&SIR::mle_sir, mySIR,  parameters);
  cout << "Final log likelihood: " << mySIR.mle_sir(tempParams) << endl;
  parameters.clear();
  parameters.push_back(log(0.001));
  parameters.push_back(log(0.1));
  parameters.push_back(log(500));
  parameters.push_back(log(10));
  parameters.push_back(log(0.001));
  parameters.push_back(log(0.05));
  parameters.push_back(log(650));
  parameters.push_back(log(50.0));
  parameters.push_back(log(0.002));
  parameters.push_back(log(0.075));
  parameters.push_back(log(450));
  parameters.push_back(log(75.0));
  
  cout << "Likelihood of actual parameters: " << mySIR.mle_sir(parameters) << endl;
  //cout << "SSE: " << mySIR.calculate_SSE(model,data) << endl;
  model = mySIR.ode_solve_combined(tempParams);
  for(unsigned int i = 0;i<tempParams.size();++i){
    tempParams[i] = exp(tempParams[i]);
  }
  printcon(tempParams);
  plotGraph(model, data, 1);
    
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







void Handler::overall_test(double targetRsq){
  EpiType newEpidemic;
  vector<double> finalParams, finalParamsTemp;
  vector<vector<double> > baseModel, combinedResults, residuals, combinedResultsTemp, temp;
  vector<vector<vector<double> > > componentResults, componentResultsTemp;
  double RSquare, RSquareK1, SSE, tempSSE, currentBestSSE;
  cout << "Here" << endl;
  newEpidemic = sir;
  add_epidemic(newEpidemic, 0);
  for(unsigned int f = 0;f<epidemics.size();++f){
    epidemics[f]->update_data(current_data);
  }
    
  cout << "Number of epidemics: " << epidemics.size() << endl;
  cout << "Data size: " << current_data.size() << endl;

  temp_data = current_data;
  finalParams.push_back(log(0.001));  
  finalParams.push_back(log(0.1));  
  //finalParams.push_back(log(0.05));  
  finalParams.push_back(log(500));
  finalParams.push_back(log(1));
    //temp = epidemics[0]->ode_solve(finalParams);
  temp=ode_solve(finalParams);
  cout << "Temp size: " << temp.size() << endl;
  //  print_vector(temp);
  print_vector(temp_data);
  //return;
  cout << calculate_SSE(temp,temp_data) << endl;
  SSE = optimiseEpidemics(finalParams, combinedResults, componentResults);
  
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

