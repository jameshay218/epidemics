using namespace std;

#include "datahandler.hpp"
#include "epidemic.hpp"
#include "sir.hpp"
#include "simplex.hpp"
#include <string>

Handler::Handler(){
}

Handler::~Handler(){
  cout << "Delete Data Handler" << endl;
}

/* ================================ FINAL LEAST SQUARES ============================= */


void Handler::realtime_fit_multi(double targetRSq){
  clock_t t1, t2; // To record the total run time
  EpiType newEpidemic;
  vector<double> finalParams, finalParamsTemp;
  vector<vector<double> > temp, baseModel, combinedResults, residuals, combinedResultsTemp;
  vector<vector<vector<double> > > componentResults, componentResultsTemp;
  double RSquare, RSquareK1, SSE, tempSSE, currentBestSSE;
  // Start of model fitting process
  t1=clock();

  // Start off with a baseline model of the mean of the first 4 points
  for(unsigned int j = 0; j <= 4; ++j){
    temp.push_back(current_data[j]);
  }
  current_model = base_model(temp);
  baseModel = current_model;

  // For each time point in the current data set, carry out the fitting procedure
  for(unsigned int i = 50;i<current_data.size();++i){
    cout << endl << "-----------------" << endl;
    cout << "Iteration number " << i << endl;
    
    // Clear temp vector
    temp.clear();
   
    // Store all data up to the current index
    for(unsigned int j = 0; j < i; ++j){
      temp.push_back(current_data[j]);
    } 
    
    
    // Create an optimised model fit for the given data (temp) and calculate RSquare from the SSE
    SSE = optimiseEpidemics(finalParams, temp, combinedResults, componentResults);
    RSquare = 1 - SSE/SStot(temp, 2);
    
    // If fit is sufficiently good with k epidemics, try fitting with k-1 epidemics. If this fit
    // is sufficient, use k-1 epidemics.
    if(RSquare > targetRSq && epidemics.size() > 0){
      cout << "Considering removal of an epidemic" << endl;
      // Store current set of epidemics
      tempEpidemics = epidemics;
      
      // Try removing each epidemic so far
      for(unsigned int j = 0;j<epidemics.size();++j){
	// Try removing the j'th epidemic
	epidemics = fewer_epidemics(j);
	currentBestSSE = 999999999.9;
	
	// Calculate the fit without this epidemic
	tempSSE = optimiseEpidemics(finalParamsTemp, temp, combinedResultsTemp, componentResultsTemp);
	RSquareK1 = 1 - tempSSE/SStot(temp,2);
	
	// If sufficient fit from removing an epidemic, permanently remove epidemic from set
	// Otherwise, restore epidemics to full number. Also track which removal was best
	if(RSquareK1 > targetRSq && tempSSE < currentBestSSE){
	  cout << "Fit improved. Removing an epidemic" << endl;
	  SSE = tempSSE;
	  currentBestSSE = tempSSE;
	  RSquare = RSquareK1;
	  finalParams = finalParamsTemp;
	  combinedResults = combinedResultsTemp;
	  componentResults = componentResultsTemp;
	}
	else epidemics = tempEpidemics;
      }
    }

    // Check if a new epidemic has started. If so, and RSquare from current K epidemics 
    // insufficient, try fitting another. If this is a significantly better fit, permanently
    // add an extra epidemic.
    residuals = get_residuals(temp, current_model, 2);
    newEpidemic = check_epidemic(residuals);
    if(newEpidemic != none && RSquare < targetRSq){
      cout << "Epidemic detected!" << endl;

      // Temporarily store the current epidemic state and add a new epidemic to the current set
      tempEpidemics = epidemics;
      add_epidemic(newEpidemic);

      tempSSE = optimiseEpidemics(finalParamsTemp, temp, combinedResultsTemp, componentResultsTemp);
      RSquareK1 = 1 - tempSSE/SStot(temp,2);
      if(RSquareK1 > targetRSq){
	cout << "Fit improved with additional epidemic. Adding epidemic" << endl;
	SSE = tempSSE;
	currentBestSSE = tempSSE;
	RSquare = RSquareK1;
	finalParams = finalParamsTemp;
	combinedResults = combinedResultsTemp;
	componentResults = componentResultsTemp;
      }
      else epidemics = tempEpidemics;
    }
   

    // Transform parameters back to normal space
    for(unsigned int j=0;j<finalParams.size();++j){
      finalParams[j] = exp(finalParams[j]);
    }
    // Plot graph
    plotGraphMulti(componentResults, combinedResults, temp, i, finalParams, RSquare, 2);    
    
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

vector<Epidemic> Handler::fewer_epidemics(int j){
  vector<Epidemic> epi;
  for(unsigned int i = 0;i<tempEpidemics.size();++i){
    if((int)i != j){
      epi.push_back(tempEpidemics[i]);
    }
  }
  return(epi);
}


void Handler::add_epidemic(EpiType _newEpidemic){
  Epidemic *additionalEpi;
  switch(_newEpidemic){
  case sir:
    additionalEpi = new SIR(current_data.size(),current_data,_newEpidemic);
    epidemics.push_back(*additionalEpi);
    break;
  default:
    additionalEpi = new SIR(current_data.size(),current_data,_newEpidemic);
    epidemics.push_back(*additionalEpi);
    break;
  }
}

double Handler::optimiseEpidemics(vector<double> &parameters, vector<vector<double> > data, vector<vector<double> > &results, vector<vector<vector<double> > > &allResults){
  double SSE = 99999999999.9;
  double tempSSE;
  vector<double> tempParams, seedParams;
  vector<vector<vector<double> > > tempAll;
  Simplex simplex;
  
  // If no epidemics have yet been detected, use the mean of the current data as the 
  // current model and return the corresponding SSE.
  if(epidemics.size() == 0){
    results = base_model(data);
    SSE = Epidemic::calculate_SSE(results, data, 2);
    parameters.clear();
    allResults.clear();
    allResults.push_back(results);
    return(SSE);
  }

  // If there are epidemics to be fitted, perform 40 random fits and keep best fitting model
  for(int index=0;index<20;index++){
    cout << "." << flush;
    
    // Clear the temporary seed parameters and seed rand
    seedParams.clear();      
    srand(clock());

    // Create a list of random seed parameters
    seedParams = generate_seed_parameters();
      
    // Get the optimised parameters from nelder mead algorithm
    tempParams = simplex.neldermead(&Handler::fitEpidemics, *this,  seedParams);
    // Store the SSE value for this
    tempSSE = fitEpidemics(tempParams);
	
    // If this SSE value is better than the previous, store it and the
    // corresponding parameters
    if(tempSSE < SSE){
      SSE = tempSSE;
      parameters=tempParams;
    }
  }

  // Get the combined values from these parameters, as well as a vector of 
  // each sub-epidemic
  results = ode_solve(parameters);
  tempAll = ode_solve_separate(parameters);
  allResults.clear();
  allResults.push_back(base_model(data));
  for(unsigned int x = 0;x<tempAll.size();++x){
    allResults.push_back(tempAll[x]);
  }
  cout << endl;
  return(SSE);

}

double Handler::fitEpidemics(vector<double> params){
  current_model.clear();
  for(unsigned int i = 0;i<epidemics.size();++i){
    tempParams.clear();
    for(unsigned int z = 0;z<epidemics[i].return_parameters().size();++z){
      tempParams.push_back(epidemics[i].return_parameters()[z]);
    }
    temp_model = epidemics[i].ode_solve(tempParams);
    current_model += combine_vectors(current_model, temp_model);
  }
  return(0.1);
}
  


vector<vector<double> > Handler::ode_solve(vector<double> params){
  vector<vector<double> > overallResults;
  for(unsigned int i = 0;i<epidemics.size();++i){
    tempParams.clear();
    temp_model.clear();
    for(unsigned int z = 0;z<epidemics[i].return_parameters().size();++z){
      tempParams.push_back(epidemics[i].return_parameters()[z]);
    }
    temp_model = epidemics[i].ode_solve(tempParams);
    overallResults += combine_vectors(overallResults, temp_model);
  }
  return(overallResults);
}


vector<vector<vector<double> > > Handler::ode_solve_separate(vector<double> params){
  vector<vector<vector<double> > > separateResults;
  for(unsigned int i = 0;i<epidemics.size();++i){
    tempParams.clear();
    temp_model.clear();
    for(unsigned int z = 0;z<epidemics[i].return_parameters().size();++z){
      tempParams.push_back(epidemics[i].return_parameters()[z]);
    }
    temp_model = epidemics[i].ode_solve(tempParams);
    separateResults.push_back(temp_model);
  }
  return(separateResults);
}

  

/* Takes the residuals as a 2xN vector. Returns true if the latest X residuals are a certain number
   of standard deviations away from the previous total-X residuals */
EpiType Handler::check_epidemic(vector<vector<double> > residuals){
  EpiType detected;
  double mean, sd;
  int smallIndex, largeIndex;
  vector<vector<double> > previousResiduals;
  previousResiduals.clear();
 
  // Take all but the latest 3 residuals
  for(unsigned int j = 0; j < residuals.size() - 3;++j){
    previousResiduals.push_back(residuals[j]);
  }

  smallIndex = previousResiduals.size();
  largeIndex = residuals.size();

  // Calculate the SD and mean of these residuals
  mean = calculate_mean(previousResiduals, 1);
  sd = calculate_sd(previousResiduals, 1);

  // For the latest 3 residuals, check if these satisfy the Z test (ie. within a certain
  // number of standard deviations from the mean of the residuals)
  for(int i = smallIndex; i < largeIndex;++i){
    //cout << "Residual " << i << ": " << residuals[i][1]<< endl;
    if(residuals[i][1] <= (mean + sd*3)){
      return none;
    }
    else if(residuals[i][1] <= (mean + sd*6)){
      detected = sir;
    }
    else{
      detected = sir;
    }
  }
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
  cout << "File imported! Number of rows: " << current_data.size() << endl << endl;
  return double(rows);
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



/* ================================= MATHEMATICAL FUNCTIONS =====================================*/




vector<double> Handler::generate_seed_parameters(){
  vector<double> params;
  for(unsigned int i = 0;i<epidemics.size();++i){
    params = concatenate_vectors(params, rand_params(epidemics[i].return_type()));
  }
  return(params);
}


// Generates a vector of random parameters
// THIS IS WHERE WE GENERATE SEED VALUES
vector<double> Handler::rand_params(EpiType _type){
  vector<double> params;
  double beta,gamma,s0,t0,alpha;
  
  setprecision(9);
  beta = (rand()%100+1)/10000.0;
  alpha = (rand()%100+1)/10000.0;
  gamma = (rand()%100+beta)/1000.0;
  s0 = rand()%1000+100;
  t0 = rand()%80+1;
  params.push_back(log(beta));
  switch(_type){
  case sir:
    params.push_back(log(gamma));
    params.push_back(log(s0));
    params.push_back(log(t0));
    return(params);
  case seir:
    params.push_back(log(alpha));
    params.push_back(log(gamma));
    params.push_back(log(s0));
    params.push_back(log(t0));
    return(params);
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

  xlab = "RSquare: " + to_string(_RSquare);  // Xlabel is the RSquare value
  name = name + _index + ".jpeg";  // Graph name

  gp << "set terminal jpeg size 1000,800 enhanced font \"Helvetica, 10\"\n";    // Here edits type, size and formatting of graph
  gp << "set xlabel '" << xlab << "'\n";
  gp << "set output '" << name << "'\n";  // Set output to the specified file name
  gp << "set termoption dash\n"; // Allow dashes

  // Firstly, plot the actual data (3rd column)
  gp << "plot '-' using 1: " << column << " with linespoints lt 19 title 'Data'";
  
  // Secondly, plot the overall model values
  gp << ", '-' using 1: " << column << " with lines lw 1 title 'Total'";
  
  // For each sub epidemic, add a dashed line to the graph and a corresponding label to the key
  label = "Baseline Level";
  gp << ", '-' using 1: " << column << " with lines lt 0 lc " << 0 << " title '" << label << "'";
 

  for(unsigned int i = 0;i<epidemics.size();++i){
    // Create the label depending on the type of epidemic
    vector<double> par1 = epidemics[i].return_parameters();
    label = "Sub-Epidemic " + to_string(i+1) + " (" + "beta: " + to_string(par1[0]);


    switch(epidemics[i].return_type()){

    case sir:
      label += "; gamma: " + to_string(par1[0]) +  "; S0: " + to_string(par1[2]) 
	+ "; T0: " + to_string(par1[3]) + ")";
      gp << ", '-' using 1: " << epidemics[i].infectedIndex << " with lines lt 0 lc " << i << " title '" << label << "'";
      break;

    default:
      label += "; gamma: " + to_string(par1[0]) +  "; S0: " + to_string(par1[2]) 
	+ "; T0: " + to_string(par1[3]) + ")";
      gp << ", '-' using 1: " << epidemics[i].infectedIndex << " with lines lt 0 lc " << i << " title '" << label << "'";
      break;
    }
  }


  gp << "\n";
  
  // Send the data to Gnuplot
  gp.send1d(data);
  gp.send1d(totalResults);
  for(unsigned int j = 0;j<finalResults.size();++j){
    gp.send1d(finalResults[j]);
  }
}


/* Older function to plot graph of only data and the combined results */
void Handler::plotGraph(vector<vector<double> > finalResults, vector<vector<double> > data, int index){
  Gnuplot gp;
  string name = "graphs/output";
  string _index = to_string(index);
  name = name + _index + ".jpeg";
 
  gp << "set terminal jpeg size 1000,800 enhanced font \"Helvetica, 20\"\n";
  gp << "set output '" << name << "'\n";
  gp << "plot '-' using 1:2 with lines title 'I', '-' using 1:3 with linespoints title 'Data'\n";
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
  SIR mySIR(double(current_data.size()), current_data, sir);
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
