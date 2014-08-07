using namespace std;

#include "datahandler.hpp"
#include "sir.hpp"
#include "simplex.hpp"
#include <string>
Handler::Handler(){
}

Handler::~Handler(){
  cout << "Delete Data Handler" << endl;
}

void Handler::likelihood_test(vector<double> &params){
  vector<double> parameters, tempParams;
  vector<vector<int> > data, model;
  SIR mySIR(double(current_data.size()), current_data);
  Simplex simplex;
  parameters.push_back(log(0.001));
  parameters.push_back(log(0.1));
  parameters.push_back(log(500));
  parameters.push_back(log(10));
  
  data = mySIR.combined_model(parameters);
  print_vector(data);
  //long temp = 222;
  //mySIR.factorial(temp);
  /*mySIR.update_data(data);

  parameters.clear();
  parameters.push_back(log(0.0005));
  parameters.push_back(log(0.11));


  cout << "Test: " << mySIR.mle_sir(parameters) << endl;

  tempParams = simplex.neldermead(&SIR::mle_sir, mySIR,  parameters);
  cout << "Final log likelihood: " << mySIR.mle_sir(tempParams) << endl;
  tempParams.push_back(log(500));
  tempParams.push_back(log(10));
  parameters.clear();
  parameters.push_back(log(0.001));
  parameters.push_back(log(0.1));
  cout << "Likelihood of actual parameters: " << mySIR.mle_sir(parameters) << endl;
  cout << "SSE: " << mySIR.calculate_SSE(model,data) << endl;
  parameters.push_back(log(500));
  parameters.push_back(log(10));
  model = mySIR.sse_sir_combined(parameters);
  for(int i = 0;i<tempParams.size();++i){
    tempParams[i] = exp(tempParams[i]);
  }
  printcon(tempParams);
  plotGraph(model, data, 1);
  */
}



/* ================================ FINAL LEAST SQUARES ============================= */


void Handler::realtime_fit_multi(vector<double> &params, int version, double targetRSq){
  clock_t t1, t2; // To record the total run time
  vector<double> finalParams, finalParamsTemp;
  vector<int> detectionTimes;
  vector<vector<double> > temp, tempModel, currentModel, combinedResults, residuals, combinedResultsTemp;
  vector<vector<vector<double> > > componentResults, componentResultsTemp;
  int epidemicCount = 1;
  double RSquare, RSquareK1, SSE;
  // Start of model fitting process
  t1=clock();

  // Start off with a baseline model of the mean of the first 4 points
  for(unsigned int j = 0; j <= 4; ++j){
    temp.push_back(current_data[j]);
  }
  currentModel = base_model(temp);
  
  // For each time point in the current data set, carry out the fitting procedure
  for(unsigned int i = 50;i<current_data.size();++i){
    cout << endl << "-----------------" << endl;
    cout << "Iteration number " << i << endl;
    
    // Clear temp vector
    temp.clear();
    tempModel.clear();

    // Store all data up to the current index and give this to SIR
    for(unsigned int j = 0; j < i; ++j){
      temp.push_back(current_data[j]);
      tempModel.push_back(currentModel[j]);
      
    } 

    
    // Create an optimised model fit for the given data (temp) and calculate RSquare from the SSE
    SSE = optimiseEpidemics(epidemicCount, finalParams, temp, combinedResults, componentResults);
    RSquare = 1 - SSE/SStot(temp, 2);
    

    // If fit is sufficiently good with k epidemics, try fitting with k-1 epidemics. If this fit
    // is sufficient, use k-1 epidemics.
    if(RSquare > targetRSq && epidemicCount > 0){
      cout << "Considering removal" << endl;
      SSE = optimiseEpidemics((epidemicCount-1),finalParamsTemp, temp, combinedResultsTemp, componentResultsTemp);
      RSquareK1 = 1 - SSE/SStot(temp,2);
      if(RSquareK1 > targetRSq){
	cout << "Removing an epidemic" << endl;
	RSquare = RSquareK1;
	finalParams = finalParamsTemp;
	combinedResults = combinedResultsTemp;
	componentResults = componentResultsTemp;
	epidemicCount--;
      }
    }

    // Check if a new epidemic has started. If so, and RSquare from current K epidemics 
    // insufficient, try fitting another. If this is a significantly better fit, permanently
    // add an extra epidemic.
    residuals = get_residuals(temp, tempModel);
    if(check_epidemic(residuals) && RSquare < targetRSq){
      cout << "Epidemic detected!" << endl;
      detectionTimes.push_back(i);
      SSE = optimiseEpidemics((epidemicCount+1),finalParamsTemp, temp, combinedResultsTemp, componentResultsTemp);
      RSquareK1 = 1 - SSE/SStot(temp,2);
      if(RSquareK1 > targetRSq){
	RSquare = RSquareK1;
	finalParams = finalParamsTemp;
	combinedResults = combinedResultsTemp;
	componentResults = componentResultsTemp;
	epidemicCount++;
      }
    }
   

    // Transform parameters back to normal space
    for(unsigned int j=0;j<finalParams.size();++j){
      finalParams[j] = exp(finalParams[j]);
    }
    // Plot graph
    plotGraphMulti(componentResults, combinedResults, temp, i, finalParams, RSquare, detectionTimes);    
    
    // Print out results for this iteration
    cout << "Final SSE was: " << SSE << endl;
    cout << "Final RSquare: " << RSquare << endl;
    cout << "Final parameters: ";
    printcon(finalParams) ;
    cout << endl << "-----------------" << endl;
    currentModel = combinedResults;
  }
  t2=clock();
  cout << "Time taken: " << (t2-t1)/CLOCKS_PER_SEC << endl << endl;
}



double Handler::optimiseEpidemics(int epiCount, vector<double> &parameters, vector<vector<double> > data, vector<vector<double> > &results, vector<vector<vector<double> > > &allResults){
  SIR mySIR(double(current_data.size()), current_data);
  double SSE = 99999999999.9;
  double tempSSE;
  vector<double> tempParams, seedParams;
  vector<vector<double> > currentModel;
  vector<vector<vector<double> > > tempAll;
  Simplex simplex;

  mySIR.update_data(data);

  // If no epidemics have yet been detected, use the mean of the current data as the 
  // current model and return the corresponding SSE.
  if(epiCount == 0){
    currentModel = base_model(data);
    SSE = mySIR.calculate_SSE(currentModel, data);
    parameters.clear();
    results = currentModel;
    allResults.clear();
    allResults.push_back(currentModel);
    return(SSE);
  }

  // If there are epidemics to be fitted, perform 40 random fits and keep best fitting model
  for(int index=0;index<20;index++){
    cout << ".";
    //cout << "Optimise " << index << endl;
    // Clear the temporary seed parameters and seed rand
    seedParams.clear();      
    srand(clock());

    // Create a list of random seed parameters
    for(int j =1; j<=epiCount;++j){
      seedParams = concatenate_vectors(seedParams, mySIR.rand_params4());
    }
    
    // Get the optimised parameters from nelder mead algorithm
    tempParams = simplex.neldermead(&SIR::sse_sir_multi, mySIR,  seedParams);
    // Store the SSE value for this
    tempSSE = mySIR.sse_sir_multi(tempParams);
	
    // If this SSE value is better than the previous, store it and the
    // corresponding parameters
    if(tempSSE < SSE){
      SSE = tempSSE;
      parameters=tempParams;
    }
  }

  // Get the combined values from these parameters, as well as a vector of 
  // each sub-epidemic
  results = mySIR.sse_sir_combined(parameters); 
  currentModel = base_model(data);
  tempAll = mySIR.sse_sir_components(parameters);
  allResults.clear();
  allResults.push_back(currentModel);
  for(unsigned int x = 0;x<tempAll.size();++x){
    allResults.push_back(tempAll[x]);
  }
  cout << endl;
  return(SSE);

}





/* Takes the residuals as a 2xN vector. Returns true if the latest X residuals are a certain number
   of standard deviations away from the previous total-X residuals */
bool Handler::check_epidemic(vector<vector<double> > residuals){
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
  /*
  cout << "Residuals mean: " << mean << endl;
  cout << "Residuals SD: " << sd << endl;
  */

  // For the latest 3 residuals, check if these satisfy the Z test (ie. within a certain
  // number of standard deviations from the mean of the residuals)
  for(int i = smallIndex; i < largeIndex;++i){
    //cout << "Residual " << i << ": " << residuals[i][1]<< endl;
    if(residuals[i][1] <= (mean + sd*3)){
      return false;
    }
    
  }
  return true;
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
vector<vector<double> > Handler::get_residuals(vector<vector<double> > data1, vector<vector<double> > data2){
  unsigned int i = 0;
  unsigned int j = 0;
  vector<vector<double> > results;
  vector<double> row;
  while(i < data1.size() && j < data2.size()){
    //If indices give same time point, use difference
    if(data1[i][0] == data2[j][0]){
      row.clear();
      row.push_back(data1[i][0]);
      row.push_back((data1[i][2] - data2[j][2]));
      results.push_back(row);
      i++;
      j++;
    }
    //If first dataset is still before second, full difference
    else if(data1[i][0] < data2[j][0]){
      row.clear();
      row.push_back(data1[i][0]);
      row.push_back(data1[i][2]);
      results.push_back(row);
      i++;
    }
    //...
    else{
      row.clear();
      row.push_back(data2[j][0]);
      row.push_back(data2[j][2]);
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
void Handler::plotGraphMulti(vector<vector<vector<double> > > finalResults, vector<vector<double> > totalResults, vector<vector<double> > data, int index, vector<double> parameters, double _RSquare, vector<int> _detected){
  Gnuplot gp;   // Need instance of the Gnuplot class to pipe commands to gnuplot
  string name = "graphs1/output"; // The save location and general name of the graph to be saved
  string _index = to_string(index);
  string label, xlab;

  xlab = "RSquare: " + to_string(_RSquare);  // Xlabel is the RSquare value
  name = name + _index + ".jpeg";  // Graph name

  gp << "set terminal jpeg size 1000,800 enhanced font \"Helvetica, 10\"\n";    // Here edits type, size and formatting of graph
  gp << "set xlabel '" << xlab << "'\n";
  gp << "set output '" << name << "'\n";  // Set output to the specified file name
  gp << "set termoption dash\n"; // Allow dashes

  // OFF
  // Add lines for detection times 
  /*for(unsigned int i =0;i<_detected.size();++i){
    gp << "set arrow from " << _detected[i] << ",0 to " << _detected[i] << "," << 500 << " nohead lc rgb 'red'\n";
    } */

  // Firstly, plot the actual data (3rd column)
  gp << "plot '-' using 1:3 with linespoints lt 19 title 'Data'";
  
  // Secondly, plot the overall model values
  gp << ", '-' using 1:3 with lines lw 1 title 'Total'";
  
  // For each sub epidemic, add a dashed line to the graph and a corresponding label to the key
  label = "Baseline Level";
  gp << ", '-' using 1:3 with lines lt 0 lc " << 0 << " title '" << label << "'";
  for(unsigned int i = 1;i<finalResults.size();++i){
    // Create the label
    label = "Sub-Epidemic " + to_string(i) + " (" + "beta: " + to_string(parameters[4*(i-1)]) 
      + "; gamma: " + to_string(parameters[1+4*(i-1)]) +  "; S0: " + to_string(parameters[2+4*(i-1)]) 
      + "; T0: " + to_string(parameters[3+4*(i-1)]) + ")";

    gp << ", '-' using 1:3 with lines lt 0 lc " << i << " title '" << label << "'";
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
  gp << "plot '-' using 1:3 with lines title 'I', '-' using 1:3 with linespoints title 'Data'\n";
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



/* ======================================== OLD ========================================= */


// OLD FUNCTION
/* Old test function to test the fitting of a single epidemic model */
void Handler::realtime_fit(vector<vector<double> > &results, vector<double> &params, int version){
  clock_t t1, t2;
  SIR mySIR(double(current_data.size()), current_data);
  double SSE;
  vector<double> resultPar;
  vector<vector<double> > temp;
  double tempSSE, RSquare;
  Simplex simplex;
  
  t1=clock();
 
  /* Added artificial t0 delay of 15 */
  for(unsigned int i = 0; i<current_data.size();++i){
    current_data[i][0] += 15;
    }
  

  for(unsigned int i = 60;i<current_data.size();++i){
    SSE = 999999999999.9;
    cout << endl << "-----------------" << endl;
    cout << "Iteration number " << i << endl;

    // Store all data up to the current index and give this to SIR
    temp.clear();
    for(unsigned int j = 0; j < i; ++j){
      temp.push_back(current_data[j]);
    } 
    mySIR.update_data(temp);
    
    // Optimise fit 50 times
    for(int index=0;index<=100;index++){
      //cout << "Optimise " << index << endl;
      srand(clock());
      vector<double> par2;
      // Optimise for 3 or 4 parameters
      switch(version){
      case 1:
	par2= mySIR.rand_params3();
	resultPar = simplex.neldermead(&SIR::sse_sir, mySIR,  par2);
	tempSSE = mySIR.sse_sir(resultPar);
	break;
      case 2:
	par2 = mySIR.rand_params4();
	resultPar = simplex.neldermead(&SIR::sse_sir_t0, mySIR,  par2);
	tempSSE = mySIR.sse_sir_t0(resultPar);
	break;
     
      }
   
      if(tempSSE < SSE){
	SSE = tempSSE;
	params=resultPar;
      }
    }
    
    cout << "Final SSE was: " << SSE << endl;
    results = mySIR.sse_sir_single(params);  
    RSquare = 1 - SSE/SStot(temp, 2);
    cout << "Final RSquare: " << RSquare << endl;
    for(unsigned int j=0;j<params.size();++j){
      params[j] = exp(params[j]);
    }
    
    cout << "Final parameters: ";
    printcon(params) ;
    cout << "-----------------" << endl;
    plotGraph(results, temp, i);
    //print_vector(temp);
  }
  t2=clock();
  cout << "Time taken: " << (t2-t1)/CLOCKS_PER_SEC << endl << endl;
}



// OLD FUNCTION
/* Older function to test whether the SIR vector addition function is working correctly */
void Handler::testAddition(vector<vector<double> > data1, vector<vector<double> > data2, double offset){
  
  SIR mySIR(double(current_data.size()), current_data);
  for(unsigned int i = 0; i<data2.size();++i){
    data2[i][0] += (int)offset;
  }
  cout << mySIR.calculate_SSE(data1,data2) << endl;
  plotGraph(data1,data2,100);

}



// OLD FUNCTION
/* Old function to test the ability to detect the start of new epidemics */
void Handler::test_detect(vector<vector<double> > &results, vector<double> &params){
  SIR mySIR(double(current_data.size()), current_data);
  vector<double> par, par1, mean, sd;
  vector<vector<double> > myData, currentData, currentModel, temp, residuals, tempResiduals;
  vector<vector<vector<double> > > temp1;
  int epidemicCount = 0;
  vector<int> detectionTimes;
  par.push_back(log(0.001));
  par.push_back(log(0.1));
  par.push_back(log(500));
  par.push_back(log(10));
  par.push_back(log(0.001));
  par.push_back(log(0.05));
  par.push_back(log(600));
  par.push_back(log(50));

  myData = mySIR.sse_sir_combined(par);

  par1.push_back(log(0.001));
  par1.push_back(log(0.1));
  par1.push_back(log(500));
  par1.push_back(log(10));

  currentData = mySIR.sse_sir_combined(par1);

  temp1.push_back(currentData);
  plotGraphMulti(temp1, myData, myData, 1, par, 1.0, detectionTimes);

  //residuals = mySIR.get_residuals(currentData,myData);
  

  

  for(unsigned int i = 4;i<100;i++){
    cout << "Iteration " << i << endl;
    temp.clear();
    tempResiduals.clear();
    for(unsigned int j = 0; j < i; ++j){
      temp.push_back(myData[j]);
    } 
    if(epidemicCount == 0){
      currentModel = currentData;
    }
    else{
      currentModel = myData;
    }
    for(unsigned int j = 0; j < i; ++j){
      tempResiduals.push_back(currentModel[j]);
    } 
    residuals.clear();
    residuals = get_residuals(temp, tempResiduals);
    //print_vector(residuals);
    if(check_epidemic(residuals) == true){
      cout << "Epidemic detected at time " << i << endl;
      epidemicCount++;
    }
    cout << endl;
  }

}
