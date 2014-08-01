#include <iostream>
#include <string>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <boost/tuple/tuple.hpp>
#include <time.h>
#include <algorithm>
#include <limits>

#include <iterator>
//#include "simplex.h"

using namespace std;

#include "datahandler.hpp"
#include "sir.hpp"
const int STEP_SIZE = 0.1;



int main() {
  
  Handler dataHandler;
  vector<vector<double> > finalResults;
  vector<double> finalParams;
  char file[80];
  int option;      

  cout << "====================================================================" << endl;
  cout << "                            SIR Fitting                             " << endl;
  cout << "====================================================================" << endl << endl;
  /*
  dataHandler.import_data("test.csv",finalResults);
  double parameters[7] = {0.001, 0.1, 1.0, 500.0, 1.0, 0.0, dataHandler.return_data().size()};  
  SIR mySIR(parameters, dataHandler.return_data());
  dataHandler.print_vector(dataHandler.return_data());
  finalResults = mySIR.combine_vectors(dataHandler.return_data(), dataHandler.return_data());
  dataHandler.print_vector(finalResults);
  return 0;*/
  while(true){
    /* User specifies a file which is store in an instance of the SIR class */
    cout << "--------------------------------------------------------------------" << endl;
    cout << "Please specify the type of epidemic that you would like to fit:     " << endl;
    cout << "          1. Single Epidemic                                        " << endl;
    cout << "          2. Multiple Epidemic                                        " << endl;
    cout << "--------------------------------------------------------------------" << endl;
    cin >> option;
    switch(option){
    case 1:
      cout << "Please specify a data file (.csv): " << endl;
      cout << "--------------------------------------------------------------------" << endl << endl;
      cin >> file;
      cout << "--------------------------------------------------------------------" << endl;
      cout << "Specified: " << file << endl;
      cout << "--------------------------------------------------------------------" << endl;
      
      cout << "--------------------------------------------------------------------" << endl;
      cout << "********************* FITTING MODEL PARAMETERS *********************" << endl;
      cout << "--------------------------------------------------------------------" << endl;
      dataHandler.import_data(file, finalResults);
      //dataHandler.plotGraph(dataHandler.return_data(), finalResults, 250);
      //dataHandler.testAddition(dataHandler.return_data(), dataHandler.return_data(), 10.31521361);
      dataHandler.realtime_fit(finalResults, finalParams, 2);
      break;
    case 2:
      cout << "Please specify a data file (.csv): " << endl;
      cout << "--------------------------------------------------------------------" << endl << endl;
      cin >> file;
      cout << "--------------------------------------------------------------------" << endl;
      cout << "Specified: " << file << endl;
      cout << "--------------------------------------------------------------------" << endl;
      
      cout << "--------------------------------------------------------------------" << endl;
      cout << "********************* FITTING MODEL PARAMETERS *********************" << endl;
      cout << "--------------------------------------------------------------------" << endl;
      dataHandler.import_data(file, finalResults);
      dataHandler.realtime_fit2(finalResults, finalParams, 2);
      break;
    }

  
  }
  return(0); 
}

