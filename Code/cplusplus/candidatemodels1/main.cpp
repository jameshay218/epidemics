#include <iostream>
#include <string>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <boost/tuple/tuple.hpp>
#include <boost/filesystem.hpp>
#include <time.h>
#include <algorithm>
#include <limits>

#include <iterator>
//#include "simplex.h"

using namespace std;

#include "datahandler.hpp"
#include "epidemic.hpp"
#include "sir.hpp"
#include "seir.hpp"

const int STEP_SIZE = 0.1;


bool check_existence(char place[80]){
  if(boost::filesystem::exists(place) && boost::filesystem::is_directory(place)){
    return true;
  }
  return false;
}
    
EpiType convert_to_epi_type(int epi){
  switch(epi){
  case 0: return none;
  case 1: return sir;
  case 2: return irsir;
  case 3: return seir;
  case 4: return serir;
  case 5: return spike;
  default: return sir;
  }
  return sir;
}


int main() {
  Handler dataHandler;
  EpiType epiType = sir;
  vector<vector<double> > finalResults;
  vector<double> finalParams;
  char file[80], saveLocation[80];
  int option, fitting, testEpiType = 0;
  double targetRSq = 0.95;
  char defaultOptions,createDir,useT0,useI0;      
  bool defaults, mle=false,boolT0=true,boolI0=false,singleEpidemic=false, savePlot=true,saveResults=false;
  ifstream file1;
  boost::filesystem::path dir;
  cout << "====================================================================" << endl;
  cout << "                Welcome to the SIR Fitting Framework                       " << endl;
  cout << "====================================================================" << endl << endl;
  
  while(true){
    /* User specifies a data file */
    cout << "--------------------------------------------------------------------" << endl;
    cout << "What would you like to do?                                          " << endl;
    cout << "          1. Fit single epidemic                                    " << endl;
    cout << "          2. Fit synthedemic                                       " << endl;
    cout << "--------------------------------------------------------------------" << endl;
    cin >> option;
    cout << "Would you like to use default options (y/n)?" << endl;
    cin >> defaultOptions;

    switch(defaultOptions){
    case 'y':
      defaults = true;
      break;
    case 'n':
      cout << "Using custom options" << endl;
      defaults=false;
      break;
    default:
      cout << "Invalid response. Using default options" << endl;
      defaults=true;
      break;
    }
      

 
    cout << "Please specify a data file (.csv): " << endl;
    cout << "--------------------------------------------------------------------" << endl;

    cin >> file;
    cout << "--------------------------------------------------------------------" << endl;
    cout << "Specified: " << file << endl;
    cout << "--------------------------------------------------------------------" << endl;
    while(dataHandler.import_data(file) == -1){
      cout << endl << "Error - file invalid. Please select another file" << endl;
      cin >> file;
	
    }
    if(!defaults){
      cout << "Please specify a save location: " << endl;
      cin >> saveLocation;
      cout << "--------------------------------------------------------------------" << endl;
      cout << "Specified: " << saveLocation << endl;
      cout << "--------------------------------------------------------------------" << endl;
      while(!check_existence(saveLocation)){
	cout << endl << "Error - no such directory exists. Would you like to create it (y/n)?" << endl;
	cin >> createDir;
	switch(createDir){
	case 'y':
	  dir = saveLocation;
	  if(boost::filesystem::create_directory(dir)) cout << "Created directory named \"" << saveLocation << "\"" << endl;
	  else{ 
	    cout << "Error creating directory. Please specify another location:" << endl;
	    cin >> saveLocation;
	  }
	  break;
	case 'n':
	  cout << "Please specify another directory:" << endl;
	  cin >> saveLocation;
	  break;
	default:
	  cout << "Invalid option" << endl;
	  cout << "Please specify another directory:" << endl;
	  cin >> saveLocation;
	  break;
	}
      }
      cout << "Use MLE or SSE based fitting (1/2)?" << endl;
      cin >> fitting;
      switch(fitting){
      case 1:
	mle=true;
	break;
      case 2 :
	mle=false;
	break;
      default:
	mle=false;
	break;
      }
      cout << "Would you like to include I0 in the optimisation (y/n)?" << endl;
      cin >> useI0;
      if(useI0 == 'y') boolI0 = true;
      
      cout << "What is your target R-Squared? (between 0 and 1, default is 0.95)" << endl;
      cin >> targetRSq;
      
      

    }
    
    switch(option){
    case 1:
      singleEpidemic = true;
      boolT0 = false;
      boolI0 = false;
      mle = false;
      if(!defaults){
	cout << "Would you like to include t0 in the optimisation (y/n)?" << endl;
	cin >> useT0;
	if(useT0 == 'y') boolT0 = true;
	else boolT0 = false;
	
	cout << "Which epidemic type would you like to test?" << endl;
	cout << "     0. Unknown (attempts to detect best fitting type)" << endl;
	cout << "     1. SIR" << endl;
	cout << "     2. irSIR" << endl;
	cout << "     3. SEIR" << endl;
	cout << "     4. SERIR" << endl;
	cout << "     5. EXP" << endl;
	
	cin >> testEpiType;

	while(testEpiType < 0 || testEpiType > 5){
	  cout << "Invalid type. Please select another: " << endl;
	  cin >> testEpiType;
	}
	epiType = convert_to_epi_type(testEpiType);
      }
      dataHandler.update_options(mle,boolT0,boolI0,singleEpidemic,savePlot,saveResults);
      dataHandler.realtime_fit_single(targetRSq, epiType);
      cout << "--------------------------------------------------------------------" << endl;
      cout << "********************* FITTING MODEL PARAMETERS *********************" << endl;
      cout << "--------------------------------------------------------------------" << endl;
      break;
    case 2:
      singleEpidemic = false;
      cout << "--------------------------------------------------------------------" << endl;
      cout << "********************* FITTING MODEL PARAMETERS *********************" << endl;
      cout << "--------------------------------------------------------------------" << endl;
      //      boolT0 = false;
      dataHandler.update_options(mle,boolT0,boolI0,singleEpidemic,savePlot,saveResults);
      dataHandler.realtime_fit_multi(0.8);
      //dataHandler.overall_test(0.95);
      break;
    }
  
  }
  return(0); 
}

