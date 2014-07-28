#include <iostream>
#include <string>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <boost/tuple/tuple.hpp>

#include "gnuplot-iostream/gnuplot-iostream.h"

using namespace std;

#include "sir.hpp"


int main() {
  Gnuplot gp;
  int option;
  char file[80];
  vector< vector<int> > my_data;

  while(true){
    
    cout << "Please select a function:" << endl << endl;
    cout << "         1. Import data." << endl;
    cout << "         2. Fit model." << endl;
    cout << "         3. Plot data." << endl;
    cout << "         0. Exit." << endl;
    cin >> option;
    
    
    switch(option){
    case 1:
      cout << "Specify a file:" << endl;
      cin >> file;
      cout << "Specified " << file << endl;
      import_data(file, my_data);
      gp << "set output 'my_graph_1.png'\n";
      gp << "plot '-' with lines\n";
      gp.send1d(my_data);
      break;
    case 2:
      cout << "Fit model." << endl;
      
      break;
    case 3:
      
      break;
    case 0:
      cout << "Exiting program" << endl;
      exit(1);
    }
    

  }
  return 0;
}
