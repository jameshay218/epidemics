#include <iostream>
#include <string>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <boost/tuple/tuple.hpp>
#include <time.h>
#include "gnuplot-iostream/gnuplot-iostream.h"

using namespace std;

#include "sir.hpp"
#include "datahandler.hpp"

const int STEP_SIZE = 0.1;

int main() {
  Gnuplot gp;
  Handler dataHandler;
  SIR* mySIR;
  vector<vector<double> > results;
  int dataRows;
  char file[80];
  vector< vector<int> > my_data;
  
  cout << "Specify a file: " << endl;
  cin >> file;
  cout << "Specified " << file << endl;
  dataRows = dataHandler.import_data(file);
  results.resize(dataRows);
  for(int i=0;i<=dataRows;++i){
    results[i].resize(4);
  }
  double parameters[7] = {0.001, 0.1, 1.0, 500.0, 1.0, 0.0, double(dataRows)};
  mySIR = new SIR(parameters);
  mySIR->Solve_Eq(results);
  cout << mySIR->add_arrays(results,dataHandler.return_data()) << endl;
  gp << "set output 'my_graph_1.png'\n";
  gp << "plot '-' using 1:3 with lines title 'I', '-' using 1:2 with lines title 'S', '-' using 1:4 with lines title 'R', '-' using 1:2 with lines title 'Data'\n";
  gp.send1d(results);
  gp.send1d(results);
  gp.send1d(results);
  gp.send1d(dataHandler.return_data());

  cout << "OK now solving" << endl;
  double par[3]={log(0.001), log(0.1), log(500)};
  clock_t t1, t2;
  t1=clock();
  for(int k =0;k<=10000;k++){
  mySIR->sse_sir(par, dataHandler.return_data());
  }
  t2=clock();
  cout << (t2-t1)/CLOCKS_PER_SEC << endl;
  delete mySIR;
  
  
  return(0); 
}

