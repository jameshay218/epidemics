#include <iostream>
#include <string>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <boost/tuple/tuple.hpp>
#include <time.h>
#include <algorithm>
#include <limits>
#include "gnuplot-iostream/gnuplot-iostream.h"
#include <iterator>
#include "simplex.h"
#include "data.h"

using namespace std;

#include "sir.hpp"
#include "datahandler.hpp"

const int STEP_SIZE = 0.1;


template<class Con>
void printcon(const Con& c){
  std::cout.precision(12);
  copy( c.begin(), 
	c.end(), 
	ostream_iterator<typename Con::value_type>(cout, "  ") );
  cout<<endl;
}

extern double sir_sse(vector<double>);

int main() {
  Gnuplot gp;
  Handler dataHandler;
  SIR* mySIR;
  vector<vector<double> > results;
  int dataRows;
  char file[80];
    
  cout << "Specify a file: " << endl;
  cin >> file;
  cout << "Specified " << file << endl;
  dataRows = dataHandler.import_data(file);
  
  results.resize(dataRows);
  for(int i=0;i<dataRows;i++){
      results[i].resize(4);
  }
  double parameters[7] = {0.001, 0.1, 1.0, 500.0, 1.0, 0.0, double(dataRows)};
  
  mySIR = new SIR(parameters, dataHandler.return_data());
  
  
  cout << "------------- Now trying to do optimisation -----------------" << endl;
  clock_t t1, t2;
  t1=clock();
  using Measurement::GetData;
  vector<double> pars;
  pars.push_back(0.001);
  pars.push_back(0.1);
  pars.push_back(1.0);
  pars.push_back(500.0);
  pars.push_back(1.0);
  pars.push_back(0.0);
  pars.push_back(double(dataRows));
  
  GetData(dataHandler.return_data(), pars);
  cout << "Got... data?" << endl;

  
  double par1[3] = {log(0.0013), log(0.013), log(701)};
  //vector<double> init(par1, par1+3);
  double SSE = 999999999.9;
  
  vector<double> resultPar = mySIR->rand_params();

  for(int index=0;index<=50;index++){
    vector<double> par2;
    par2 = mySIR->rand_params();
    //cout << "Starting parameters: ";
    /*for(unsigned int i = 0; i<par2.size();i++){
      cout << (double)exp(par2[i]) << ' ';
    }
    cout << endl;
    */
    using BT::Simplex;
    
    vector<double> result =  Simplex(Measurement::sse_sir, par2);

    for(unsigned int j=0;j<result.size();++j){
      par1[j] = result[j];
    }

    if(mySIR->sse_sir(par1) < SSE){
      SSE = mySIR->sse_sir(par1);
      resultPar=result;
    }
  
  }
    
  cout << "Final SSE was: " << SSE << endl;
  for(unsigned int j=0;j<resultPar.size();++j){
    par1[j] = resultPar[j];
    resultPar[j] = exp(resultPar[j]);
  }
  cout << "Final parameters: ";
  printcon(resultPar) ; cout << endl;
  
  results = mySIR->sse_sir2(par1);  
  t2=clock();

  gp << "set output 'my_graph_1.png'\n";
  gp << "plot '-' using 1:3 with lines title 'I', '-' using 1:3 with lines title 'Data'\n";
  gp.send1d(results);
  gp.send1d(dataHandler.return_data());
  

  cout << "Time taken: " << (t2-t1) << endl;
  
  
  delete mySIR;
   
  return(0); 
}

