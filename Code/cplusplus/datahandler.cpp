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


double Handler::import_data(const char* file, vector<vector<double> > _results){
  ifstream in_stream;
  string temp;
  current_data.clear();
  int rows = 0;
  cout << "Importing file... " << file << endl;
  current_data.clear();
  in_stream.open(file);
  if (in_stream.fail()) {
    cout << "Error opening file" << endl;
    return -1;
  }
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
  
  _results.resize(rows);
  for(int i=0;i<rows;i++){
      _results[i].resize(4);
  }

  return double(rows);
}

void Handler::print_vector(vector< vector<double> > my_data){
  for(vector< vector<double> >::const_iterator i = my_data.begin(); i !=my_data.end(); ++i){
    for(vector<double>::const_iterator j = i->begin(); j!=i->end();++j){
      cout << *j << ' ';
    }
    cout << endl;  
  }
}
  
void Handler::realtime_fit(vector<vector<double> > &results, vector<double> &params, int version){
  clock_t t1, t2;
  double parameters[7] = {0.001, 0.1, 1.0, 500.0, 1.0, 0.0, double(current_data.size())};  
  SIR mySIR(parameters, current_data);
  double SSE;
  vector<double> resultPar;
  vector<vector<double> > temp;
  double tempSSE;
  Simplex simplex;
  
  t1=clock();
  
  for(int i = 0; i<current_data.size();++i){
    current_data[i][0] += 15;
  }

  for(unsigned int i = 90;i<current_data.size();++i){
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
    for(int index=0;index<=50;index++){
      cout << "Optimise " << index << endl;
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
      
      cout << "Current SSE: " << tempSSE << endl;
      cout << "Current params: ";
      printcon(resultPar); cout << endl;
      if(tempSSE < SSE){
	SSE = tempSSE;
	params=resultPar;
      }
    }
    
    cout << "Final SSE was: " << SSE << endl;
    results = mySIR.sse_sir2(params);  
    
    for(unsigned int j=0;j<params.size();++j){
      params[j] = exp(params[j]);
    }
    
    cout << "Final parameters: ";
    printcon(params) ;
    cout << "-----------------" << endl;
    plotGraph(results, temp, i);
  }
  t2=clock();
  cout << "Time taken: " << (t2-t1)/CLOCKS_PER_SEC << endl << endl;
}


  

  
void Handler::plotGraph(vector<vector<double> > finalResults, vector<vector<double> > data, int index){
  Gnuplot gp;
  string name = "graphs/output";
  string _index = to_string(index);
  name = name + _index + "jpeg";
 
  gp << "set terminal jpeg size 1000,800 enhanced font \"Helvetica, 20\"\n";
  gp << "set output '" << name << "'\n";
  gp << "plot '-' using 1:3 with lines title 'I', '-' using 1:3 with lines title 'Data'\n";
  gp.send1d(finalResults);
  gp.send1d(data);
  
}

void Handler::testAddition(vector<vector<double> > data1, vector<vector<double> > data2, double offset){
  double parameters[7] = {0.001, 0.1, 1.0, 500.0, 1.0, 0.0, double(current_data.size())};  
  SIR mySIR(parameters, current_data);
  for(int i = 0; i<data2.size();++i){
    data2[i][0] += (int)offset;
  }
  cout << mySIR.add_arrays(data1,data2) << endl;
  plotGraph(data1,data2,100);

}
