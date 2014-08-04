#ifndef DATAHANDLER_HPP
#define DATAHANDLER_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iterator>
#include "gnuplot-iostream/gnuplot-iostream.h"

using namespace std;


template<class Con>
void printcon(const Con& c){
  std::cout.precision(12);
  copy( c.begin(), 
	c.end(), 
	ostream_iterator<typename Con::value_type>(cout, "  ") );
  cout<<endl;
}


class Handler{
private:
  vector<vector<double> > current_data;

public:
  Handler();
  ~Handler();
  double import_data(const char* file, vector<vector<double> > _results);
  void print_vector(vector< vector<double> > my_data);
  void plotGraph(vector<vector<double> > finalResults, vector<vector<double> > data, int index);
  void plotGraph2(vector<vector<vector<double> > > finalResults, vector<vector<double> > totalResults, vector<vector<double> > data, int index, vector<double> parameters, double _RSquare);

  void realtime_fit(vector<vector<double> > &results, vector<double> &params, int version);
  void realtime_fit2(vector<vector<double> > &results, vector<double> &params, int version);
  vector<vector<double> > return_data(){
    return(current_data);
  }
  vector<double> concatenate_vectors(vector<double> a, vector<double> b);
  double SStot(vector<vector<double> > data, int column);
  double calculate_mean(vector<vector<double> > data, int column);
  void testAddition(vector<vector<double> > data1, vector<vector<double> > data2, double offset);
};


#endif
