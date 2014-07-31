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
  void realtime_fit(vector<vector<double> > &results, vector<double> &params, int version);

vector<vector<double> > return_data(){
    return(current_data);
  }

  void testAddition(vector<vector<double> > data1, vector<vector<double> > data2, double offset);
};


#endif
