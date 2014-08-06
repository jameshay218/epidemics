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
 
  void realtime_fit_multi(vector<double> &params, int version);
  bool check_epidemic(vector<vector<double> > residuals);
  double optimiseEpidemics(int epiCount, vector<double> &parameters, vector<vector<double> > data, vector<vector<double> > &results, vector<vector<vector<double> > > &allResults);
  // Data handling fucntions
  vector<vector<double> > return_data(){
    return(current_data);
  }
  vector<double> concatenate_vectors(vector<double> a, vector<double> b);
  double import_data(const char* file);
  void print_vector(vector< vector<double> > my_data);

  // Maths functions
  vector<vector<double> > base_model(vector<vector<double> > data);
  double SStot(vector<vector<double> > data, int column);
  double calculate_mean(vector<vector<double> > data, int column);
  double calculate_sd(vector<vector<double> > data, int column);
  vector<vector<double> > get_residuals(vector<vector<double> > data1, vector<vector<double> > data2);
  
  // Graph plotting functions
  void plotGraph(vector<vector<double> > finalResults, vector<vector<double> > data, int index);
  void plotGraphMulti(vector<vector<vector<double> > > finalResults, vector<vector<double> > totalResults, vector<vector<double> > data, int index, vector<double> parameters, double _RSquare, vector<int> _detection);


  // Old functions
  void testAddition(vector<vector<double> > data1, vector<vector<double> > data2, double offset);
  void test_detect(vector<vector<double> > &results, vector<double> &params);
  void realtime_fit(vector<vector<double> > &results, vector<double> &params, int version);

};


#endif
