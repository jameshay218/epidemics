#ifndef DATAHANDLER_HPP
#define DATAHANDLER_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

using namespace std;

class Handler{
private:
  vector<vector<double> > current_data;

public:
  Handler();
  ~Handler();
  int import_data(const char* file);
  void print_vector(vector< vector<double> > my_data);
  vector<vector<double> > array_to_vector(double** data);
  void plot_graph(vector<vector<double> > results);
  vector<vector<double> > return_data(){
    return(current_data);
  }
};


#endif
