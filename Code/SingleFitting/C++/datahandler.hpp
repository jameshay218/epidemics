#ifndef DATAHANDLER_HPP
#define DATAHANDLER_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iterator>
#include "gnuplot-iostream/gnuplot-iostream.h"

#include "epidemic.hpp"

using namespace std;

// Function to print a vector of type Con
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
  vector<double> tempParams;
  vector<vector<double> > current_data, current_model, temp_model, empty_model, temp_data, baseModel;
  vector<Epidemic*> epidemics, tempEpidemics;
  bool useMLE, optimT0, optimI0, singleEpi, plot, save;
  string saveLocation = "graphs/";

public:
  Handler();
  ~Handler();
  
  
  // Optimisation functions
  EpiType check_epidemic(vector<vector<double> > residuals);
  double optimiseEpidemics(vector<double> &parameters, vector<vector<double> > &results, vector<vector<vector<double> > > &allResults, int& itr);
  double fitEpidemics(vector<double> params);
  double optimise_single(vector<double> &parameters, vector<vector<double> > &results, vector<vector<vector<double> > > &allResults, int& itr);
  void realtime_fit_single(double targetRSq, EpiType _epi);
  vector<vector<double> > ode_solve(vector<double> params);
  vector<vector<vector<double> > > ode_solve_separate(vector<double> params);
  double dpois(vector<vector<double> > model, vector<vector<double> > data);
  double poisson_pmf(double k, double lambda);
  double fitEpidemicsMLE(vector<double> params);


  // Parameter functions
  vector<double> rand_params(EpiType _type);
  vector<double> generate_seed_parameters();
  void sense_check(vector<vector<double> > &model);
  bool params_check(vector<double> pars, EpiType epi);

  // Data handling functions
  vector<vector<double> > return_data(){
    return(current_data);
  }
  vector<vector<double> > combine_vectors(vector<vector<double> > model, vector<vector<double> > data);
  vector<vector<double> > create_empty_data_vector(int _rows);
  void remove_epidemic(Epidemic* remove);
  vector<Epidemic*> fewer_epidemics(int j);
  Epidemic* new_epidemic(EpiType _newEpidemic, int time, double infected);
  vector<double> concatenate_vectors(vector<double> a, vector<double> b);
  double import_data(const char* file);
  void print_vector(vector< vector<double> > my_data);
  void print_vector(vector< vector<int> > my_data);
  void update_options(bool mle, bool useT0, bool useI0, bool _singleEpi, bool savePlot, bool saveResults,string location);
  void print_epidemic_type(EpiType epi);

  // Maths functions
  double aicc(double sse, int n, int k);
  vector<vector<double> > base_model(vector<vector<double> > data);
  double SStot(vector<vector<double> > data, int column);
  double calculate_mean(vector<vector<double> > data, int column);
  double calculate_sd(vector<vector<double> > data, int column);
  vector<vector<double> > get_residuals(vector<vector<double> > data1, vector<vector<double> > data2, int index);
  
  // Graph plotting functions
  void plotGraph(vector<vector<double> > finalResults, vector<vector<double> > data, int index);
  void plotGraphMulti(vector<vector<vector<double> > > finalResults, vector<vector<double> > totalResults, vector<vector<double> > data, int index, vector<double> parameters, double _RSquare, int column);
  string create_label(Epidemic* epi, vector<double> parameters, int& i);

  // Old functions
  void testAddition(vector<vector<double> > data1, vector<vector<double> > data2, double offset);
  void test_detect(vector<vector<double> > &results, vector<double> &params);
  void realtime_fit(vector<vector<double> > &results, vector<double> &params, int version);
  void reset_epidemics();
  bool check_already_tested(vector<EpiType> tested, EpiType toCheck);
  vector<Epidemic*> copy_epidemics(vector<Epidemic*> epi);
  double calculate_SSE(vector<vector<double> > data1, vector<vector<double> > data2);
};


#endif
