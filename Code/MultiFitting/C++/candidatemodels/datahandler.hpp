#ifndef DATAHANDLER_HPP
#define DATAHANDLER_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iterator>
#include "../gnuplot-iostream/gnuplot-iostream.h"

#include "epidemic.hpp"

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
  vector<double> tempParams;
  vector<int> epidemicSizes, tempSizes;
  vector<vector<double> > baseModel, current_data, current_model, temp_model, empty_model, temp_data, future_data;
  vector<Epidemic*> epidemics, tempEpidemics;
  bool useMLE, optimT0, optimI0, singleEpi, plot, save,useLogistic;
  bool activeBool = false;
  string saveLocation = "graphs/";
  double betaupper, betalower, gammalower, gammaupper;
  double betaSeirUp,betaSeirDown,gammaSeirDown,gammaSeirUp,alphaSeirUp,alphaSeirDown;
  double betaIrsirUp,betaIrsirDown,gammaIrsirDown,gammaIrsirUp;
  double betaSerirUp,betaSerirDown,gammaSerirDown,gammaSerirUp,alphaSerirUp,alphaSerirDown,muSerirUp,muSerirDown;
  double gammaExpUp,gammaExpDown;
public:
  Handler();
  ~Handler();
 

  // Optimisation Functions
  vector<vector<double> > combine_vectors(vector<vector<double> > model, vector<vector<double> > data);
  void likelihood_test(vector<double> &params);
  void realtime_fit_multi(double targetRSp);
  EpiType check_epidemic(vector<vector<double> > residuals);
  double optimiseEpidemics(vector<double> &parameters, vector<vector<double> > &results, vector<vector<vector<double> > > &allResults, int& itr);
  double fitEpidemics(vector<double> params);
  vector<vector<double> > ode_solve(vector<double> params);
  vector<vector<vector<double> > > ode_solve_separate(vector<double> params);
  double dpois(vector<vector<double> > model, vector<vector<double> > data);
  double fitEpidemicsMLE(vector<double> params);

  
  // Data handling functions
  vector<vector<double> > return_data(){
    return(current_data);
  }
  
  double import_data(const char* file);
  void update_options(bool mle, bool useT0, bool useI0, bool _singleEpi, bool savePlot, bool saveResults,string location,bool useLogist);
  void save_all_results(vector<vector<double> > saveResults, vector<EpiType> selectedTypes);
  vector<vector<double> > create_empty_data_vector(int _rows);


  // Display functions
  void print_epidemic_type(EpiType epi);
  void print_vector(vector< vector<double> > my_data);
  void print_vector(vector< vector<int> > my_data);



  // Maths functions
  double calculate_median(vector<vector<double> > data);
  double calculate_MAD(vector<vector<double> > data, double median);
  vector<vector<double> > base_model(vector<vector<double> > data);
  double SStot(vector<vector<double> > data, int column);
  double calculate_mean(vector<vector<double> > data, int column);
  double calculate_sd(vector<vector<double> > data, int column);
  vector<vector<double> > get_residuals(vector<vector<double> > data1, vector<vector<double> > data2, int index);
  double aicc(double sse, int n, int k);
  double poisson_pmf(const double k, const double lambda);


  // Parameter and Epidemic handling functions
  vector<double> convert_all_params_back(vector<double> pars);
  void sense_check(vector<vector<double> > &model);
  double logistic(double x, double xmin, double xmax);
  double logit(double y, double xmin, double xmax);
  void remove_epidemic(Epidemic* remove);
  vector<Epidemic*> fewer_epidemics(int j);
  Epidemic* new_epidemic(EpiType _newEpidemic, int time);
  vector<double> concatenate_vectors(vector<double> a, vector<double> b);
  vector<double> generate_seed_parameters();
  void deactivate_epidemics();
  void activate_last_epidemic();
  vector<double> rand_params(EpiType _type);
  void update_latest_epidemic(vector<double> finalParams);
  void initial_param_bound();
  vector<int> fewer_sizes(int j);
  vector<double> convert_parameters_back(vector<double> pars, int number, double detection, double lowerTime, double upperTime, EpiType _type);
  vector<double> convert_parameters_forward(vector<double> pars, int number, double time,double lowerTime, double upperTime, EpiType _type);


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
