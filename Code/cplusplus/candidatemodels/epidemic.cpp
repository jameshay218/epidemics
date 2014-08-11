
using namespace std;

#include "epidemic.hpp"

Epidemic::~Epidemic(){
}

Epidemic::Epidemic(double _tmax, vector<vector<double> > x, EpiType _type){
  if(_type==sir || _type==irsir){
    noPops=3;
  }
  else if(_type==spike){
    noPops=1;
  }
  else{
    noPops=4;
  }
  dPop.resize(noPops);
  dPop1.resize(noPops);
  dPop2.resize(noPops);
  dPop3.resize(noPops);
  dPop4.resize(noPops);
  tmpPop.resize(noPops);
  initialPop.resize(noPops);
  populations.resize(noPops);
  diffIndex = 0;
  step = 1.0;
  tmax = _tmax;
  current_data = x;
  temp_model.resize(x.size());
  // Create a vector with the appropriate first column (ie. time values)
  for(unsigned int i=0;i<x.size();++i){
    temp_model[i].resize(4);
    fill(temp_model[i].begin(),temp_model[i].end(),0.0);
    temp_model[i][0] = i;
  }
  total_model = temp_model;
}

/* ============================= RUNGE KUTTA INTEGRATION ALGORITHM =========================*/
void Epidemic::Runge_Kutta(){
  /* Integrates the equations one step, using Runge-Kutta 4
     Note: we work with arrays rather than variables to make the
     coding easier */
  diffIndex = 0;
  for(diffIndex=0;diffIndex<noPops;++diffIndex){
    initialPop[diffIndex] = populations[diffIndex];
  }

  Diff(initialPop);

  for(diffIndex=0;diffIndex<noPops;++diffIndex){
    dPop1[diffIndex]=dPop[diffIndex];
    tmpPop[diffIndex]=initialPop[diffIndex]+step*dPop1[diffIndex]/2;
  }

  Diff(tmpPop);
  for(diffIndex=0;diffIndex<noPops;++diffIndex){
    dPop2[diffIndex]=dPop[diffIndex];
    tmpPop[diffIndex]=initialPop[diffIndex]+step*dPop2[diffIndex]/2;
  }

  Diff(tmpPop);
  for(diffIndex=0;diffIndex<noPops;++diffIndex){
    dPop3[diffIndex]=dPop[diffIndex];
    tmpPop[diffIndex]=initialPop[diffIndex]+step*dPop3[diffIndex];
  }
  
  Diff(tmpPop);
  for(diffIndex=0;diffIndex<noPops;++diffIndex){
    dPop4[diffIndex]=dPop[diffIndex];
    tmpPop[diffIndex]=initialPop[diffIndex]+(dPop1[diffIndex]/6 + dPop2[diffIndex]/3 + dPop3[diffIndex]/3 + dPop4[diffIndex]/6)*step;
  }
  for(diffIndex=0;diffIndex<noPops;++diffIndex){
    populations[diffIndex] = tmpPop[diffIndex];
  }
}


/* Solves the set of differential equations with the current SIR parameters and
   saves these to the passed _results vector. Note that this only takes place up to 
   the current data size. */
void Epidemic::Solve_Eq_t0(vector<vector<double> >& _results){
  t=0;
  int i=0;
  int j=1;
  do{
    Runge_Kutta();
    _results[i][0] = t + int(t0);
    for(j=1;j<=noPops;++j){
      _results[i][j] = populations[j-1];
    }
    t+=step; // Step size
    i++;
  }
  while(t<current_data.size());
}

/* Solves the ODEs for the current SIR parameters and saves these to _results. Difference
   to above is that this is carried out for the entire data range (tmax) */
void Epidemic::Solve_Eq_total(vector<vector<double> >& _results){
   t=0;
  int i=0;
  int j=1;
  do{
    Runge_Kutta();
    _results[i][0] = t + int(t0);
    for(j=1;j<=noPops;++j){
      _results[i][j] = populations[j-1];
    }
    t+=step; // Step size
    i++;
  }
  while(t<tmax);
}



/* Returns the negative log likelihood of a model given a set of data */
double Epidemic::dpois(vector<vector<double> > model, vector<vector<double> > data){
  double logLikelihood;
  int N = data.size();
  logLikelihood = 0.0;
  for(int i=0;i<N;++i){
    if(model[i][2] == 0 && data[i][2] == 0){
      logLikelihood += log(1);
    }
    else{
      //logLikelihood += log(gsl_ran_poisson_pdf(model[i][2], data[i][2]));
      logLikelihood += log(poisson_pmf(model[i][2], data[i][2]));
    }
  }
  return(-logLikelihood);
}

double Epidemic::poisson_pmf(const double k, const double lambda) {
  return exp(k * log(lambda) - lgamma(k + 1.0) - lambda);
}


/* Goes through to sets of data and adds all values (apart from the first column, which is the time
   values). Returns the combined result. NOTE PASSED VECTORS HAVE 4 COLUMNS */
vector<vector<double> > Epidemic::combine_vectors(vector<vector<double> > data1, vector<vector<double> > data2){
  unsigned int i = 0;
  unsigned int j = 0;
  unsigned int k = 0;
  int l = 1;
 
  vector<vector<double> > final = data1;
  
  while(i < data1.size() && j < data2.size()){
    if(data1[i][0] == data2[j][0]){
      final[k][0] = data1[i][0];
      for(l=1;l<=noPops;++l){
	final[k][l] = data1[i][l] + data2[j][l];
      }
      i++;
      j++;
      k++;
    }
    //If first dataset is still before second, full difference
    else if(data1[i][0] < data2[j][0]){
      final[k][0] = data1[i][0];
      for(l=1;l<=noPops;++l){
	final[k][l] = data1[i][l];
      }
      i++;
      k++;
    }
    //...
    else{
      final[k][0] = data2[j][0];
      for(l=1;l<=noPops;++l){
	final[k][l] = data2[j][l];
      }
      j++;
      k++;
    }
  }
  return(final);
}


/* Calculates and returns the sum of squared errors between to sets of data. Checks 
   the first column to see if on the same time point, and finds squared difference if
   so. If one dataset is at an earlier time point, use the entire squared value as
   the squared residual.
*/
double Epidemic::calculate_SSE(vector<vector<double> > data1, vector<vector<double> > data2, int index){
  unsigned int i = 0;
  unsigned int j = 0;
  double sse = 0;
  while(i < data1.size() && j < data2.size()){
    //If indices give same time point, find sse
    if(data1[i][0] == data2[j][0]){
      sse += pow((data1[i][index] - data2[j][index]),2.0);
      i++;
      j++;
    }
    //If first dataset is still before second, full difference
    else if(data1[i][0] < data2[j][0]){
      sse += pow(data1[i][index],2.0);
      i++;
    }
    //...
    else{
      sse += pow(data2[j][index],2.0);
      j++;
    }
  }
  return(sse);
}





/* ============================== HOUSEKEEPING FUNCTIONS ================================== */

vector<double> Epidemic::return_parameters(){
  return(pars);
}

void Epidemic::reset_models(int size){
  temp_model.clear();
  total_model.clear();
  temp_model.resize(size);
  // Create a vector with the appropriate first column (ie. time values)
  for(int i=0;i<size;++i){
    temp_model[i].resize(4);
    fill(temp_model[i].begin(),temp_model[i].end(),0.0);
    temp_model[i][0] = i;
  }
  total_model= temp_model;
}



// Updates SIR data
void Epidemic::update_data(vector<vector<double> > x){
  current_data.clear();
  current_data = x;
}
