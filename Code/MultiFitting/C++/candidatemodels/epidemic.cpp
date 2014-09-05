
using namespace std;

#include "epidemic.hpp"

Epidemic::~Epidemic(){
}

// Constructor
Epidemic::Epidemic(double _tmax, vector<vector<double> > x, EpiType _type, int detection){
  // Determine the type of epidemic and adjust population size stores accordingly
  if(_type==sir || _type==irsir){
    noPops=3;
    parSize= 4;
  }
  else if(_type==spike){
    noPops=1;
    parSize= 3;
  }
  else if (_type==seir){
    noPops = 4;
    parSize= 5;
  }
  else if (_type == serir){
    noPops = 5;
    parSize= 6;
  }
  else{
    noPops=3;    
    parSize= 4;
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
  
  // Step size for Runge Kutta. Decrease to increase accuracy, increase to improve speed.
  step = 0.5;
  
  // Desired size of total solved model
  tmax = _tmax;

  // Set current data
  current_data = x;
  temp_model.resize(x.size()/step);

  // Create a vector with the appropriate first column (ie. time values)
  for(unsigned int i=0;i<(x.size()/step);++i){
    temp_model[i].resize(2);
    fill(temp_model[i].begin(),temp_model[i].end(),0.0);
    temp_model[i][0] = i;
  }
  total_model = temp_model;
  type = _type;
  detectionTime = detection;
  detectionTimeTemp = detection;
  active = true;
  optimTime = detection;
  //if(minTime <=0) minTime = 0;
  //  if(seedTime <=0) seedTime = 0;

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
void Epidemic::Solve_Eq_t0(vector<vector<double> >& _results, int index){
  t=0;
  int i=0;
  _results[i][0] = t + t0;
  _results[i][1] = populations[index];
  t+=step;
  i++;
  do{
    Runge_Kutta();
    _results[i][0] = t + t0;
    _results[i][1] = populations[index];
    t+=step; // Step size
    i++;
  }
  while((current_data.size()-t) > step);
}

/* Solves the ODEs for the current SIR parameters and saves these to _results. Difference
   to above is that this is carried out for the entire data range (tmax) */
void Epidemic::Solve_Eq_total(vector<vector<double> >& _results, int index){
  t=0.0;
  int i=0;
  _results[i][0] = t + t0;
  _results[i][1] = populations[index];
  t+=step;
  i++;
  do{
    Runge_Kutta();
    _results[i][0] = t + t0;
    _results[i][1] = populations[index];
    t+=step; // Step size
    i++;
  }
  while((tmax - t) > step);
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
    temp_model[i].resize(2);
    fill(temp_model[i].begin(),temp_model[i].end(),0.0);
    temp_model[i][0] = i*step;
  }
  total_model= temp_model;
}



// Updates SIR data
void Epidemic::update_data(vector<vector<double> > x){
  current_data.clear();
  current_data = x;
}
