
using namespace std;

#include "sir.hpp"

SIR::SIR(double parameters[7], vector<vector<double> > x){

  beta = parameters[0];
  gamma =parameters[1];
  step = parameters[2];

  S = parameters[3];
  I = parameters[4];
  R = parameters[5];

  tmax = parameters[6];
  
  current_data = x;
}

SIR::~SIR(){
  cout <<"Delete SIR"<<endl;
}


void SIR::Diff(double Pop[3]) {
  // The differential equations
  dPop[0] = - beta*Pop[0]*Pop[1];              // dS/dt
  dPop[1] = beta*Pop[0]*Pop[1] - gamma*Pop[1];   // dI/dt
  dPop[2] = gamma*Pop[1];                    // dR/dt
}

void SIR::Runge_Kutta(){
  int i;
  double dPop1[3], dPop2[3], dPop3[3], dPop4[3];
  double tmpPop[3], initialPop[3];

  /* Integrates the equations one step, using Runge-Kutta 4
     Note: we work with arrays rather than variables to make the
     coding easier */

  initialPop[0]=S; initialPop[1]=I; initialPop[2]=R;

  Diff(initialPop);
  for(i=0;i<3;i++){
    dPop1[i]=dPop[i];
    tmpPop[i]=initialPop[i]+step*dPop1[i]/2;
  }

  Diff(tmpPop);
  for(i=0;i<3;i++){
    dPop2[i]=dPop[i];
    tmpPop[i]=initialPop[i]+step*dPop2[i]/2;
  }

  Diff(tmpPop);
  for(i=0;i<3;i++){
    dPop3[i]=dPop[i];
    tmpPop[i]=initialPop[i]+step*dPop3[i];
  }
  
  Diff(tmpPop);
  
  for(i=0;i<3;i++){
    dPop4[i]=dPop[i];
    tmpPop[i]=initialPop[i]+(dPop1[i]/6 + dPop2[i]/3 + dPop3[i]/3 + dPop4[i]/6)*step;
  }


  S=tmpPop[0]; I=tmpPop[1]; R=tmpPop[2];


}

void SIR::Solve_Eq(vector<vector<double> >& _results){
  t=0;
  int i=0;
  do{
    Runge_Kutta();
    _results[i][0] = t;
    _results[i][1] = S;
    _results[i][2] = I;
    _results[i][3] = R;
    t+=step;
    i++;
  }
  while(t<tmax);
  }

double SIR::sse_sir(double parameters[3]){
  vector<vector<double> > tempData;
  tempData.resize(current_data.size());
  for(unsigned int i=0;i<current_data.size();++i){
    tempData[i].resize(4);
  }
  double sse = 0;
  beta = exp(parameters[0]);
  gamma = exp(parameters[1]);
  S = exp(parameters[2]);
  I = 1.0;
  R = 0.0;
  Solve_Eq(tempData);
  sse = add_arrays(current_data,tempData);
  return(sse);
}

vector<vector<double> > SIR::sse_sir2(double parameters[3]){
  vector<vector<double> > tempData;
  tempData.resize(current_data.size());
  for(unsigned int i=0;i<current_data.size();++i){
    tempData[i].resize(4);
  }
  double sse = 0;
  beta = exp(parameters[0]);
  gamma = exp(parameters[1]);
  S = exp(parameters[2]);
  I = 1.0;
  R = 0.0;
  Solve_Eq(tempData);
  sse = add_arrays(current_data,tempData);
  return(tempData);
}

double SIR::operator()(vector<double> param){
  vector<vector<double> > tempData;
  tempData.resize(current_data.size());
  for(unsigned int i=0;i<current_data.size();++i){
    tempData[i].resize(4);
  }
  double sse = 0;
  beta = exp(param[0]);
  gamma = exp(param[1]);
  S = exp(param[2]);
  I = 1.0;
  R = 0.0;
  Solve_Eq(tempData);
  sse = add_arrays(current_data,tempData);
  return(sse);
}


double SIR::add_arrays(vector<vector<double> > data1, vector<vector<double> > data2){
  unsigned int i = 0;
  unsigned int j = 0;
  double sse = 0;
  while(i < data1.size()-1 && j < data2.size()-1){
    //cout << "I: " << i << endl;
    //cout << "J: " << j << endl;
    
    if(data1[i][0] == data2[j][0]){
      sse += pow((data1[i][2] - data2[j][2]),2.0);
      i++;
      j++;
    }
    else if(data1[i][0] < data2[j][0]){
      sse += pow(data1[i][2],2.0);
      i++;
    }
    else{
      sse += pow(data2[i][2],2.0);
      j++;
    }
  }
  //cout << "SSE: " << sse << endl;
  return(sse);
}

void SIR::update_params(double parameters[7]){
  beta = parameters[0];
  gamma = parameters[1];
  step = parameters[2];
  S = parameters[3];
  I = parameters[4];
  R = parameters[5];
  tmax = parameters[6];
}
  
    
void SIR::user_params(){
  double parameters[7];
  cout << "Enter your starting parameters (beta, gamma, step size, S0, I0, R0 and max time): " << endl;
  for(int i=0;i<7;i++){
    cin >> parameters[i];
  }
  update_params(parameters);
}

void SIR::update_data(vector<vector<double> > x){
  current_data = x;
}

vector<double> SIR::rand_params(){
  vector<double> params;
  double beta,gamma,s0;
  srand(clock());
  setprecision(9);
  beta = (rand()%100+1)/10000.0;
  gamma = (rand()%100+beta)/1000.0;
  s0 = rand()%1000+100;
  params.push_back(log(beta));
  params.push_back(log(gamma));
  params.push_back(log(s0));
  return(params);
}

