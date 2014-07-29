/* Solve simple SIR Model
   Author: Jiansen Lu
   Date: June 25, 2010
   Description:
   dS/dt = - beta*S*I
   dI/dt = beta*S - gamma*I
   dR/dt = gamma*I
   where S: Susceptibel
   I: Infectious
   R: Recovered
   beta: transmission rate
   gamma: recover rate:
   set initial condition S(0)=1e-6, I(0)=1.0-S(0)
   beta = 1.5, gamma = 0.2
   Algorithm: using Runge-Kutta method
   dy/dt = f(t,y), y(t_0) = y_0
   t_(n+1) = t_n+h
   y_(n+1) = y_n + (k_1+2*k_2+2*_k3+k_4)/6
   where: k_1=f(t_n,y_n), k_2=f(t_n+h/2,y_n+h*k_1/2)
   k_3=f(t_n+h/2,y_n+h*k_2/2), k_4=f(t_n+h,y_n+h*k_3)
*/

#include<iostream>
#include<vector>
#include<cmath>
#include "gnuplot-iostream/gnuplot-iostream.h"

using namespace std;

class SIR{
private:
  
  double t, step;
  double S,I,R,Pop[3];
  double dPop[3];

  double beta, gamma;
  double tmax;
public:


  SIR();
  SIR(double beta0, double gamma0, double step0, double S00, \
      double I00, double R00, double tmax0);
  ~SIR();

  void Diff(double Pop[3]);
  void Runge_Kutta();

  void Solve_Eq(double _results[][3]);
};
// Initialise the equations and Runge-Kutta integration
SIR::SIR(double beta0, double gamma0, double step0,double S00, \
	 double I00, double R00, double tmax0)
{

  beta = beta0;
  gamma =gamma0;
  step = step0;

  S = S00;
  I = I00;
  R = R00;

  tmax = tmax0;
}
SIR::~SIR(){
  cout <<"delete SIR"<<endl;
}

void SIR::Diff(double Pop[3])
{

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
  for(i=0;i<3;i++)
    {

      dPop1[i]=dPop[i];
      tmpPop[i]=initialPop[i]+step*dPop1[i]/2;
    }

  Diff(tmpPop);
  for(i=0;i<3;i++)
    {


      dPop2[i]=dPop[i];
      tmpPop[i]=initialPop[i]+step*dPop2[i]/2;
    }

  Diff(tmpPop);
  for(i=0;i<3;i++)
    {

      dPop3[i]=dPop[i];
      tmpPop[i]=initialPop[i]+step*dPop3[i];
    }

  Diff(tmpPop);

  for(i=0;i<3;i++)
    {

      dPop4[i]=dPop[i];
      tmpPop[i]=initialPop[i]+(dPop1[i]/6 + dPop2[i]/3 + dPop3[i]/3 + dPop4[i]/6)*step;
    }


  S=tmpPop[0]; I=tmpPop[1]; R=tmpPop[2];


}

void SIR::Solve_Eq(double _results[][3]){
  t=0;
  int i=0;
  cout <<"t    S    I       R"<<endl;
  do
    {
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
int main(int argc, char** argv)
{
  Gnuplot gp;
 
  double beta0 = 0.001;
  double gamma0 =0.1;

  int I00 = 1;
  int S00 = 500;
  int R00 = 0;
  double tmax0 = 50;

  /* Find a suitable time-scale for outputs */
  double step0 = 0.01;
  int dataSize = int(tmax0/step0);
  
  double results[dataSize][3];

  SIR mySIR(beta0, gamma0,step0,S00,  I00, R00,  tmax0);

  mySIR.Solve_Eq(results);
  
  vector<vector<double> > results1;

  for(int i = 0; i<int(tmax0/step0);++i){
    vector<double> row;
    row.push_back(results[i][0]);
    row.push_back(results[i][2]);
    results1.push_back(row);
    }
  /*
  for(vector< vector<double> >::const_iterator i = results1.begin(); i !=results1.end(); ++i){
    for(vector<double>::const_iterator j = i->begin(); j!=i->end();++j){
      cout << *j << ' ';
    }
    cout << endl;  
    }
  */
  gp << "set output 'my_graph_1.png'\n";
  gp << "plot '-' with lines\n";
  gp.send1d(results1);
  
  
  return(0);

}
