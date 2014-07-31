
using namespace std;

#include "simplex.hpp"

Simplex::Simplex(){
}

Simplex::~Simplex(){
}

void Simplex::print_out(){
  cout << "Hi" << endl;
}
template<class D>
vector<D> Simplex::tester(vector<D> init){
  cout << "hello" << endl;
  vector<D> hi;
  return(hi);
}



template<class D, class OP, class X>
vector<D> Simplex::neldermead(X(OP::*f)(vector<D>),   //target function of object
				   OP& obj,                //object to work from
				   vector<D> init    //initial guess of the parameters
				   ){    //iteration step number
		    

  int N=init.size();                         //space dimension
  const double a=1.0, b=1.0, g=0.5, h=0.5;   //coefficients
  //a: reflection  -> xr  
  //b: expansion   -> xe 
  //g: contraction -> xc
  //h: full contraction to x1
  vector<D> xcentroid_old(N,0);   //simplex center * (N+1)
  vector<D> xcentroid_new(N,0);   //simplex center * (N+1)
  vector<D> vf(N+1,0);            //f evaluated at simplex vertexes       
  int x1=0, xn=0, xnp1=0;         //x1:   f(x1) = min { f(x1), f(x2)...f(x_{n+1} }
  //xnp1: f(xnp1) = max { f(x1), f(x2)...f(x_{n+1} }
  //xn:   f(xn)<f(xnp1) && f(xn)> all other f(x_i)
  int cnt=0; //iteration step number

  D tol=1E8*numeric_limits<D>::epsilon(); //termination criteria
  vector<vector<D> > x =  vector<vector<D> >(); //x: The Simplex
  int iterations=1E5;    //iteration step number

  if(x.size()==0) //if no initial simplex is specified
    { //construct the trial simplex
      //based upon the initial guess parameters
      vector<D> del( init );
      transform(del.begin(), del.end(), del.begin(), 
		     bind2nd( divides<D>() , 20) );//'20' is picked 
                                                             //assuming initial trail close to true
      
      for(int i=0; i<N; ++i){
	vector<D> tmp( init );
	tmp[i] +=  del[i];
	x.push_back( tmp );
      }
      x.push_back(init);//x.size()=N+1, x[i].size()=N
      
      //xcentriod
      transform(init.begin(), init.end(), 
		     xcentroid_old.begin(), bind2nd(multiplies<D>(), N+1) );
    }//constructing the simplex finished
    
  //optimization begins
  for(cnt=0; cnt<iterations; ++cnt){

    for(int i=0;i<N+1;++i){
      vf[i]= (obj.*f)(x[i]);
    }
      
    x1=0; xn=0; xnp1=0;//find index of max, second max, min of vf.
      
    for(unsigned int i=0;i<vf.size();++i){
      if(vf[i]<vf[x1]){
	x1=i;
      }
      if(vf[i]>vf[xnp1]){
	xnp1=i;
      }
    }
      
    xn=x1;
      
    for(unsigned int i=0; i<vf.size();++i){ 
      if( vf[i]<vf[xnp1] && vf[i]>vf[xn] )
	xn=i;
    }
    //x1, xn, xnp1 are found

    vector<D> xg(N, 0);//xg: centroid of the N best vertexes
    for(unsigned int i=0; i<x.size(); ++i){
      if(i!=xnp1)
	transform(xg.begin(), xg.end(), x[i].begin(), xg.begin(), plus<D>() );
    }
    transform(xg.begin(), xg.end(), 
		   x[xnp1].begin(), xcentroid_new.begin(), plus<D>());
    transform(xg.begin(), xg.end(), xg.begin(), 
		   bind2nd(divides<D>(), N) );
    //xg found, xcentroid_new updated

    //termination condition
    D diff=0;          //calculate the difference of the simplex centers
    //see if the difference is less than the termination criteria
    for(int i=0; i<N; ++i)     
      diff += fabs(xcentroid_old[i]-xcentroid_new[i]);

    if (diff/N < tol) break;              //terminate the optimizer
    else xcentroid_old.swap(xcentroid_new); //update simplex center
      
    //reflection:
    vector<D> xr(N,0); 
    for( int i=0; i<N; ++i)
      xr[i]=xg[i]+a*(xg[i]-x[xnp1][i]);
    //reflection, xr found
      
    D fxr=(obj.*f)(xr);//record function at xr
      
    if(vf[x1]<=fxr && fxr<=vf[xn])
      copy(xr.begin(), xr.end(), x[xnp1].begin() );
      
    //expansion:
    else if(fxr<vf[x1]){
      vector<D> xe(N,0);
      for( int i=0; i<N; ++i)
	xe[i]=xr[i]+b*(xr[i]-xg[i]);
      if( (obj.*f)(xe) < fxr )
	copy(xe.begin(), xe.end(), x[xnp1].begin() );
      else
	copy(xr.begin(), xr.end(), x[xnp1].begin() );
    }//expansion finished,  xe is not used outside the scope
      
    //contraction:
    else if( fxr > vf[xn] ){
      vector<D> xc(N,0);
      for( int i=0; i<N; ++i)
	xc[i]=xg[i]+g*(x[xnp1][i]-xg[i]);
      if( (obj.*f)(xc) < vf[xnp1] )
	copy(xc.begin(), xc.end(), x[xnp1].begin() );
      else{
	  
	for(unsigned int i=0; i<x.size(); ++i ){
	  if( i!=x1 ){ 
	    for(int j=0; j<N; ++j) 
	      x[i][j] = x[x1][j] + h * ( x[i][j]-x[x1][j] );
	  }
	}
	  
      }
    }//contraction finished, xc is not used outside the scope

  }//optimization is finished

  if(cnt==iterations){//max number of iteration achieves before tol is satisfied
    cout<<"Iteration limit achieves, result may not be optimal"<<endl;
  }
  return x[x1];
}


