#include <vector>
namespace Measurement{
 
   void GetData(std::vector<std::vector<double> > _data, std::vector<double> params);
   void PrintParams();
   void Diff(double Pop[3]);
   void Runge_Kutta();
   void Solve_Eq(std::vector<std::vector<double> >& data);
   double sse_sir(std::vector<double> parameters);
   double add_arrays(std::vector<std::vector<double> > data1, std::vector<std::vector<double> > data2);


}
