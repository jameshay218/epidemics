#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

#include "sir.hpp"

int import_data(const char* file, vector<vector<int> > & my_data){
  ifstream in_stream;
  int temp1;
  int temp2;
    
  cout << "Importing file " << file << endl;

  in_stream.open(file);
  if (in_stream.fail()) {
    cout << "Error opening file" << endl;
    return 1;
  }


  while(!in_stream.eof()){
    vector<int> row;
    in_stream >> temp1;
    row.push_back(temp1);
    in_stream >> temp2;
    row.push_back(temp2);
    my_data.push_back(row);
  }
  
  for(vector< vector<int> >::const_iterator i = my_data.begin(); i !=my_data.end(); ++i){
    for(vector<int>::const_iterator j = i->begin(); j!=i->end();++j){
      cout << *j << ' ';
    }
    cout << endl;  
  }
  
  return 0;

}
