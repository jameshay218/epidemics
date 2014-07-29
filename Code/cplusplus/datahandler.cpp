using namespace std;

#include "datahandler.hpp"

Handler::Handler(){
}

Handler::~Handler(){
  cout << "Delete Data Handler" << endl;
}


int Handler::import_data(const char* file){
  ifstream in_stream;
  string temp;
  current_data.clear();
  int rows = 0;
  cout << "Importing file " << file << endl;
  current_data.clear();
  in_stream.open(file);
  if (in_stream.fail()) {
    cout << "Error opening file" << endl;
    return -1;
  }
  while(!in_stream.eof()){
    if(!getline(in_stream,temp)) break;
    istringstream ss(temp);
    vector<double> record;
    while(ss){
      double temp1;
      if(!getline(ss,temp,',')) break;
      stringstream convert(temp);
      if(!(convert >> temp1)) temp1 = 0.0;
      record.push_back(temp1);
    }
    current_data.push_back(record);
    rows++;
  }
  cout << "File imported!" << endl;
  cout << "Size: " << current_data.size() << endl;
  return rows;
}

void Handler::print_vector(vector< vector<double> > my_data){
  for(vector< vector<double> >::const_iterator i = my_data.begin(); i !=my_data.end(); ++i){
    for(vector<double>::const_iterator j = i->begin(); j!=i->end();++j){
      cout << *j << ' ';
    }
    cout << endl;  
  }
}

vector<vector<double> > Handler::array_to_vector(double** data){
  vector<vector<double> > results;
  vector<double> row;
  for(int i = 0; i < 100;++i){
    //    cout << i << endl;
    //row.clear();
    //cout << "Here?" << endl;
    cout << data[i][0] << ',' << data[i][2] << endl;
    //    row.push_back(data[i][0]);
    //cout << "Or here?" << endl;
    //row.push_back(data[i][2]);

    //results.push_back(row);
  }
  return(results);
}
  


