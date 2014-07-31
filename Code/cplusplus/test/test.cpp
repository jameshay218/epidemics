#include <iostream>
#include <vector>

using namespace std;

class A {
public:
  void methA(){
    cout << "Function called A" << endl;
  }
};

class C {
private:
  int i;
public:
  C(int _i){
    i = _i;
  }
  int methC(int j) {
    cout << "Function called C" << endl;
    return(i*j);
  }
};

class B {
public:
  template<class OP, class D>
  D methB(D(OP::*function)(), OP& obj){
    (obj.*function)();
  }
  template<class OP, class D, class X>
  D methB(D(OP::*function)(X), OP& obj, X x){
    (obj.*function)(x);
  }
};

int main() {
  A a;
  C c(4);
  B b;
  int l = 3;
  b.methB(&A::methA, a);
  int i = b.methB(&C::methC, c, l);
  cout << i << endl;
  return(0);
}
