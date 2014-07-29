CXXFLAGS = -Wall -g
CXX = g++

sir: sirSolver.o 
	$(CXX) $(CXXFLAGS) sirSolver.o -o sir -lboost_system -lboost_iostreams -lboost_filesystem

sirSolver.o: sirSolver.cpp
	$(CXX) $(CXXFLAGS) -c sirSolver.cpp -o sirSolver.o

clean:
	rm -f *.o
