CXXFLAGS = -Wall -g -std=c++11
CXX = g++

main: main.o sir.o datahandler.o
	$(CXX) $(CXXFLAGS) main.o sir.o datahandler.o -o main -lboost_system -lboost_iostreams -lboost_filesystem

main.o: main.cpp sir.hpp datahandler.hpp simplex.hpp
	$(CXX) $(CXXFLAGS) -c main.cpp -o main.o

sir.o: sir.cpp sir.hpp
	$(CXX) $(CXXFLAGS) -c sir.cpp -o sir.o

datahandler.o: datahandler.cpp datahandler.hpp
	$(CXX) $(CXXFLAGS) -c datahandler.cpp -o datahandler.o

clean:
	rm -f *.o
