CXXFLAGS = -Wall -g -std=c++11  -fexceptions
CXX = g++

main: main.o epidemic.o sir.o seir.o spike.o datahandler.o
	$(CXX) $(CXXFLAGS) main.o epidemic.o sir.o seir.o spike.o datahandler.o -o main -lboost_system -lboost_iostreams -lboost_filesystem -lgsl -lgslcblas -lm

main.o: main.cpp epidemic.hpp sir.hpp seir.hpp spike.hpp datahandler.hpp simplex.hpp
	$(CXX) $(CXXFLAGS) -c main.cpp -o main.o

epidemic.o: epidemic.cpp epidemic.hpp
	$(CXX) $(CXXFLAGS) -c epidemic.cpp -o epidemic.o 


sir.o: sir.cpp sir.hpp epidemic.hpp
	$(CXX) $(CXXFLAGS) -c sir.cpp -o sir.o

seir.o: seir.cpp seir.hpp epidemic.hpp
	$(CXX) $(CXXFLAGS) -c seir.cpp -o seir.o	

spike.o: spike.cpp spike.hpp epidemic.hpp
	$(CXX) $(CXXFLAGS) -c spike.cpp -o spike.o		


datahandler.o: datahandler.cpp datahandler.hpp
	$(CXX) $(CXXFLAGS) -c datahandler.cpp -o datahandler.o

clean:
	rm -f *.o
