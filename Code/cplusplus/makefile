CXXFLAGS = -Wall -g
CXX = g++

main: main.o sir.o datahandler.o data.o
	$(CXX) $(CXXFLAGS) main.o sir.o datahandler.o data.o -o main -lboost_system -lboost_iostreams -lboost_filesystem

main.o: main.cpp sir.hpp datahandler.hpp
	$(CXX) $(CXXFLAGS) -c main.cpp -o main.o

sir.o: sir.cpp sir.hpp
	$(CXX) $(CXXFLAGS) -c sir.cpp -o sir.o

data.o: data.cpp data.h
	$(CXX) $(CXXFLAGS) -c data.cpp

datahandler.o: datahandler.cpp datahandler.hpp
	$(CXX) $(CXXFLAGS) -c datahandler.cpp -o datahandler.o

clean:
	rm -f *.o
