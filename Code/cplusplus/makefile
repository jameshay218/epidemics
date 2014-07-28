CXXFLAGS = -Wall -g
CXX = g++

main: main.o sir.o
	$(CXX) $(CXXFLAGS) main.o sir.o -o main -lboost_system -lboost_iostreams -lboost_filesystem

main.o: main.cpp sir.hpp
	$(CXX) $(CXXFLAGS) -c main.cpp -o main.o

sir.o: sir.cpp sir.hpp
	$(CXX) $(CXXFLAGS) -c sir.cpp -o sir.o

clean:
	rm -f *.o
