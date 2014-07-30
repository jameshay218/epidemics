spectrum_nlsq_fit: main.o data.o spectrum_RSS.o
	g++ main.o data.o spectrum_RSS.o -o spectrum_nlsq_fit.out
main.o: main.cpp data.h
	g++ -c main.cpp
data.o: data.cpp data.h
	g++ -c data.cpp
spectrum_RSS.o: spectrum_RSS.cpp data.h
	g++ -c spectrum_RSS.cpp