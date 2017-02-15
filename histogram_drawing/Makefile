all : getResults makeHistogram


getResults : getResults.o
	g++ -o getResults getResults.o
makeHistogram : makeHistogram.o
	g++ -o makeHistogram makeHistogram.o

getResults.o : getResults.cpp
	g++ -c getResults.cpp
makeHistogram.o : makeHistogram.cpp
	g++ -c makeHistogram.cpp
clean : 
	rm *.o
