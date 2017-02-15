all : main

main : main.o roadnetwork.o
	g++ -o main main.o

main.o : main.cpp
	g++ -c main.cpp

roadnetwork.o : roadnetwork.cpp
	g++ -c roadnetwork.cpp

clean : 
	rm *.o
