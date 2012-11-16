CC       = g++
DEBUG    = -g -Wall
LIBFLAGS = -llapack -lblas -larmadillo
CFLAGS   = $(DEBUG) -O2 -c -fopenmp
LFLAGS   = $(DEBUG) -fopenmp

Jacobi : main.o classJacobi.o
	$(CC) $(LFLAGS) main.o classJacobi.o $(LIBFLAGS) -o $@

main.o : main.cpp
	$(CC) $(CFLAGS) main.cpp -o $@

classJacobi.o : classJacobi.cpp classJacobi.hpp
	$(CC) $(CFLAGS) classJacobi.cpp -o $@

clean:
	rm *.o Jacobi
