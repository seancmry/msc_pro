
CC = gcc
#CC = gcc -pg

CFLAGS = -Wall -g -std=c99
LDFLAGS = -lm -lgsl -lgslcblas

CFLAGS=-I/home/support/spack/spack-0.16.1/spack/opt/spack/linux-scientific7-x86_64/gcc-9.3.0/gsl-2.5-c2jmkgwsye3xruhcuieuw5f7vzmsqi74/include -I/home/support/apps/apps/openmpi/rhel-7/3.1.6-gcc9.3.0/include
LDFLAGS+=-L/home/support/spack/spack-0.16.1/spack/opt/spack/linux-scientific7-x86_64/gcc-9.3.0/gsl-2.5-c2jmkgwsye3xruhcuieuw5f7vzmsqi74/lib

#Uncomment for demo
#OBJECTS = main.o pso.o utils.o
OBJECTS = main.o utils.o path.o pso.o

TARGET = prog

#Serial profiling (gprof)
all: $(OBJECTS)
	$(CC) $(CFLAGS) -o $(TARGET) $^ $(LDFLAGS)
		
test: all
#	./path.sh
#	./demo.sh

# .PHONY: test clean

clean: 
	$(RM) $(OBJECTS) $(TARGET) *.dat *.out 
