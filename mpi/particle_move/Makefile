CC = gcc

CLFAGS = -Wall -g

#MPICCFLAGS = -Wall -g
#PREP = scorep
#PREP_CFLAGS = -L/home/support/apps/cports/rhel-6.x86_64/gnu/papi/5.6.0/lib
#LDFLAGS = -lm -lgsl -lgslcblas
LDFLAGS = -lm

#CFLAGS=-I/home/support/spack/spack-0.16.1/spack/opt/spack/linux-scientific7-x86_64/gcc-9.3.0/gsl-2.5-c2jmkgwsye3xruhcuieuw5f7vzmsqi74/include -I/home/support/apps/apps/openmpi/rhel-7/3.1.6-gcc9.3.0/include
#CFLAGS = -I/home/support/apps/apps/openmpi/rhel-7/3.1.6-gcc9.3.0/include
#LDFLAGS+=-L/home/support/spack/spack-0.16.1/spack/opt/spack/linux-scientific7-x86_64/gcc-9.3.0/gsl-2.5-c2jmkgwsye3xruhcuieuw5f7vzmsqi74/lib

OBJECTS = main.o funcs.o

TARGET = prog

all: $(OBJECTS)
	$(CC) $(CFLAGS) -o $(TARGET) $^ $(LDFLAGS)

test: all
	./test.sh

# .PHONY: test clean

clean: 
	$(RM) $(OBJECTS) $(TARGET)  
