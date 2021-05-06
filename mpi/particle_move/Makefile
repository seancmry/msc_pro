CC = mpicc

CFLAGS= -Wall -g -I/home/support/spack/spack-0.16.1/spack/opt/spack/linux-scientific7-x86_64/gcc-9.3.0/gsl-2.5-c2jmkgwsye3xruhcuieuw5f7vzmsqi74/include -I/home/support/apps/apps/openmpi/rhel-7/3.1.6-gcc9.3.0/include

LDFLAGS = -lm 

OBJECTS = main.o funcs.o

TARGET = prog

all: $(OBJECTS)
	$(CC) $(CFLAGS) -o $(TARGET) $^ $(LDFLAGS)

test: all
	./test.sh

# .PHONY: test clean

clean: 
	$(RM) $(OBJECTS) $(TARGET)  
