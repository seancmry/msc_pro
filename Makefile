
CC = gcc
#CC = gcc -pg
MPICC = mpicc

CFLAGS = -Wall -g -std=c99
MPICCFLAGS = -Wall -g
#PREP = scorep
#PREP_CFLAGS = -L/home/support/apps/cports/rhel-6.x86_64/gnu/papi/5.6.0/lib
LDFLAGS = -lm -lgsl -lgslcblas

#OBJECTS = main.o pso.o utils.o
OBJECTS = main.o utils.o path.o pso.o

TARGET = prog

#Serial profiling (gprof)
all: $(OBJECTS)
	#$(CC) $(CFLAGS) -o $(TARGET) $^ $(LDFLAGS)
	$(MPICC) $(MPICCFLAGS) -o $(TARGET) $^ $(LDFLAGS)
		
#Parallel profiling (Score-P, Vampir)
#all: $(OBJECTS)
#ifeq ($(PROFILING),no)
#	$(CC) $(CFLAGS) -o $(TARGET) $^ $(LDFLAGS)
#else
#	$(PREP) $(CC) $(CFLAGS) -o $(TARGET) $^ $(LDFLAGS) $(PREP_CFLAGS)
#endif

test: all
	./test.sh

.PHONY: test clean

clean: 
	$(RM) $(OBJECTS) $(TARGET) *.dat *.out 

