CC = gcc
#CC = gcc -pg

CFLAGS = -Wall -g -std=c99
#CFLAGS = -Wall -std=c99
#PREP = scorep
#PREP_CFLAGS = -L/home/support/apps/cports/rhel-6.x86_64/gnu/papi/5.6.0/lib
LDFLAGS = -lm -lgsl -lgslcblas

OBJECTS = main.o path.o pso.o utils.o

TARGET = prog

#Serial profiling (gprof)
all: $(OBJECTS)
	$(CC) $(CFLAGS) -o $(TARGET) $^ $(LDFLAGS)

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

