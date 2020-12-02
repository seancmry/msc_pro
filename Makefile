CC = mpicc

CFLAGS = -W -Wall -g

LDFLAGS = -lm -lgsl -lgslcblas

OBJECTS: main.o path.o pso.o utils.o

SOURCES: main.c 

TARGET = prog

all: $(OBJECTS)
	#mkdir -p results timings
	$(CC) $(SOURCES) $(CFLAGS) -o $(TARGET) $(LDFLAGS)

test: all
	./test.sh

.PHONY: clean

clean: 
	$(RM) $(OBJECTS) $(TARGET) *.csv a.out

