CC = gcc

CFLAGS = -Wall -g -std=c99

LDFLAGS = -lm -lgsl -lgslcblas

OBJECTS = main.o pso.o utils.o

TARGET = prog

all: $(OBJECTS)
	#mkdir -p results timings
	$(CC) $(CFLAGS) -o $(TARGET) $^ $(LDFLAGS)

test: all
	./test.sh

.PHONY: test clean

clean: 
	$(RM) $(OBJECTS) $(TARGET) 

