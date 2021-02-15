CC = gcc

#CFLAGS = -Wall -g -std=c99
CFLAGS = -Wall -std=c99

LDFLAGS = -lm -lgsl -lgslcblas

OBJECTS = main.o path.o pso.o utils.o

TARGET = prog

all: $(OBJECTS)
	$(CC) $(CFLAGS) -o $(TARGET) $^ $(LDFLAGS)

test: all
	./test.sh

.PHONY: test clean

clean: 
	$(RM) $(OBJECTS) $(TARGET) 

