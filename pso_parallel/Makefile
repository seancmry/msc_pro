
CC = mpicc

CFLAGS = -Wall -g

LDFLAGS = -lm

DEPS = pso.h

SOURCE_FILES = $(shell find . -name '*.c')

OBJECTS = $(SOURCE_FILES:.c=.o) 

TARGET = pso.o

all: %.c $(DEPS)
	$(CC) -c -o $(TARGET) $< $(CFLAGS)
	
demopso: $(OBJECTS)
	$(CC) $^ -o $@ $(LDFLAGS)

test: all
	./bash.sh

.PHONY: clean

clean:
	rm -f demopso $(OBJECTS)
