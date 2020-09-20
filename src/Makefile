
CC = gcc


CFLAGS = -Wall -g
DEPS = pso.h  

SOURCE_FILES = *.c
OBJ_FILES = *.o

LDFLAGS = -L/usr/local/lib -lgsl -lgslcblas -lm 
INCDIRS = -I/usr/local/include

all: 
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(INCDIRS)
	
demopso: $(OBJ_FILES)
	$(CC) $(CFLAGS) -o $@ $(LDFLAGS) $(INCDIRS)

test: all
	./demopso

.PHONY: clean

clean:
	rm -f demopso $(OBJ_FILES)


