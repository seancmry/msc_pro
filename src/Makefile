
CC = gcc
CFLAGS = -Wall -g
DEPS = pso.h path.h utils.h

SOURCE_FILES = $(shell find . -name '*.c')
OBJ_FILES = $(SOURCE_FILES:.c=.o) 
LDFLAGS = -L/usr/local/lib -lgsl -lgslcblas -lm 
INCDIRS = -I/usr/local/include

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) $(INCDIRS) -o $@ $^ $(LDFLAGS)
	
demopso: $(OBJ_FILES)
	$(CC) $^ -o $@ $(LDFLAGS)


.PHONY: clean
clean:
	rm -f demopso $(OBJ_FILES)
