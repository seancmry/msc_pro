
#include <stdio.h>

#ifndef FUNCS_H_
#define FUNCS_H_

typedef struct{
	int dim;
	int size;
}list_a_t;


typedef struct{
	double *solution;
	double error;
}list_b_t;


list_a_t *list_new(int dim);
void list(list_a_t *first, list_b_t *second);

#endif
