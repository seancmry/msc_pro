
#include <stdio.h>
#include "mpi.h"

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

void calc_dims(int nproc, int *ndims);
list_a_t *list_new(int dim);
void list(list_a_t *first, list_b_t *second);
void *memmove_mpi(void* dest, const void* src, unsigned int n);
void list_mpi(list_a_t *first, list_b_t *second, MPI_Comm cart_comm);

#endif
