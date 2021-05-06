#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>

#include "mpi.h"
#include "funcs.h"


//Allocate matrices
double **matrix_new(int size, int dim){
	int i;
	double **m = (double **)malloc(size *sizeof(double *));
	for (i = 0; i<size; i++) {
		m[i] = (double *)malloc(dim * sizeof(double));
	}
	return m;
}

void matrix_free(double **m, int size){ 
	int i;
	for(i=0; i<size; i++) {
		free(m[i]);
	}
	free(m);
}




list_a_t *list_new(int dim) {

	list_a_t *first = (list_a_t *)malloc(sizeof(list_a_t));
	if (first == NULL) {return NULL;}
	
  	// set some default values
  	first->dim = dim;
	first->size = 20;
	return first;
}	




void list(list_a_t *first, list_b_t *second){

	int i,j;
	double **a = matrix_new(first->size, first->dim);
	
	srand(time(NULL));

	second->error = DBL_MAX;

	//Serial	
	for (i=0;i<first->size;i++){
		for(j=0;j<first->dim;j++){
			a[i][j] = rand() % 20;
		}
	
		//Ordinary memmove - from b to solution vector, the space for which is allocated in main:
		memmove((void *)second->solution, (void *)a[i], sizeof(double) * first->dim);
		printf("Solution: %f\n", second->solution[i]);	
	}
	
	matrix_free(a, first->size);
}


void list_mpi(list_a_t *first, list_b_t *second){
	int i, j;
	double **b = matrix_new(first->size, first->dim);
	int rank;

	srand(time(NULL));
	second->error = DBL_MAX;

	//Parallel
	if (rank == 0){
		//Parallel	
		for (i=0;i<first->size;i++){
			for(j=0;j<first->dim;j++){
				b[i][j] = rand() % 20;
			}
		}
		//DO move in parallel here somewhere with an exchange func
	}				
	matrix_free(b, first->size);
}






