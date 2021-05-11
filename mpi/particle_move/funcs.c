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
	double **m = (double **)malloc(size *sizeof(double *));
	for (int i = 0; i<size; i++) {
		m[i] = (double *)malloc(dim * sizeof(double));
	}
	return m;
}

void matrix_free(double **m, int size){ 
	for(int i=0; i<size; i++) {
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


/*
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
		//memmove((void *)second->solution, (void *)a[i], sizeof(double) * first->dim);
		memmove_mpi((void *)second->solution, (void *)a[i], sizeof(double) * first->dim);
		printf("Solution: %f\n", second->solution[i]);	
	}
	
	matrix_free(a, first->size);
}
*/


void list_mpi(list_a_t *first, list_b_t send, list_b_t recv){
	
	int i, j;
	double **b = matrix_new(first->size, first->dim);
	srand(time(NULL));
	send.error = DBL_MAX;
    	int rank = 0;
	MPI_Status stat;
    
	/* create a type for struct car */
    	const int nitems=1;
    	int blocklengths[2] = {1,1};
    	MPI_Datatype types[1] = {MPI_DOUBLE};
    	MPI_Datatype mpi_move_type;
    	MPI_Aint offsets[1];

    	offsets[0] = offsetof(list_b_t, solution);

    	MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_move_type);
    	MPI_Type_commit(&mpi_move_type);

	//Parallel	
	for (i=0;i<first->size;i++){
		send.solution[i] = 0.0;
		for(j=0;j<first->dim;j++){
			b[i][j] = rand() % 20;
		}
		if(rank == 0){
			send.solution = b[i];
			//printf("Send.solution: %f\n", send.solution[i]);
			MPI_Send(&send, 1, mpi_move_type, 1, 0, MPI_COMM_WORLD);
			printf("Rank %d: send sendbuf\n", rank);
		}
		else if(rank == 1) {
			//recv.solution[i] = 0.0;
			//list_b_t recv;
			MPI_Recv(&recv, 1, mpi_move_type, 0, 0, MPI_COMM_WORLD, &stat);
			printf("Rank %d: recvbuf: b[i] = %f\n", rank, recv.solution[i]);
		}
	}
   
	MPI_Type_free(&mpi_move_type);
	matrix_free(b,first->size);
}

