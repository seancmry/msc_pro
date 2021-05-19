#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>

#include "mpi.h"
#include "funcs.h"

#define N 2

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
	first->size = 10;
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



void list_mpi(list_a_t *first, list_b_t *second, MPI_Comm cart_comm){
	
	int i, j;
	double **b = matrix_new(first->size, first->dim);
	srand(time(NULL));
    	int rank, nproc = N;
	const int tag = 13;
	const int dest = 1, src = 0;
	MPI_Status stat;
	MPI_Datatype strip;
	double x = DBL_MAX;
	int strip_size;
	double *stripdata, **strip_m;


	//Parallel setup	
	MPI_Comm_rank(cart_comm, &rank);
	
	if (rank == src){
		
		//Calculate strip size
		strip_size = first->dim / nproc;
		
		//Bcast row and column size as well as strip size
		MPI_Bcast(&first->size,1,MPI_INT,0,cart_comm);
		MPI_Bcast(&first->dim,1,MPI_INT,0,cart_comm);
		MPI_Bcast(&strip_size,1,MPI_INT,0,cart_comm);	

		//define datatype for submatrix
		MPI_Type_vector(strip_size,first->size,first->dim,MPI_DOUBLE,&strip);
		MPI_Type_commit(&strip);
	
		stripdata = (double *)malloc(sizeof(double)*strip_size*first->size);
		strip_m = (double **)malloc(sizeof(double*)*strip_size);	
		for (i=0;i<strip_size;i++){
			strip_m[i] = &(stripdata[i*first->size]);
		}

		//Init second list
		list_b_t *second = malloc(sizeof(list_b_t));
		// allocate memory for the best position buffer
    		second->solution = (double *)malloc(first->dim * sizeof(double));
		second->error = x;
		MPI_Bcast(&second->error,1,MPI_DOUBLE,0,cart_comm);		

		for(i=0;i<first->size;i++){
			for(j=0;j<first->dim;j++){
				b[i][j] = rand() % 20;
				printf("%f ", b[i][j]);
			}
			printf("\n");
			memmove((void *)second->solution, (void *)b[i], sizeof(double) * first->dim);
		
			MPI_Send(second->solution,first->dim,MPI_DOUBLE,dest,tag,cart_comm);
			//MPI_Send(&second->error, first->dim, MPI_DOUBLE,dest,tag,MPI_COMM_WORLD);
		}
		free(second->solution);
		free(second);
	}
	else if (rank == dest) {
		list_b_t *second = malloc(sizeof(list_b_t));
		second->error = x;
		MPI_Bcast(&second->error,1,MPI_DOUBLE,0,cart_comm);		

		second->solution = (double *)malloc(first->dim * sizeof(double));
		MPI_Recv(second->solution,first->dim,MPI_DOUBLE,src,tag,cart_comm, &stat);
		//MPI_Recv(&second->error,first->dim,MPI_DOUBLE,src,tag,MPI_COMM_WORLD, &stat);
		
		for (i=0;i<first->dim;++i){
			printf("Received solution %f from rank %d\n",second->solution[i],rank);
		} 	
	}	
	
	MPI_Type_free(&strip);
	free(strip_m);
	free(stripdata);
	matrix_free(b,first->size);
}

