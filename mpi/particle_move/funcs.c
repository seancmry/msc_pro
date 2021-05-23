#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>

#include "mpi.h"
#include "funcs.h"

#define N 4

//Calculate dims
void calc_dims(int nproc, int *ndims){
	int root = (int)sqrt(nproc);
	while(nproc % root != 0)
		root--;
	ndims[0] = nproc/root;
	ndims[1] = root;
}


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
	srand(time(NULL));
	double **b, *b_dat;
	int src = 0, dest = 1;
    	int rank;
	int nproc = N;
	const int tag = 13;
	MPI_Status stat;
	double x = DBL_MAX;
	MPI_Datatype strip;
	int strip_size;
	double **strip_m, *stripdata;
	double sum = 0.0, *sum_buf = NULL;

	//Parallel setup	
	MPI_Comm_rank(cart_comm, &rank);

	if (rank == src){	
		
		//Calculate strip size
		strip_size = first->dim / nproc;	
		
		//Init matrix
		b_dat = (double *)malloc(sizeof(double)*first->size*first->dim);
		b = (double **)malloc(first->size *sizeof(double *));
		for (int i = 0; i<first->size; i++) {
			b[i] = &(b_dat[i*first->dim]);
		}


		for(i=0;i<first->size;i++){
			for(j=0;j<first->dim;j++){
				b[i][j] = rand() % 20;
				//printf("%f ", b[i][j]);
			}
		}
	}
         	
	//TEST: Split the calculation into submatrices
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

	//Scatter to all procs
	MPI_Scatter(b_dat,1,strip, &(strip_m[0][0]),1,strip,0,cart_comm);
	MPI_Barrier(cart_comm);

	for (i=0;i<strip_size;i++){
		if(i == 0){
			printf("rank = %d\n",rank);
		}
		for(j=0;j<first->dim;j++){
			printf(" %lf  ",strip_m[i][j]);
		}
		printf("\n");
	}
	
	//Do send and recv
	if (rank == src){

		//Initialise sum buffer
		sum_buf = (double *)malloc(sizeof(double) * nproc);
	
		//Init second list
		list_b_t *second = malloc(sizeof(list_b_t));
		// allocate memory for the best position buffer
   		second->solution = (double *)malloc(first->dim * sizeof(double));
		second->error = x;
		MPI_Bcast(&second->error,1,MPI_DOUBLE,0,cart_comm);	
	
		//Move data to buffer	
		for(i=0;i<first->size;i++){
	     		memmove((void *)second->solution, (void *)b_dat, sizeof(double) * first->dim);	
			//Send
			MPI_Send(second->solution,first->dim,MPI_DOUBLE,dest,tag,cart_comm);
			//MPI_Send(&second->error, first->dim, MPI_DOUBLE,dest,tag,MPI_COMM_WORLD);
			
			//Sum for Gather
			sum += second->solution[i];
		}	
		free(second->solution);
		free(second);
	}
	else if (rank == dest) {

		//Initialise sum buffer stuff
		sum_buf = (double *)malloc(sizeof(double) * nproc);
	
		list_b_t *second = malloc(sizeof(list_b_t));
		second->error = x;
		MPI_Bcast(&second->error,1,MPI_DOUBLE,0,cart_comm);		
	
		second->solution = (double *)malloc(first->dim * sizeof(double));
		MPI_Recv(second->solution,first->dim,MPI_DOUBLE,src,tag,cart_comm, &stat);
		//MPI_Recv(&second->error,first->dim,MPI_DOUBLE,src,tag,MPI_COMM_WORLD, &stat);
		
		for (i=0;i<first->dim;++i){
			printf("Received solution %f from rank %d\n",second->solution[i],rank);
		}

		//Sum what is received
		for(i=0;i<first->size;i++){
			sum += second->solution[i];
		}		
	}
	//Gather sum
	//MPI_Allgather(&sum, 1, MPI_DOUBLE, sum_buf, 1, MPI_DOUBLE,cart_comm);	
	if(rank == src || rank == dest){
		for(i=0;i<nproc;i++) printf("This sum comes from rank %d : %f\n",rank,sum_buf[i]);
	}

	if(rank == src){
		free(sum_buf);
		MPI_Type_free(&strip);
		free(strip_m);
		free(stripdata);	
		free(b_dat);
		free(b);
	}
}

