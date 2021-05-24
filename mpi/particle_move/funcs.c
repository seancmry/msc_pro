#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>

#include "mpi.h"
#include "funcs.h"

#define N 8

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
	double **global, **local;
    	int rank;
	int nproc = N;
	const int gridsize = N;
	const int procgridsize = 8;
	//const int tag = 13;
	//MPI_Status stat;
	//double x = DBL_MAX;

	//Parallel setup	
	MPI_Comm_rank(cart_comm, &rank);

	if (rank == 0){	
	
		//Init matrix
		global = matrix_new(first->size, first->dim);
		printf("Global matrix is: \n");
		for(i=0;i<first->size;i++){
			for(j=0;j<first->dim;j++){
				global[i][j] = rand() % 20;
				printf("%2f ", global[i][j]);
			}
			printf("\n");
		}
	}
	
	//create local matrix
	local = matrix_new(first->size/procgridsize, first->dim/procgridsize);  
	//Create datatype to describe the submatrices of the global matrix
	int sizes[2] = {gridsize,gridsize};
	int subsize[2] = {gridsize/procgridsize,gridsize/procgridsize};
	int starts[2] = {0,0};
	MPI_Datatype grid, subgrid;
	MPI_Type_create_subarray(2, sizes, subsize, starts, MPI_ORDER_C, MPI_DOUBLE, &grid);
	MPI_Type_create_resized(grid, 0, gridsize/procgridsize*sizeof(int), &subgrid);
	MPI_Type_commit(&subgrid); 
 
	double *globalptr = NULL;
	if (rank == 0){
		globalptr = &(global[0][0]);
	}

	//Scatter to all procs
	int sendcount[procgridsize*procgridsize];
	int disp[procgridsize*procgridsize];

	if (rank == 0){
		for(i=0;i<procgridsize*procgridsize; i++){
			sendcount[i] = 1;
		}
		int d = 0;
		for(i=0;i<procgridsize;i++){
			for(j=0;j<procgridsize;j++){
				disp[i*procgridsize+j] = d;
				d += 1;
			}
			d += ((gridsize/procgridsize)-1)*procgridsize;
		}
	}
	
	MPI_Scatterv(globalptr, sendcount, disp, subgrid, &(local[0][0]),
			gridsize*gridsize/(procgridsize*procgridsize),MPI_DOUBLE,0,cart_comm);

	//Procs print their local data
	int p;
	for (p=0;p<nproc;p++){     
		if(rank == p){
			printf("Local proc on rank %d is: \n", rank);
			for (i=0;i<gridsize/procgridsize;i++){
				putchar('|');
				for (j=0;j<gridsize/procgridsize;j++){
					printf("%10f ", local[i][j]);
				}
				printf("|\n");
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	
 
		//Now process local grids before gathering back to 0
		for (i=0;i<gridsize/procgridsize;i++){
			for(j=0;j<gridsize/procgridsize;j++){
				local[i][j] += 1;
			}
		}

		/*
		//Do send and recv
		if (rank == p){
	
			//Init second list
			list_b_t *second = malloc(sizeof(list_b_t));
			// allocate memory for the best position buffer
   			second->solution = (double *)malloc(first->dim * sizeof(double));
			second->error = x;
			MPI_Bcast(&second->error,1,MPI_DOUBLE,0,cart_comm);	
	
			//Move data to buffer	
			for(i=0;i<first->size;i++){
	     			memmove((void *)second->solution, (void *)global[i], sizeof(double) * first->dim);	
				//Send
				MPI_Send(second->solution,first->dim,MPI_DOUBLE,1,tag,cart_comm);
			}	
			free(second->solution);
			free(second);
		}
		else {
	
			list_b_t *second = malloc(sizeof(list_b_t));
			second->error = x;
			MPI_Bcast(&second->error,1,MPI_DOUBLE,0,cart_comm);		
	
			second->solution = (double *)malloc(first->dim * sizeof(double));
			MPI_Recv(second->solution,first->dim,MPI_DOUBLE,0,tag,cart_comm, &stat);
			//MPI_Recv(&second->error,first->dim,MPI_DOUBLE,src,tag,MPI_COMM_WORLD, &stat);
		
			for (i=0;i<first->dim;++i){
				printf("Received solution %f from rank %d\n",second->solution[i],rank);
			}	
		}
		*/
	}
	MPI_Gatherv(&(local[0][0]), gridsize*gridsize/(procgridsize*procgridsize), MPI_DOUBLE, globalptr, sendcount, disp, subgrid, 0, cart_comm);
	
	//Free type resource
	MPI_Type_free(&subgrid);
	
	//Display new grid
	if(rank == 0){
		for (i=0;i<gridsize;i++){
			for(j=0;j<gridsize;j++){
				printf("%2f ", global[i][j]);
			}
		printf("\n");
		}		
		matrix_free(local, first->size/procgridsize);
		matrix_free(global, first->size);
	}
}

