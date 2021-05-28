#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>

#include "mpi.h"
#include "funcs.h"

#define NP 4
#define GRID_SIZE 2

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



void list_mpi(list_a_t *first, list_b_t *second, int row, int col, MPI_Comm cart_comm, MPI_Comm row_comm, MPI_Comm col_comm){
	
	int iter;
	srand(time(NULL));
	double **global, *local, *local_array, *solution_array;
    	int rank;
	MPI_Status stat;
	double x = DBL_MAX;
	int matsize[2], matsize_local;
	int nrows, ncols;
	int irow, icol;
	int pproc = 0;
	int iproc, jproc, myid;
	int nrows_local, ncols_local;
	int global_row_id, global_col_id, local_id;	
	int src, dest, recv_tag, send_tag;

	//Parallel setup	
	MPI_Comm_rank(cart_comm, &rank);

	if(rank == 0){
		//Set matrix size
		matsize[0] = nrows = first->size;
		matsize[1] = ncols = first->dim;

		//Init global matrix
		global = (double **)malloc (nrows * sizeof(double *));
		for (irow=0; irow<nrows; irow++){
			global[irow] = (double *)malloc(ncols * sizeof(double));
			for(icol=0; icol<ncols; icol++){
				global[irow][icol] = rand() % 20;
				printf("%2f ", global[irow][icol]);
			}
			printf("\n");
		}
	}
	MPI_Barrier(cart_comm);
	
	//Bcast matrix size to all procs
	MPI_Bcast(matsize, 2, MPI_INT, 0, cart_comm);
	nrows = matsize[0];
	ncols = matsize[1];

	//Stopping condition for unequal proportions
	if (nrows != ncols){
		MPI_Finalize();
		if (rank == 0){
			printf("Matrix not proportional");
		}
		exit(-1);
	}
	
	nrows_local = nrows /pproc;
	ncols_local = ncols /pproc;

	matsize_local = nrows_local * ncols_local;
	
	//create local matrices
	local = (double *)malloc(matsize_local * sizeof(double));
	//Put data in 1 dim arrays prior to scatter
	local_array = (double *)malloc(sizeof(double) * nrows * ncols);
	//Rarrange them in the appropriate order
	if (rank == 0){
		for (iproc=0; iproc<pproc; iproc++){
			for (jproc=0; jproc<pproc; jproc++){
				myid = iproc * pproc + jproc;
					global_row_id = iproc * nrows_local * irow;
					for (icol=0; icol < ncols_local; icol++){
						local_id = (myid * matsize_local) +
							(irow * ncols_local) + icol;
						global_col_id = jproc * ncols_local + icol;
						local_array[local_id] = global[global_row_id][global_col_id];
					}
			}
		}
	}
	MPI_Barrier(cart_comm);
	//Scatter data to procs
	MPI_Scatter(local_array, matsize_local, MPI_DOUBLE, local,
		matsize_local,MPI_DOUBLE,0,cart_comm);

	//Arrange the cols and rows of the matrices using sendrecv_replace so that
	//the data that is being sent and replaced in the same buffer
	if (row != 0){
		src = (col + row) % pproc;
		dest = (col + pproc - row) % pproc;
		recv_tag = 0;
		send_tag = 0;
		MPI_Sendrecv_replace(local, matsize_local, MPI_DOUBLE, dest, send_tag, src, recv_tag, row_comm, &stat);
	}	
	if (col != 0){
		src = (row + col) % pproc;
		dest = (row + pproc - col) % pproc;
		recv_tag = 0;
		send_tag = 0;
		MPI_Sendrecv_replace(local, matsize_local, MPI_DOUBLE, dest, send_tag, src, recv_tag, col_comm, &stat);
	}
	
	//Main loop
	if (rank == 0){
		//Allocate space for solution buffer
		//Init second list
		//list_b_t second = malloc(sizeof(list_b_t));
		// allocate memory for the best position buffer
   		second->solution = (double *)malloc(nrows_local * sizeof(double));
		//Memory for output global matrix in the form of array
		solution_array = (double *)malloc(sizeof(double) * nrows);
		second->error = x;
		MPI_Bcast(&second->error,1,MPI_DOUBLE,0,cart_comm);	
	
		//Main loop for moving data
		for(iter=0; iter<pproc; iter++){
			for (irow=0; irow<nrows_local;irow++){
					//Move data to buffer	
	     				memmove(&solution_array[irow], &local_array[irow], sizeof(double) * nrows);
				MPI_Send(&solution_array,nrows,MPI_DOUBLE,second->solution[irow],0,cart_comm);
			}
		}		
		//Free resources
		free(second->solution);
		//free(second);
	} else {
		//Allocate space for solution buffer
		//Init second list
		//list_b_t second = malloc(sizeof(list_b_t));
		// allocate memory for the best position buffer
   		second->solution = (double *)malloc(nrows_local * sizeof(double));
		//Memory for output global matrix in the form of array
		solution_array = (double *)malloc(sizeof(double) * nrows);
		second->error = x;
		MPI_Bcast(&second->error,1,MPI_DOUBLE,0,cart_comm);	
	
		//Main loop for recv
		for(iter=0; iter<pproc; iter++){
			MPI_Recv(&second->solution,nrows,MPI_DOUBLE,solution_array[irow],0,cart_comm, &stat);
		}		
	}

	MPI_Barrier(cart_comm);

	//Gather output blocks at process 0
	MPI_Gather(local_array, nrows, MPI_DOUBLE, solution_array, nrows, MPI_DOUBLE, 0, cart_comm);

	//Memory for output global array for solution after rearrangement
	if (rank == 0){
		second->solution = (double *)malloc(nrows * sizeof(double));
	} 

	//Arrange the output in order
	if (rank == 0){
		for (iproc=0; iproc<pproc; iproc++){
			for (jproc=0; jproc<pproc; jproc++){
				myid = iproc * pproc + jproc;
				for (irow=0; irow<nrows_local; irow++){
					global_row_id = iproc * nrows_local * irow;
					for (icol=0; icol < ncols_local; icol++){
						second->solution[global_row_id] = solution_array[local_id];
					}
				}
			}
		}
		//Print results
		printf("RESULTS: \n");
		printf(" Process %d, Global matrix : dimension %d * %d : \n", rank, nrows, ncols);
		for(irow = 0; irow < nrows; irow++){
			for(icol=0; icol<ncols; icol++){
				printf("%5f ", global[irow][icol]);
				printf("\n");
			}
			printf("\n");
		}
		printf(" Process %d, Solution array : dimension %d * %d : \n", rank, nrows, ncols);
		for(irow = 0; irow < nrows; irow++){
				printf("%5f ", second->solution[irow]);
			}
			printf("\n");
		}

}


