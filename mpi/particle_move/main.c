#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <getopt.h>
#include <stddef.h>

#include "funcs.h"
#include "mpi.h"

#define N 1
#define M 1
#define UP 0
#define DOWN 1
#define LEFT 2
#define RIGHT 3

int main(int argc, char **argv){

	int nproc, rank;
	//int ndims[2] = {N,M};
	//int pbc[2] = {0,0};
	//int coords[2], nbr[4];
	//int nrows, ncols;
	//MPI_Comm cart_comm; 
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    	if (nproc < 2) {
        	fprintf(stderr,"Requires at least two processes.\n");
        	exit(-1);
    	}

	//Allocate buffer space for second receiver list
	list_a_t *first = NULL;
	list_b_t send;

	first = list_new(30);	

	/*
	if (nproc == N){	
		nrows = ncols;
		//Set grid dimensions
		if(nrows != ncols){
			printf("ERROR: Not possible with given nrows/ncols\n");
			MPI_Finalize();
			return -1;
		}
		
		//Create dims
		MPI_Dims_create(nproc, 2, ndims);

		//Create Cartesian communicator
		MPI_Cart_create(MPI_COMM_WORLD, 2, ndims, pbc, 0, &cart_comm);
		MPI_Cart_coords(cart_comm, rank, 2, coords);
		//get neighbours
		MPI_Cart_shift(cart_comm,0,1,&nbr[UP],&nbr[DOWN]);
		MPI_Cart_shift(cart_comm,1,1,&nbr[LEFT],&nbr[RIGHT]);

		printf("%d, RpC: %d, CpC: %d\n", rank, nrows, ncols);
	}
	*/

	//Space for solution buffer
	send.solution = (double *)malloc(first->dim * sizeof(double));

	//RUN ALGS
	//printf("This is the result of a standard memmove: \n");
	//list(first, &second);	
	
	printf("This is the result of the MPI version: \n");
	list_mpi(first,send);
		
	free(send.solution);
	
	//Free first struct	
	free(first);

	//Free resources
	//MPI_Comm_free(&cart_comm);
	MPI_Finalize();

	return 0;
}
