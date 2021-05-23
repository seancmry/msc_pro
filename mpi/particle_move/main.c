#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <getopt.h>
#include <stddef.h>
#include <stdbool.h>

#include "funcs.h"
#include "mpi.h"


int main(int argc, char **argv){

	int nproc, rank;
	int ndims[2], coords[2];
	int c;
	bool serial = false;
	bool mpi = false;
	int periods[2] = {0, 0}; //Make both dims periodic
	int reorder = true; //Let MPI assign arbitrary ranks if it deems it necessary	
	int nbrup, nbrdown, nbrleft, nbrright;

	//Handle arguments
	while((c = getopt(argc,argv, "a:b")) != -1){
		switch(c){
			case 'a':
				serial = true;
				break;
			case 'b':
				mpi = true;
				break;
			default:
				fprintf(stderr, "Invalid option given\n");
				return -1;
			}
	}
	//Allocate buffer space for second receiver list
	list_a_t *first = NULL;
	list_b_t second;

	first = list_new(10);	

	//RUN ALGS
	if (serial == true){
		// allocate memory for the best position buffer
    		second.solution = (double *)malloc(first->dim * sizeof(double));
	
		//Run alg
		printf("This is the result of a standard memmove: \n");
		list(first, &second);
		
		//Free buffer
		free(second.solution);	
	}
	

	if (mpi == true){
		
		MPI_Comm cart_comm;
		MPI_Init(&argc, &argv);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	
    		if (nproc != 4 && nproc != 8) {
        		fprintf(stderr,"Requires at least four processes.\n");
        		exit(-1);
    		} 

		// Calculate dims and create a communicator given the 2D torus topology.
		calc_dims(nproc,ndims);
		
		if (rank == 0){
			printf("Process decomposition is %3d %3d\n", ndims[0], ndims[1]);
		}

		//Set periods for wraparound connection
		periods[0] = periods[1] = 1;
	
        	MPI_Cart_create(MPI_COMM_WORLD, 2, ndims, periods, reorder, &cart_comm); 
                          
		//Find new ranks
		MPI_Comm_rank(cart_comm, &rank);
    
		//Get my coords in the new communicator
    		MPI_Cart_coords(cart_comm, rank, 2, coords);	

		//Get ranks in neighbouring procs
		MPI_Cart_shift(cart_comm, 0, 1, &nbrleft, &nbrright);
		MPI_Cart_shift(cart_comm, 1, 1, &nbrdown, &nbrup);
	
		//Determine dims of local matrix block
		//
	
		// Print my location in the 2D torus.
   		printf("[MPI process %d] I am located at (%d, %d).\n", rank, coords[0],coords[1]);

		//Execute list_mpi
		printf("This is the result of the MPI version: \n");
		list_mpi(first,&second,cart_comm);
		
		MPI_Comm_free(&cart_comm);	
		MPI_Finalize();
	}	
	
	//Free first struct	
	free(first);

	return 0;
}
