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

#define DIMS 2


int main(int argc, char **argv){

	int nproc, row, col, rank;
	int ndims[2], subdims[2], coords[2];
	int c;
	bool serial = false;
	bool mpi = false;
	int periods[2] = {0, 0}; //Make both dims periodic


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
		
		printf("This is the result of the MPI version: \n");
		double t2, t1;	
		MPI_Comm cart_comm, row_comm, col_comm;
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
	
        	MPI_Cart_create(MPI_COMM_WORLD, DIMS, ndims, periods, 1, &cart_comm); 
                           
		//Get my coords in the new communicator
    		MPI_Cart_coords(cart_comm, rank, DIMS, coords);	

		//Create and set new row and column communicators
		row = coords[0];
		col = coords[1];
		subdims[0] = 0;
		subdims[1] = 1;
		
		//Outputs column communicator containing processes
		MPI_Cart_sub(cart_comm, subdims, &row_comm);
		
		subdims[0] = 1;
		subdims[1] = 0;
		
		//Outputs row communicator
		MPI_Cart_sub(cart_comm, subdims, &col_comm);


		t1 = MPI_Wtime();
		//Execute list_mpi
		list_mpi(first,&second,row,col,cart_comm,row_comm,col_comm);
		t2 = MPI_Wtime();
		printf("Time taken: %2.2f\n",t2-t1);

			
		MPI_Comm_free(&cart_comm);	
		MPI_Finalize();
	}	
	
	//Free first struct	
	free(first);

	return 0;
}
