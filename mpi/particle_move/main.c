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
	int c;
	bool serial = false;
	bool mpi = false;
	int ndims[2] = {0, 0};
	int periods[2] = {0, 0}; //Make both dims periodic
	int reorder = true; //Let MPI assign arbitrary ranks if it deems it necessary
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	
	MPI_Dims_create(nproc, 2, ndims);

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
 
        // Create a communicator given the 2D torus topology.
        MPI_Comm cart_comm;
        MPI_Cart_create(MPI_COMM_WORLD, 2, ndims, periods, reorder, &cart_comm);
                           
	//Declare neighbours
	enum DIRECTIONS {DOWN,UP,LEFT,RIGHT};
	char* neighbours_names[4] = {"down", "up", "left", "right"};
	int neighbours_ranks[4];
    
	//Let consider dims[0] = X, so the shift tells us our left and right neighbours
	MPI_Cart_shift(cart_comm, 0, 1, &neighbours_ranks[LEFT], &neighbours_ranks[RIGHT]);
	      
	// Let consider dims[1] = Y, so the shift tells us our up and down neighbours
	MPI_Cart_shift(cart_comm, 1, 1, &neighbours_ranks[DOWN], &neighbours_ranks[UP]);

        // My rank in the new communicator
        MPI_Comm_rank(cart_comm, &rank);
	
	//Get my coords in the new communicator
	int my_coords[2];
    	MPI_Cart_coords(cart_comm, rank, 2, my_coords);
 
	for(int i = 0; i < 4; i++)
    	{
        if(neighbours_ranks[i] == MPI_PROC_NULL){
            printf("[MPI process %d] I have no %s neighbour.\n", rank, neighbours_names[i]);
        }else{
            printf("[MPI process %d] I have a %s neighbour: process %d.\n", rank, neighbours_names[i], neighbours_ranks[i]);
		}
    	}

    	// Print my location in the 2D torus.
   	printf("[MPI process %d] I am located at (%d, %d).\n", rank, my_coords[0],my_coords[1]);
        
	
	//Allocate buffer space for second receiver list
	list_a_t *first = NULL;
	list_b_t second;

	first = list_new(30);	

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
    		if (nproc < 2) {
        		fprintf(stderr,"Requires at least two processes.\n");
        		exit(-1);
    		}

		printf("This is the result of the MPI version: \n");
		list_mpi(first,&second,cart_comm);
	}
		
	//Free first struct	
	free(first);

	MPI_Finalize();

	return 0;
}
