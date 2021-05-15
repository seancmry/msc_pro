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

	int nproc;
	int c;
	bool serial = false;
	bool mpi = false;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	
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
	
    	if (nproc < 2) {
        	fprintf(stderr,"Requires at least two processes.\n");
        	exit(-1);
    	}

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
		printf("This is the result of the MPI version: \n");
		list_mpi(first,&second);
	}
		
	//Free first struct	
	free(first);

	MPI_Finalize();

	return 0;
}
