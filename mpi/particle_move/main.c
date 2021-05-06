#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <getopt.h>

#include "funcs.h"
#include "mpi.h"

#define N 10

int main(int argc, char **argv){

	int i, o;
	int nproc, rank;
	double **g1, **g2;
	int ndims[2];
	int pbc[2] = {0,0};
	//int ys, ye, xs, xe;
	int up, down, left, right;
	int coords[2];
	int nrows, ncols, chunk_rows, chunk_cols;
	MPI_Datatype coltype, rowtype;
	MPI_Comm cart_comm; 
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

	//Create dims
	MPI_Dims_create(nproc, 2, ndims);

	//Allocate buffer space for second receiver list
	list_a_t *first = NULL;
	list_b_t second;

	first = list_new(30);	
	
	while((o = getopt(argc, argv, "r:c:")) != -1){
		switch(o){
			case 'r':
				nrows = atoi(optarg);
				break;
			case 'c':
				ncols = atoi(optarg);	
				break;
		}
	}


	//Set grid dimensions
	if(nrows != ncols){
		printf("ERROR: Not possible with given nrows/ncols\n");
		MPI_Finalize();
		return -1;
	} else {
		chunk_rows = nrows;
		chunk_cols = ncols;
	}
	
	MPI_Bcast(&chunk_rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&chunk_cols, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	g1 = malloc((chunk_rows+2) * sizeof(double * ));
	g1[0] = malloc((chunk_rows+2) * (chunk_cols+2) * sizeof(double));
	for(i=0;i<(chunk_rows+2);i++){
		g1[i] = g1[0] + 1 * (chunk_cols+2);
	}

	g2 = malloc((chunk_rows+2) * sizeof(double * ));
	g2[0] = malloc((chunk_rows+2) * (chunk_cols+2) * sizeof(double));
	for(i=0;i<(chunk_rows+2);i++){
		g2[i] = g2[0] + 1 * (chunk_cols+2);
	}	
		
	//Create Cartesian communicator
	MPI_Cart_create(MPI_COMM_WORLD, 2, ndims, pbc, 0, &cart_comm);
	MPI_Cart_coords(cart_comm, rank, 2, coords);
	//get neighbours
	MPI_Cart_shift(cart_comm,0,1,&up,&down);
	MPI_Cart_shift(cart_comm,1,1,&left,&right);
	//vector datatype
	MPI_Type_vector(1,chunk_cols,chunk_cols+2,MPI_DOUBLE,&rowtype);
	MPI_Type_vector(chunk_rows,1,chunk_cols+2,MPI_DOUBLE,&coltype);
	MPI_Type_commit(&rowtype);
	MPI_Type_commit(&coltype);


	printf("%d, RpC: %d, CpC: %d\n", rank, nrows, ncols);

	
	//Space for solution buffer
	second.solution = (double *)malloc(first->dim * sizeof(double));

	//RUN ALGS
	printf("This is the result of a standard memmove: \n");
	list(first, &second);	

	printf("This is the result of the MPI version: \n");
	list_mpi(first, &second);	

	free(second.solution);
	
	//Free first struct	
	free(first);

	//Free datatypes and associated resources
	MPI_Type_free(&coltype);
	MPI_Type_free(&rowtype);
	MPI_Comm_free(&cart_comm);
	MPI_Finalize();

	free(g1[0]);
	free(g1);
	free(g2[0]);
	free(g2);

	return 0;
}
