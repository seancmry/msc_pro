#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>
#include <ctype.h>

#include "mpi.h"
#include <omp.h>

#include "utils.h"
#include "pso.h"
#include "defs.h"
//#include "path.h"


//==============================================================
//                  BENCHMARK FUNCTIONS
//==============================================================


double pso_ackley(double *vec, int dim, void *params) {
	double c = 2*M_PI;
	double b = 0.2;
	double a = 20;
	double sum1 = 0;
	double sum2 = 0;
	int i;
	for (i=0; i<dim; i++) {
		sum1 = sum1 + gsl_pow_2(vec[i]);
		sum2 = sum2 + cos(c*vec[i]);
	}
	double term1 = -a * exp(-b*sqrt(sum1/dim));
	double term2 = -exp(sum2/dim);
	
	return term1 + term2 + a + M_E;
}



double pso_sphere(double *vec, int dim, void *params) {

    double sum = 0;
    int i;
    for (i=0; i<dim; i++)
        sum+=(vec[i]*vec[i]);

    return sum;
}



double pso_rosenbrock(double *vec, int dim, void *params) {

    double sum = 0;
    int i;
    for (i=0; i<dim-1; i++)
        sum += 100 * pow((vec[i+1] - pow(vec[i], 2)), 2) +	\
            pow((1 - vec[i]), 2);

    return sum;

}


double pso_griewank(double *vec, int dim, void *params) {

    double sum = 0.;
    double prod = 1.;
    int i;
    for (i=0; i<dim;i++) {
        sum += pow(vec[i], 2);
        prod *= cos(vec[i] / sqrt(i+1));
    }

    return sum / 4000 - prod + 1;

}




int main(int argc, char **argv) {

	double **g1, **g2;
	int ndims[2] = {0,0};
	int pbc[2] = {0,0};
	int ys, ye, xs, xe;
	int up, down, left, right;
	int coords[2];
	int chunk_rows, chunk_cols;
	int nrows, ncols;
	int i;
    	MPI_Comm cart_comm;
	MPI_Datatype coltype, rowtype;
		
	//Set up MPI environment
	int nproc, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	
	//Create dims
	//MPI_Dims_create(nproc, 2, ndims);

	//Parse arguments and print options
	parse_arguments(argc,argv);	
			
	//Initialise PSO settings
	pso_settings_t *settings = NULL; 
	/*		
	//Set grid dimensions
	//MPI_Bcast(&nrows, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	//nrows = ncols;

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
	

	//INITIALIZE DIRICHLET BOUNDARY CONDITIONS
	//init_range(g1,rank,chunk_rows+2,chunk_cols+2,&xe,&xs,&ye,&ys);
	//init_range(g2,rank,chunk_rows+2,chunk_cols+2,&xe,&xs,&ye,&ys);
			
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

	//Test for periodicity
	//if (myid == 0){
	//	printf("===========================================================\n");
	//	printf("Initial grid with periodic boundary conditions: \n\n");
	//}

	//init_range(unew, uold, f, xs, xe, ys, ye, nx, ny, &fone, &fone, &fone, &fone);
	//print_grid(unew, nx, ny, xs, xe, ys, ye, coords, settings->dim, cart_comm);
	*/
	
	//Initialise timer
	//struct timing_report* stats = malloc(sizeof(double));
			
	//Begin timer
	//start_timer(&(stats->parallel_time));
	
	//Execute
	pso_parallel(settings,argc,argv);

	//Stop timer
	//end_timer(&(stats->parallel_time));

	//MPI_Barrier(MPI_COMM_WORLD);

	//Print timing
	//print_elapsed_time((char*) "PARALLEL ", stats->parallel_time.start, stats->parallel_time.finish);

	//Free timer
	//free(stats);
	
	//MPI_Comm_free(&cart_comm);
	MPI_Finalize();

	return 0;    
}


void pso_parallel(pso_settings_t *settings, int argc, char **argv) {
	
		//Initialise function settings
		pso_obj_fun_t obj_fun = NULL;

    		// parse command line argument (function name) - NOTE: Booleans should dictate the use of both serial and demo functions, so argc == 3
    		if (argc == 2) {
        		if (strcmp(argv[1], "ackley") == 0) {
				obj_fun = pso_ackley;
				settings = pso_settings_new(1000, -32.8, 32.8);
				printf("Optimising function: ackley (dim=%d, swarm size=%d)\n",
					settings->dim, settings->size);
			} else if (strcmp(argv[1], "rosenbrock") == 0) {
            			obj_fun = pso_rosenbrock;
            			settings = pso_settings_new(1000, -2.048, 2.048);
            			printf("Optimizing function: rosenbrock (dim=%d, swarm size=%d)\n",
                   			settings->dim, settings->size);
        		} else if (strcmp(argv[1], "griewank") == 0) {
            			obj_fun = pso_griewank;
            			settings = pso_settings_new(1000, -600, 600);
            			printf("Optimizing function: griewank (dim=%d, swarm size=%d)\n",
                   			settings->dim, settings->size);
        		} else if (strcmp(argv[1], "sphere") == 0) {
            			obj_fun = pso_sphere;
            			settings = pso_settings_new(1000, -100, 100);
            			printf("Optimizing function: sphere (dim=%d, swarm size=%d)\n",
                   			settings->dim, settings->size);
        		} else {
            			printf("Unsupported objective function: %s", argv[1]);
            			return;
        		}
    		} else if (obj_fun == NULL || settings == NULL) {
        		obj_fun = pso_sphere;
        		settings = pso_settings_new(1000, -100, 100);
        		printf("Optimizing function: sphere (dim=%d, swarm size=%d)\n",
                   		settings->dim, settings->size);
    		}

    		// set some general PSO settings
    		settings->goal = 1e-5;
    		settings->size = 90;
    		settings->nhood_size = 10;
		settings->nhood_strategy = PSO_NHOOD_RING;
    		settings->w_strategy = PSO_W_LIN_DEC;

    		// initialize GBEST solution
    		pso_result_t solution;
    
		// allocate memory for the best position buffer
    		solution.gbest = (double *)malloc(settings->dim * sizeof(double));

    		// run optimization algorithm
    		pso_solve(obj_fun, NULL, &solution, settings);

    		// free the gbest buffer
    		free(solution.gbest);
		
		//free settings
		pso_settings_free(settings);

}


