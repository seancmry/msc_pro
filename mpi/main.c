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

	double *old_val, *new_val, *f_val;
	double **old, **new, **f;
	int myid, nprocs;
	int dims = 2, ndims[2];
	int periods[2] = {1,0};
	int ys, ye, xs, xe;
	int nbrup, nbrdown, nbrleft, nbrright;
	int coords[2];
	int nx, ny;
    	MPI_Comm cart_comm;
	MPI_Datatype coltype;
	
	//Parse arguments and print options
	parse_arguments(argc,argv);	
	
	//Initialise PSO settings
	pso_settings_t *settings = NULL; 
			
	//Set up MPI environment
	int rank, size;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	//Set grid dimensions
	MPI_Bcast(&nx, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	ny = nx;

	//Allocate arrays to store grids
	old_val = (double*)calloc((nx+2)*(ny+2),sizeof(double));
	new_val = (double*)calloc((nx+2)*(ny+2),sizeof(double));
	f_val = (double*)calloc((nx+2)*(ny+2),sizeof(double));
	old = (double**)malloc((nx+2)*sizeof(double*));
	new = (double**)malloc((nx+2)*sizeof(double*));
	f = (double**)malloc((nx+2)*sizeof(double*));

	//init_arr(nx+2, ny+2, old_val, uold);
	//init_arr(nx+2, ny+2, new_val, unew);
	//init_arr(nx+2, ny+2, f_val, f);

	//Get dimensions of the Cartesian grid and set up the communicator
	calculate_dims(nprocs, ndims, settings);
	MPI_Cart_create(MPI_COMM_WORLD, dims, ndims, periods, 0, &cart_comm);

	MPI_Cart_coords(cart_comm, myid, 2, coords);
	decomp2d(nx, ny, ndims[0], ndims[1], coords, &xs, &xe, &ys, &ye); 

	//Test for periodicity
	if (myid == 0){
		printf("===========================================================\n");
		printf("Initial grid with periodic boundary conditions: \n\n");
	}

	//init_range(unew, uold, f, xs, xe, ys, ye, nx, ny, &fone, &fone, &fone, &fone);
	//print_grid(unew, nx, ny, xs, xe, ys, ye, coords, settings->dim, cart_comm);

	//Get neighbours through the shift
	MPI_Cart_shift(cart_comm, 0, 1, &nbrleft, &nbrright);
	MPI_Cart_shift(cart_comm, 1, 1, &nbrdown, &nbrup);
		
	//Committing vector datatype to send columns
	MPI_Type_vector((xe-xs+1), 1, ny+2, MPI_DOUBLE, &coltype);
	MPI_Type_commit(&coltype);
		
	//Initialise timer
	struct timing_report* stats = malloc(sizeof(double));
			
	//Begin timer
	start_timer(&(stats->parallel_time));
		
	//Execute
	pso_parallel(settings,argc,argv);

	//Stop timer
	end_timer(&(stats->parallel_time));

	MPI_Barrier(MPI_COMM_WORLD);

	//Print timing
	print_elapsed_time((char*) "PARALLEL ", stats->parallel_time.start, stats->parallel_time.finish);

	//Free timer
	free(stats);
		
	//Free values
	free(f_val);
	free(new_val);
	free(old_val);
	free(f);
	free(new);
	free(old);

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
				settings = pso_settings_new(100, -32.8, 32.8);
				printf("Optimising function: ackley (dim=%d, swarm size=%d)\n",
					settings->dim, settings->size);
			} else if (strcmp(argv[1], "rosenbrock") == 0) {
            			obj_fun = pso_rosenbrock;
            			settings = pso_settings_new(100, -2.048, 2.048);
            			printf("Optimizing function: rosenbrock (dim=%d, swarm size=%d)\n",
                   			settings->dim, settings->size);
        		} else if (strcmp(argv[1], "griewank") == 0) {
            			obj_fun = pso_griewank;
            			settings = pso_settings_new(100, -600, 600);
            			printf("Optimizing function: griewank (dim=%d, swarm size=%d)\n",
                   			settings->dim, settings->size);
        		} else if (strcmp(argv[1], "sphere") == 0) {
            			obj_fun = pso_sphere;
            			settings = pso_settings_new(100, -100, 100);
            			printf("Optimizing function: sphere (dim=%d, swarm size=%d)\n",
                   			settings->dim, settings->size);
        		} else {
            			printf("Unsupported objective function: %s", argv[1]);
            			return;
        		}
    		} else if (obj_fun == NULL || settings == NULL) {
        		obj_fun = pso_sphere;
        		settings = pso_settings_new(30, -100, 100);
        		printf("Optimizing function: sphere (dim=%d, swarm size=%d)\n",
                   		settings->dim, settings->size);
    		}

    		// set some general PSO settings
    		settings->goal = 1e-5;
    		settings->size = 20;
    		settings->nhood_strategy = PSO_NHOOD_RING;
    		settings->nhood_size = 10;
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


