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
#include "path.h"
#include "defs.h"

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


	
	//Set up MPI environment
	int nproc;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

	//Parse arguments and print options
	parse_arguments(argc,argv);	
			
	//Initialise PSO settings - comment out when executing path application
	//pso_settings_t *settings = NULL; 
		

	//Initialise timer
	//struct timing_report* stats = malloc(sizeof(double));
			
	//Begin timer
	//start_timer(&(stats->parallel_time));
	
	//Execute
	//pso_parallel(settings,argc,argv);
	pso_parallel_path(argc,argv);

	//Stop timer
	//end_timer(&(stats->parallel_time));

	MPI_Barrier(MPI_COMM_WORLD);

	//Print timing
	//print_elapsed_time((char*) "PARALLEL ", stats->parallel_time.start, stats->parallel_time.finish);

	//Free timer
	//free(stats);
	
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
    		settings->size = 40;
    		settings->nhood_size = 5;
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


void pso_parallel_path(int argc, char **argv) {
   
		//Initial PSO settings
	        parse_arguments(argc, argv);

		//Get weighting and topology	
    		pso_w_strategy = getPSOParam_w_strategy(pso_w_strategy_select);
    		pso_nhood_topology = getPSOParam_nhood_topology(pso_nhood_topology_select);

		//Print info
    		printf ("Dimension = (%f,%f), Start = (%f,%f), Target = (%f,%f)\n", 
            		inHorizonX, inHorizonY, inStartX, inStartY, inEndX, inEndY);
    		printf ("Map File = sample map\n");
    		printf ("PSO: c1 = %f, c2 = %f, weight strategy = %d, neighborhood topology = %d\n", 
           		pso_c1, pso_c2, pso_w_strategy, pso_nhood_topology);
    		if (pso_w_strategy == PSO_W_LIN_DEC)
        		printf ("\tweight min = %f, weight max = %f\n", pso_w_min, pso_w_max);
    		if (pso_nhood_topology_select == PSO_NHOOD_RANDOM)
        		printf("\tneighborhood size = %d\n", pso_nhood_size);

    		//Read occupancy map
    		int **map = readMap(inFileHandlePtr, inHorizonY, inHorizonX);
		
    		//Initialize uav
    		uav_t * uav = initUav(inUavID, inStartX, inStartY, inEndX, inEndY, inStepSize, inVelocity);
    		printUav(uav);

    		//Init pso objecttive function params 
    		pso_params_t *pso_params;
    		pso_params = malloc (sizeof (pso_params_t) * 1);
    		pso_params->env = initEnv(inOriginX, inOriginY, inHorizonX, inHorizonY, map);
    		pso_params->start[0] = uav->position_coords[0];
    		pso_params->start[1] = uav->position_coords[1];
    		pso_params->stop[0] = uav->target_coords[0];
    		pso_params->stop[1] = uav->target_coords[1];
    		pso_params->c1 = pso_c1;
    		pso_params->c2 = pso_c2;
    		pso_params->w_strategy = pso_w_strategy;
    		pso_params->w_min = pso_w_min;
    		pso_params->w_max = pso_w_max;
    		pso_params->nhood_topology = pso_nhood_topology;
    		pso_params->nhood_size = pso_nhood_size;
    		printEnv(pso_params->env);
    			
		//Initialise settings
		pso_settings_t settings;
	        pso_serial_settings(&settings);	
		// Define objective function and path settings
    		pso_obj_fun_t obj_fun = pso_path;   
		pso_set_path_settings(&settings, pso_params, pso_params->env, uav, waypoints);
		
		// Set the problem specific settings
		//settings.size = popSize;
		//settings.nhood_strategy = PSO_NHOOD_RING;
    		settings.dim = waypoints * 100;
		printf("\t Number of waypoints = %d\n", settings.dim);
		//settings.nhood_size = 10;
		//settings.w_strategy = PSO_W_LIN_DEC;
    		settings.steps = 100001;
    		settings.print_every = 1000;
    
		// Init global best solution
    		pso_result_t solution;
    
 		// allocate memory for the best position buffer
    		solution.gbest = malloc(settings.dim * sizeof(double));

    		// Run pso algorithm
    		pso_solve(obj_fun, pso_params, &solution, &settings);
    	
		// Display best result - WILL BE GROUPED IN PRINT FUNCITON
    		int i, count = 0;
    		printf ("Solution waypoints:\n");
    		for(i=0;i<settings.dim/2;i++){
        		printf ("(%f, %f)\n", solution.gbest[count], solution.gbest[count + 1]);
        		count = count + 2;
    		}
    		printf("Solution distance: %f\n", solution.error);


    		int obstacles = pso_path_countObstructions(solution.gbest, settings.dim, pso_params);
    		printf ("obstacles: %d\n", obstacles);

		// Free global best buffer
    		free(solution.gbest);

 		
}


 
