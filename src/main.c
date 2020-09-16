
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h> // for printf
#include <getopt.h>
#include <sys/time.h>

#include "pso.h"
#include "path.h"
#include "utils.h"

//==============================================================
//                  BENCHMARK FUNCTIONS
//==============================================================
//Will be revised and new ones implemented. This will form part of
//a discussion section. Check sources, with similar studies, and 
//what functions they have used. Some of these may include Ackley,
//Griewank, Rastrigin, Rosenbrock, Schaffer, Schwefel, etc. The ones
//detailed here are pretty standard for the evaluation of PSO in the 
//context of path planning.
double pso_sphere(double *vec, int dim, void *params) {

    	double sum = 0;
    	int i;
    	for (i=0; i<dim; i++)
        	sum += pow(vec[i], 2);
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


int main (int argc, char **argv){

    	/* Start nasty hard coded segment */
	//int inRoboID = 0;
    	//double inStartX = 70.0;
    	//double inStartY = 70.0;
    	//double inEndX = 136.0;
    	//double inEndY = 127.0;
    	//double inStepSize = 1;
    	//double inVelocity = 2;
    	//double inOriginX = 0;
    	//double inOriginY = 0;
    	//double inHorizonX = 200;
    	//double inHorizonY = 200;  // 70
    	//char inFileHandle[20] = "maps/sampleMap4.dat\0";
    	//char inFileHandle[] = "Berlin52.txt";
    	//int waypoints = 5;
    	/* End nasty hard coded segment */

    	/* PSO parameters */
    	//double pso_c1 = -1.0;
    	//double pso_c2 = -1.0;
    	//double pso_w_max = -1.0;
    	//double pso_w_min = -1.0;
    	//int pso_w_strategy_select = -1;
    	//int pso_nhood_size = -1;
    	//int pso_nhood_topology_select = -1;

    	//int pso_w_strategy = -1;
    	//int pso_nhood_topology = -1;

    	/* Option parsing */
    	//int verbose = 0;
    	//char *inFileHandlePtr = NULL;
	parse_args(argc, argv);

    	// PSO options from user selection
    	pso_w_strategy = getPSOParam_w_stategy(pso_w_strategy_select);
    	pso_nhood_topology = getPSOParam_nhood_topology(pso_nhood_topology_select);

    	// Print argument options
    	printf ("Dimension = (%f,%f), Start = (%f,%f), Target = (%f,%f)\n", 
            	inHorizonX, inHorizonY, inStartX, inStartY, inEndX, inEndY);
    	printf ("Map File = %s\n", inFileHandlePtr);
    	printf ("PSO: c1 = %f, c2 = %f, weight strategy = %d, neighborhood topology = %d\n", 
            	pso_c1, pso_c2, pso_w_strategy, pso_nhood_topology);
    	if (pso_w_strategy == PSO_W_LIN_DEC)
        	printf ("\tweight min = %f, weight max = %f\n", pso_w_min, pso_w_max);
    	if (pso_nhood_topology_select == PSO_NHOOD_RANDOM)
        	printf("\tneighborhood size = %d\n", pso_nhood_size);

    	/* Read occupancy map */
    	int **map = read_map (inFileHandlePtr, inHorizonY, inHorizonX);

    	/* Initialize robot */
    	uav_t *uav = init_uav(inRoboID, inStartX, inStartY, inEndX, inEndY, inStepSize, inVelocity);
    	print_uav(uav);

    	/* Init pso objecttive function params */
    	pso_params_t * pso_params;
    	pso_params = malloc (sizeof (pso_params_t) * 1);
    	pso_params->env = init_env(inOriginX, inOriginY, inHorizonX, inHorizonY, map);
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
    	print_env(pso_params->env);

    	/* Init pso settings */
    	pso_settings_t settings;
    	// Set the default settings
   	pso_set_default_settings(&settings);
    
    	/* PSO settings */
    	//int maxIterations = 500;
    	//int pop_size = 100;

    	/* Use pso library */
    	
	// Define objective function
    	pso_obj_fun_t obj_fun = pso_path;
    	
	// Set the problem specific settings
    	pso_set_path_settings(&settings, pso_params, pso_params->env, uav, waypoints);
	//settings.size = popSize;
	//settings.nhood_strategy = PSO_NHOOD_RING;
    	settings.dim = waypoints * 2;
	//settings.nhood_size = 10;
	//settings.w_strategy = PSO_W_LIN_DEC;
    	settings.steps = 100000;
    	settings.print_every = 10;
    	
	// Init global best solution
    	pso_result_t solution;
    	
	// Allocate mem for best position buffer
   	solution.gbest = malloc (settings.dim * sizeof(double));
    	
	// Run pso algorithm
    	pso_solve(obj_fun, pso_params, &solution, &settings);
    
    	// Display best result
    	int i, count = 0;
    	printf ("Solution waypoints:\n");
    	for(i=0;i<settings.dim/2;i++){
        	printf ("(%f, %f)\n", solution.gbest[count], solution.gbest[count + 1]);
        	count = count + 2;
    	}
    	printf("Solution distance: %f\n", solution.error);


    	int obstacles = pso_path_countObstructions(solution.gbest, settings.dim, pso_params);
    	printf ("obstacles: %d\n", obstacles);

	
	/* 
	 *
	 * DEMO 
	 *
	 *
	 *
    	if (argc == 2) {
        	if (strcmp(argv[1], "rosenbrock") == 0) {
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
            		return 1;
        	}
    	} else if (argc > 2) {
        	printf("Usage: demo [PROBLEM], where problem is optional with values [sphere|rosenbrock|griewank]\n ");
        	return 1;
    	}
 
    	// handle the default case (no argument given)
    	if (obj_fun == NULL || settings == NULL) {
        	obj_fun = pso_sphere;
        	settings = pso_settings_new(1000, -100, 100);
        	printf("Optimizing function: sphere (dim=%d, swarm size=%d)\n",
                   	settings->dim, settings->size);
    	}

  	*/

    	// initialize GBEST solution - global best
    	pso_result_t solution;
    	// allocate memory for the best position buffer - this is key in the parallel version.
    	solution.gbest = (double *)malloc(settings->dim * sizeof(double));

    	// run optimization algorithm
    	pso_solve(obj_fun, NULL, &solution, &settings);	

    	// free the settings
    	pso_settings_free(&settings);

    	// Free global best buffer
    	free(solution.gbest);

    	return 0;
}


