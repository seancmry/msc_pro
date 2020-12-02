#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <unistd.h>
#include <stdbool.h>
#include <ctype.h>
#include "mpi.h"

#include "path.h"
#include "pso.h"
#include "utils.h"


//==============================================================
//                  BENCHMARK FUNCTIONS
//==============================================================

//GET ACKLEY


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




int main(int argc, char **argv) {

/*
    pso_settings_t *settings = NULL;
    pso_obj_fun_t obj_fun = NULL;

    // parse command line argument (function name)
    if (argc == 2) {
        if (strcmp(argv[1], "rosenbrock") == 0) {
            obj_fun = pso_rosenbrock;
            settings = pso_settings_new(30, -2.048, 2.048);
            printf("Optimizing function: rosenbrock (dim=%d, swarm size=%d)\n",
                   settings->dim, settings->size);
        } else if (strcmp(argv[1], "griewank") == 0) {
            obj_fun = pso_griewank;
            settings = pso_settings_new(30, -600, 600);
            printf("Optimizing function: griewank (dim=%d, swarm size=%d)\n",
                   settings->dim, settings->size);
        } else if (strcmp(argv[1], "sphere") == 0) {
            obj_fun = pso_sphere;
            settings = pso_settings_new(30, -100, 100);
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
        settings = pso_settings_new(30, -100, 100);
        printf("Optimizing function: sphere (dim=%d, swarm size=%d)\n",
                   settings->dim, settings->size);
    }

    // set some general PSO settings
    settings->goal = 1e-5;
    settings->size = 30;
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

    // free the settings
    pso_settings_free(settings);
*/



	//See mpi pso if unsure
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    //Parse arguments and print options
    parse_arguments(argc,argv);
    options();

    //Get weighting and topology	
    pso_w_strategy = getPSOParam_w_stategy(pso_w_strategy_select);
    pso_nhood_topology = getPSOParam_nhood_topology(pso_nhood_topology_select);

    /* Read occupancy map */
    int **map = readMap (inFileHandlePtr, inHorizonY, inHorizonX);

    /* Initialize robot */
    robot_t * robot = initRobot(inRoboID, inStartX, inStartY, inEndX, inEndY, inStepSize, inVelocity);
    printRobot(robot);

    /* Init pso objecttive function params */
    pso_params_t *pso_params;
    pso_params = malloc (sizeof (pso_params_t) * 1);
    pso_params->env = initEnv(inOriginX, inOriginY, inHorizonX, inHorizonY, map);
    pso_params->start[0] = robot->position_coords[0];
    pso_params->start[1] = robot->position_coords[1];
    pso_params->stop[0] = robot->target_coords[0];
    pso_params->stop[1] = robot->target_coords[1];
    pso_params->c1 = pso_c1;
    pso_params->c2 = pso_c2;
    pso_params->w_strategy = pso_w_strategy;
    pso_params->w_min = pso_w_min;
    pso_params->w_max = pso_w_max;
    pso_params->nhood_topology = pso_nhood_topology;
    pso_params->nhood_size = pso_nhood_size;
    printEnv(pso_params->env);

    /* Init pso settings */
    pso_settings_t settings;
    // Set the default settings
    pso_set_default_settings(&settings);

	if(serial) {



	}
	

	if(parallel) {






	}
	
    //PARALLEL HERE

    /* Use pso library */
    // Define objective function
    pso_obj_fun_t obj_fun = pso_path;
    // Set the problem specific settings    
    pso_set_path_settings(&settings, pso_params, pso_params->env, robot, waypoints);
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
  

    //PARALLEL HERE  
  
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

	

    // Free global best buffer
    free(solution.gbest);

    return 0;
}


int empty_file(FILE* file) {
	size_t size;

	fseek(file, 0, SEEK_END);
	size = ftell(file);

	return size ? 0 : 1;

}

void timing_to_file(char *fname) {
	FILE *fptr;
	
	fptr = fopen(fname, "a+");
	if (!fptr)
		printf("Could not open file %s\n", fname);
	
	if (empty_file(fptr))
		fprintf(fptr, "Data printed below - change this row to headings\n");
	if (parallel) {
		fprintf(fptr, "parallel and serial data here\n");
	}
	else {
		fprintf(fptr, "serial data here\n");
	}
	fclose(fptr);
}


