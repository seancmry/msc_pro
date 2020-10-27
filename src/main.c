
#include <math.h>
#include <stdlib.h>
#include <stdio.h> // for printf
#include <getopt.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <string.h>
#include <sys/time.h>
#include <mpi.h>
#include <omp.h>

#include "pso.h"
#include "path.h"
#include "utils.h"

struct timeval start, finish;
long _start, _end;
double overhead;

//==============================================================
//                  BENCHMARK FUNCTIONS
//==============================================================
//Will be revised and new ones implemented. This will form part of
//a discussion section. Check sources, with similar studies, and 
//what functions they have used. Some of these may include Ackley,
//Griewank, Rastrigin, Rosenbrock, Schaffer, Schwefel, etc. The ones
//detailed here are pretty standard for the evaluation of PSO in the 
//context of path planning.

double ackley(double x[], double nDimensions) {
	double c = 2*M_PI;
	double b = 0.2;
	double a = 20;
	double sum1 = 0;
	double sum2 = 0;
	int i;
	for (i=0; i<nDimensions; i++) {
		sum1 = sum1 + gsl_pow_2(x[i]);
		sum2 = sum2 + cos(c*x[i]);
	}
	double term1 = -a * exp(-b*sqrt(sum1/nDimensions));
	double term2 = -exp(sum2/nDimensions);
	return term1 + term2 + a + M_E;
}


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

	gettimeofday(&start, &finish);
	
	/*RNG setup*/
	gsl_rng_env_setup();
	gsl_rng * r = gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(r, time(0));
  
	int inUavID = 0;
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
    	char *inFileHandlePtr = NULL;
	parse_args(argc, argv);

	/*Setting up MPI environment */
	int size, myrank;
	int particles;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	double result[particles];
	
	/*Positions for MPI implementation can be vectorised*/
	double A = (double)malloc(sizeof(double)*N+1);
	double B = (double)malloc(sizeof(double)*N+1);


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
    	uav_t *uav = init_uav(inUavID, inStartX, inStartY, inEndX, inEndY, inStepSize, inVelocity);
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
    	int maxIter = 500;
    	int pop_size = 100;

    	/* Use pso library */
    	
	// Define objective function
    	pso_obj_fun_t obj_fun = pso_path;
    	
	// Set the problem specific settings
    	pso_set_path_settings(&settings, pso_params, pso_params->env, uav, waypoints);
	settings.size = pop_size;
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

	//Free global best buffer
	free(solution.gbest);

	/* 
	 *
	 * DEMO 
	 *
	 *
	 *
    	if (argc == 2) {
		if (strcmp(argv[1], "ackley") == 0) {
			obj_fun = pso_ackley;
			settings = pso_settings_new(1000, x_lo, x_hi);
			printf("Optimizing function: ackley (dim=%d, swarm size = %d)\n",
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
*/


/*
 *
 * MPI VERSION
 *
 *
 */

	if(myrank == 0) {
		particles = pso_nhood_size/size;
		printf("%d particles were distributed\n", particles);
	}
	MPI_Bcast(particles,1,MPI_INT,0,MPI_COMM_WORLD);
	if(myrank == 0) {
		particles += pso_nhood_size%size;
		printf("%d particles were distributed \n", particles);
	}

	/*Initialise particles */

	#pragma omp parallel for private(A,B) reduction(min:gbestfitness)
	for (i=0; i<particles; i++) {
		// #pragma omp parallel for private (A,B);
		for (j=0; j<ndims; j++) {
			A = x_lo + (x_hi - x_lo) * gsl_rng_uniform(r);
	                B = x_lo + (x_hi - x_lo) * gsl_rng_uniform(r);
				positions[i][j] = a;
				pbestpositions[i][j] = a;
				velocities[i][j] = (a-b)/2;
		}
		pbestfitness[i] = ackley(positions[i], ndims);

	if (pBestFitness[i] < gBestFitness) {
            	memmove((void *)gBestPosition, (void *)&positions[i], sizeof(double) * nDimensions);
            	gBestFitness = pBestFitness[i];
        	} 
    	}

//FIXME
    	//actual calculation
    	for (step=0; step<nIterations; step++) {
        	#pragma omp parallel num_threads(4) shared(min)
        	{

            	#pragma omp for private(a,b) 
            	for (i=0; i<distributed_particles; i++) {
                 
                    	for (j=0; j<nDimensions; j++) {
                        	// calculate stochastic coefficients
                        	rho1 = c1 * gsl_rng_uniform(r);
                        	rho2 = c2 * gsl_rng_uniform(r);
                        	// update velocity
                        	velocities[i][j] = w * velocities[i][j] + \
                        	rho1 * (pBestPositions[i][j] - positions[i][j]) +  \
                        	rho2 * (gBestPosition[j] - positions[i][j]);
                        	// update position
                        	positions[i][j] += velocities[i][j];

                        	if (positions[i][j] < x_min) {
                            		positions[i][j] = x_min;
                            		velocities[i][j] = 0;
                        	} else if (positions[i][j] > x_max) {
                           		positions[i][j] = x_max;
                            		velocities[i][j] = 0;
                        	}

                    	}

                	// update particle fitness
                	fit = ackley(positions[i], nDimensions);
                	// update personal best position?
                	if (fit < pBestFitness[i]) {
                    		pBestFitness[i] = fit;
                    	// copy contents of positions[i] to pos_b[i]
                    	memmove((void *)&pBestPositions[i], (void *)&positions[i],
                        	sizeof(double) * nDimensions);
                    	}   
                	// update gbest??
               
                	}

                	#pragma omp for reduction(min:gBestFitness)
                	for(i=0;i<(int)distributed_particles;i++)	
	        	if (pBestFitness[i] < gBestFitness) {
                            	// update best fitness
                            	gBestFitness = pBestFitness[i];
                            	// copy particle pos to gbest vector
                            	}
                
                	#pragma omp  for  
                 	for(i=0;i<(int)distributed_particles;i++) {
				if (gBestFitness==pBestFitness[i])
                    			min=i;  
                   	}
            
		}   
            	memmove((void *)gBestPosition, (void *)&pBestPositions[min],sizeof(double) * nDimensions);
            	for(int k=0;k<(int)nDimensions;k++)
                	sendingdata[k]=gBestPosition[k]; 
            		//#pragma omp single
                	sendingdata[(int)nDimensions]=gBestFitness;
                
			MPI_Gather(&sendingdata,nDimensions+1, MPI_INT,&recievingdata,nDimensions+1, MPI_INT, 0, MPI_COMM_WORLD);
                	if(myrank == 0) {
                    		int min=gBestFitness;
                    		int pos=-1;
                    		for(int k=0;k<size;k++)	{ 
					//printf("%d\n",recievingdata[k*((int)nDimensions+1)+((int)nDimensions)] );
                        		if(min>=recievingdata[k*((int)nDimensions+1)+((int)nDimensions)])
                            	{
                                	min=recievingdata[k*((int)nDimensions+1)+((int)nDimensions)];
                                	pos=k*((int)nDimensions+1);
                            	}   
                    	}
                    	gBestFitness=min;
                    	int k=0;
                    
                    	for(k=pos;k<(int)nDimensions+pos;k++)
                        	gBestPosition[k-pos]=recievingdata[k];                  
                	}
                	MPI_Bcast(&gBestPosition,nDimensions,MPI_INT,0,MPI_COMM_WORLD);
                
        
    	}

    	if(myrank == 0) {

        	printf("Result: %f\n", gBestFitness);
        	gettimeofday(&TimeValue_Final, &TimeZone_Final);
        	time_start = TimeValue_Start.tv_sec * 1000000 + TimeValue_Start.tv_usec;
        	time_end = TimeValue_Final.tv_sec * 1000000 + TimeValue_Final.tv_usec;
        	time_overhead = (time_end - time_start)/1000000.0;
       		printf("\n Time in Seconds (T) : %lf\n",time_overhead);
	}   
    	gsl_rng_free(r);
    	MPI_Finalize();

//mpicc -fopenmp mpi_omp.c -lm -lgsl -lgslcblas -o mpi_omp

    	return 0;
}


