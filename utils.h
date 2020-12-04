

#ifndef _UTILS_H
#define _UTILS_H

/*Initial PSO settings */
extern int popSize = 100, maxIterations = 500; 

/* Serial and parallel option */
extern bool serial, parallel, demo; //for benchmark functions

/* Path options */
extern int inRoboID;
extern double inStartX, inStartY, inEndX, inEndY;
extern double inStepSize, inVelocity;
extern double inOriginX, inOriginY, inHorizonX, inHorizonY;  // 70
//char inFileHandle[20] = "maps/sample_map_OpenRooms.txt\0";
extern char inFileHandle[] = "sample_map_OpenRooms.txt";
extern int waypoints;

/* PSO parameters */
extern double pso_c1, pso_c2, pso_w_max, pso_w_min;
extern int pso_w_strategy_select, pso_nhood_size, pso_nhood_topology_select;
extern int pso_w_strategy, pso_nhood_topology;

/* Option parsing */
extern int verbose;
extern char *inFileHandlePtr;

// PSO SETTINGS
typedef struct{

	int dim; // problem dimensionality
	double x_lo; // lower range limit
	double x_hi; // higher range limit
	
	/*
	 * double *range_lo; //lower range limit (array of length DIM)
	 * double *range_hi; //higher range limit (array of length DIM)
	 */
	
	double goal; // optimization goal (error threshold)

	double **limits; // lower and higher ranges for each X value.

	int size; // swarm size (number of particles)
	int print_every; // ... N steps (set to 0 for no output)
	int steps; // maximum number of iterations
	int step; // current PSO step
	double c1; // cognitive coefficient
	double c2; // social coefficient
	double w_max; // max inertia weight value
	double w_min; // min inertia weight value

	int numset; // Set of numbers to use as X values. Default = DECIMAL. 

	int clamp_pos; // whether to keep particle position within defined bounds (TRUE)
	// or apply periodic boundary conditions (FALSE)
	int nhood_strategy; // neighborhood strategy (see PSO_NHOOD_*)
	int nhood_size; // neighborhood size
	int w_strategy; // inertia weight strategy (see PSO_W_*)

	void *rng; // pointer to random number generator (use NULL to create a new RNG)
	long seed; // seed for the generator

}pso_settings_t;

//TIMER
typedef struct{
	float time_pso; //Timing for serial
	float time_omp; // Timing for OpenMP optimized version
	float time_mpipso; //Timing for MPI-only version
	float time_mpiomp; //Timing for MPI-OpenMP optimized version
}pso_timer;

int parse_arguments(int argc, char **argv);
void options();
void pso_set_path_settings(pso_settings_t *settings, pso_params_t *params, env_t *env, robot_t *robot, int waypoints);
void print_usage();

#endif
