#include <stdio.h>
#include "utils.h"

/*Initial PSO settings */
int popSize = 100;
int maxIterations = 500; 

/* Serial and parallel option */
int serial = 0;
int parallel = 0;

/* Start nasty hard coded segment */
int inRoboID = 0;
double inStartX = 70.0;
double inStartY = 70.0;
double inEndX = 136.0;
double inEndY = 127.0;
double inStepSize = 1;
double inVelocity = 2;
double inOriginX = 0;
double inOriginY = 0;
double inHorizonX = 200;
double inHorizonY = 200;  // 70
//char inFileHandle[20] = "maps/sampleMap4.dat\0";
char inFileHandle[] = "sample_map_OpenRooms.txt";
int waypoints = 5;
/* End nasty hard coded segment */

/* PSO parameters */
double pso_c1 = -1.0;
double pso_c2 = -1.0;
double pso_w_max = -1.0;
double pso_w_min = -1.0;
int pso_w_strategy_select = -1;
int pso_nhood_size = -1;
int pso_nhood_topology_select = -1;

int pso_w_strategy = -1;
int pso_nhood_topology = -1;

/* Option parsing */
int verbose = 0;
char *inFileHandlePtr = NULL;


int parse_arguments(int argc, char **argv) {
    int c;
    opterr = 0;

    while ((c = getopt (argc, argv, "a:b:c:d:e:f:g:n:m:p:q:r:s:t:w:x:v:z")) != -1)
        switch (c) {
            case 'v':
                verbose = 1;
                break;
            case 'z':
		serial = 1;
		break;
	    case 'g':
		parallel = 1;
		break;
            case 'a':
                sscanf(optarg, "%lf", &inHorizonX);
                break;
            case 'b':
                sscanf(optarg, "%lf", &inHorizonY);
                break;
            case 'c':
                sscanf(optarg, "%lf", &inStartX);
                break;
            case 'd':
                sscanf(optarg, "%lf", &inStartY);
                break;
            case 'e': 
                sscanf(optarg, "%lf", &inEndX);
                break;
            case 'f':
                sscanf(optarg, "%lf", &inEndY);
                break;
            case 'n':
                sscanf(optarg, "%d", &waypoints);
                break;
            case 'm':
                inFileHandlePtr = optarg;
                break;
            case 'p': /* PSO c1 */
                sscanf(optarg, "%lf", &pso_c1);
                break;
            case 'q': /* PSO c2 */
                sscanf(optarg, "%lf", &pso_c2);
                break;
            case 'r': /* PSO w_max */
                sscanf(optarg, "%lf", &pso_w_max);
                break;
            case 's': /* PSO w_min */
                sscanf(optarg, "%lf", &pso_w_min);
                break;
            case 't': /* PSO w_strategy */
                sscanf(optarg, "%d", &pso_w_strategy_select);
                break;
            case 'w': /* PSO nhood_strategy */
                sscanf(optarg, "%d", &pso_nhood_topology_select);
                break;
	    case 'x': /*PSO nhood_size */
		sscanf(optarg, "%d", &pso_nhood_size);
		break;
	    case 'z': //PopSize
		popSize = 100;
		break;
	    default:
		abort();
}


void options() {
    printf ("Dimension = (%f,%f), Start = (%f,%f), Target = (%f,%f)\n", 
            inHorizonX, inHorizonY, inStartX, inStartY, inEndX, inEndY);
    printf ("Map File = %s\n", inFileHandlePtr);
    printf ("PSO: c1 = %f, c2 = %f, weight strategy = %d, neighborhood topology = %d\n", 
            pso_c1, pso_c2, pso_w_strategy, pso_nhood_topology);
    if (pso_w_strategy == PSO_W_LIN_DEC)
        printf ("\tweight min = %f, weight max = %f\n", pso_w_min, pso_w_max);
    if (pso_nhood_topology_select == PSO_NHOOD_RANDOM)
        printf("\tneighborhood size = %d\n", pso_nhood_size);
}



void pso_set_path_settings(pso_settings_t *settings, pso_params_t *params, env_t *env, robot_t *robot, int waypoints) {
    /* WARNING */
    // Only valid if a square environment with same start and stop 
    // EX: (0, 0) to (100, 100) because the pso lib
    // only considers each as an 'x-value', not knowing that
    // we are using a vector where odds are 'x' and evens are 'y'
    settings->x_lo = env->mins[0];
    settings->x_hi = env->maxs[0];

    settings->limits = pso_autofill_limits(settings->x_lo, settings->x_hi, settings->dim);

    settings->dim = waypoints * 2;
    // Set odd values limits using X-min and X-max
    // and even values limits using Y-min and Y-max
    int i;
    int count = 0;
    for (i = 0; i < waypoints; i++) {
        settings->limits[0][count] = env->mins[0];
        settings->limits[1][count] = env->maxs[0];
        count++;
        settings->limits[0][count] = env->mins[1];
        settings->limits[1][count] = env->maxs[1];
        count++;
    }

    pso_print_limits(settings->limits, settings->dim);

    // Set parameters, if not default
    if (params->c1 >= 0)
        settings->c1 = params->c1;
    if (params->c2 >= 0)
        settings->c2 = params->c2;
    if (params->w_min >= 0)
        settings->w_min = params->w_min;
    if (params->w_max >= 0)
        settings->w_max = params->w_max;
    if (params->w_strategy >= 0)
        settings->w_strategy = params->w_strategy;
    if (params->nhood_size >= 0)
        settings->nhood_size = params->nhood_size;
    if (params->nhood_topology >= 0)
        settings->nhood_strategy = params->nhood_topology;

    settings->goal = 1e-5;
    settings->numset = INTEGER;
}


void print_usage() {


}

