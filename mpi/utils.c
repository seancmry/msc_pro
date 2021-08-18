#include <stdio.h>
#include <getopt.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>

#include "mpi.h"
#include "utils.h"


/*Initial PSO settings */
//int maxIterations = 500; 
int popSize = 100;
/* Serial and parallel option */
//bool serial = true;
//bool parallel = true;
//bool demo = true; //for benchmark functions
//int timing = 0;

/* Path options */
int inUavID = 0;
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
int waypoints = 5;

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
    
    while ((c = getopt (argc, argv, "a:b:c:d:e:f:n:m:p:q:r:s:t:w:x:v:z")) != -1){
        switch (c) {
            case 'v':
  		verbose = 1;
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
            case 't': //timing  /* PSO w_strategy */
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
		printf ("Invalid argumentation, fail\n");
		return -1;
		}
	}
	return 0;
}

//Timer for serial code
void start_timer(struct Timer* timing) {
	timing->start = clock();
}

void end_timer(struct Timer* timing) {
	timing->finish = clock();
}


double elapsed_time(clock_t start, clock_t finish) {
	double time_in_ms = (double)(finish - start) / (CLOCKS_PER_SEC/1000);
	return time_in_ms;
}


void print_elapsed_time(char* fn_name, clock_t start, clock_t finish) {
	printf("%s: %fms \n", fn_name, elapsed_time(start,finish));
}








