
#include <stdbool.h>

#ifndef _UTILS_H
#define _UTILS_H

/*Initial PSO settings */
extern int popSize; 
extern int maxIterations; 

/* Serial and parallel option */
extern bool demo; //for benchmark functions
extern bool serial; //for serial path example
//extern int timing;
/* Path options */
extern int inUavID;
extern double inStartX, inStartY, inEndX, inEndY;
extern double inStepSize, inVelocity;
extern double inOriginX, inOriginY, inHorizonX, inHorizonY;  // 70
extern char inFileHandle[];
//extern char inFileHandle[] = "sample_map_OpenRooms.txt";
extern int waypoints;

/* PSO parameters */
extern double pso_c1, pso_c2, pso_w_max, pso_w_min;
extern int pso_w_strategy_select, pso_nhood_size, pso_nhood_topology_select;
extern int pso_w_strategy, pso_nhood_topology;

/* Option parsing */
extern int verbose;
extern char *inFileHandlePtr;
int parse_arguments(int argc, char **argv);

//TIMER
struct Timer{
	clock_t start;
	clock_t finish;
	float duration_ms;
};

struct timing_report{
	struct Timer demo_time;
	struct Timer serial_time;
};

void start_timer(struct Timer* timing);
void end_timer(struct Timer* timing);
double elapsed_time(clock_t start, clock_t finish);
void print_elapsed_time(char* fn_name, clock_t start, clock_t finish);

//void options();
//void print_usage();

#endif
