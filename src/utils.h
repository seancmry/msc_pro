#ifndef _UTILS_H_
#define _UTILS_H_

//int NUM_OBSTRUCT = 10

extern int verbose, inRoboId;
extern double inHorizonX, inHorizonY, inStartX, inStartY, inEndX, inEndY, inStepSize, inVelocity, inOriginX, inOriginY;
extern char inFileHandle[];
extern int waypoints;

/* PSO parameters */
extern double pso_c1, pso_c2, pso_w_max, pso_w_min;
extern int pso_w_strategy_select, pso_nhood_size, pso_nhood_topology_select, pso_w_strategy, pso_nhood_topology;

//Function defs
int parse_args(int argc, char **argv);


#endif
