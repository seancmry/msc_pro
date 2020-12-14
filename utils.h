
#ifndef _UTILS_H
#define _UTILS_H

/*Initial PSO settings */
extern int popSize; 
//extern int maxIterations = 500; 

/* Serial and parallel option */
extern int demo; //for benchmark functions
//extern bool serial, parallel;
extern int timing;
/* Path options */
//extern int inRoboID;
//extern double inStartX, inStartY, inEndX, inEndY;
//extern double inStepSize, inVelocity;
//extern double inOriginX, inOriginY, inHorizonX, inHorizonY;  // 70
//char inFileHandle[20] = "maps/sample_map_OpenRooms.txt\0";
//extern char inFileHandle[] = "sample_map_OpenRooms.txt";
//extern int waypoints;

/* PSO parameters */
//extern double pso_c1, pso_c2, pso_w_max, pso_w_min;
//extern int pso_w_strategy_select, pso_nhood_size, pso_nhood_topology_select;
//extern int pso_w_strategy, pso_nhood_topology;

/* Option parsing */
extern int verbose;
extern char *inFileHandlePtr;


//TIMER
struct Timer{
	clock_t start;
	clock_t finish;
	float duration_ms;
};

struct timing_report{
	struct Timer demo_time;
};

void print_elapsed_time(char* fn_name, clock_t start, clock_t finish);
double elapsed_time(clock_t start, clock_t finish);
void start_timer(struct Timer* timing);
void end_timer(struct Timer* timing);

int parse_arguments(int argc, char **argv);
//void options();
//void pso_set_path_settings(pso_settings_t *settings, pso_params_t *params, env_t *env, robot_t *robot, int waypoints);
//void print_usage();

#endif
