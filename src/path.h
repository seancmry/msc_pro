#include <math.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdbool.h>
#include <ctype.h>

typedef struct {
    	int ID; 
    	double position_coords[2];
    	double target_coords[2];
    	double distToTarget;
    	double stepSize;
    	double velocity;
} uav_t; 

typedef struct {
    	double *waypoints;
    	int numWaypoints;
} path_t;

typedef struct {
    	double mins[2];
    	double maxs[2];
    	int ** map;
} env_t;

typedef struct {
    	double start[2];
    	double stop[2];
    	double c1;
    	double c2;
    	int w_strategy;
    	double w_min;
    	double w_max;
    	int nhood_topology;
    	double nhood_size;
    	env_t * env;
} pso_params_t;


float distance_from_obstacle(/*coord inputs, */ float centre_x, float centre_y, float);
int **read_map(char *fhandle, int height, int width, int length);
int *rules(int x, int y, char* prevIter);
void print_map (int **map, int height, int width);
double euclideanDistance(double xi, double yi, double xj, double yj);
bool circle(float x, float y, float centre_x, float centre_y);

uav_t * init_uav(int ID, double xInit, double yInit, double xTarget, double yTarget, double stepSize, double velocity);
env_t *init_env(double xMin, double yMin, double xMax, double yMax, int **map);

void print_uav(uav_t *uav);
void print_env(env_t *env);
int line2 (int x0, int y0, int x1, int y1, int ** map, int xLimit, int yLimit);
int pso_path_countObstructions(double *vec, int dim, void *params);
double pso_path_penalty(double *vec, int dim, void *params);
double pso_path(double *vec, int dim, void *params); 
void pso_set_path_settings(pso_settings_t *settings, pso_params_t *params, env_t *env, robot_t *robot, int waypoints);


int getPSOParam_w_stategy(int code);
int getPSOParam_nhood_topology(int code);



