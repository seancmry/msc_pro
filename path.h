
#include <stdlib.h>
#include <stdio.h>

typedef struct {
    int ID; 
    double position_coords[2];
    double target_coords[2];
    double distToTarget;
    double stepSize;
    double velocity;
} robot_t;

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


int ** readMap (char * fhandle, int height, int width);
void printMap (int **map, int height, int width);

double euclideanDistance(double xi, double yi, double xj, double yj);
double ChessboardDistance(double xi, double yi, double xj, double yj);
 
robot_t * initRobot(int ID, double xInit, double yInit, double xTarget, double yTarget, double stepSize, double velocity);
env_t * initEnv(double xMin, double yMin, double xMax, double yMax, int ** map); 

void printRobot(robot_t *robot);
void printEnv(env_t *env);

int line2 (int x0, int y0, int x1, int y1, int ** map, int xLimit, int yLimit);
int pso_path_countObstructions(double *vec, int dim, void *params);
double pso_path_penalty(double *vec, int dim, void *params);
double pso_path(double *vec, int dim, void *params);
int getPSOParam_w_stategy(int code);
int getPSOParam_nhood_topology(int code);



