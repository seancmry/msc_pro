#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include "pso.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdbool.h>
#include <ctype.h>

#include "path.h"


int ** readMap(char * fhandle, int height, int width) {

    FILE *file;
    file = fopen(fhandle, "r");
    size_t count; 
    char *line = (char *)malloc(sizeof(char) * (width + 1));
    if (line == NULL) {
	printf("unable to allocate buffer\n");
	exit(1);
    }
    int i = 0, j = 0, ylim = 0;
    int ** map = (int **) malloc (sizeof (int *) * height);
    while (getline (&line, &count, file) != -1 && ylim < height){
        map[i] = (int *) malloc (sizeof (int) * width);
        for (j = 0; j < width; j++) {
            map[i][j] = line[j] - '0';
        }
        i++;
        ylim++;
    }
    fclose(file);
    return map;
}

void printMap (int **map, int height, int width){
    int i = 0, j = 0;
    for (i = 0; i < height; i++){
        for (j = 0; j < width; j++){
            printf ("%d", map[i][j]);
        }
        printf ("\n");
    }

}

double euclideanDistance(double xi, double yi, double xj, double yj) {
    return pow( pow(xi - xj, 2) + pow(yi - yj, 2), 0.5);
}


int line2 (int x0, int y0, int x1, int y1, int ** map, int xLimit, int yLimit) { 
// Source: https://github.com/ssloy/tinyrenderer/wiki/Lesson-1:-Bresenham%E2%80%99s-Line-Drawing-Algorithm

    int count = 0;
    float t = 0.0;
	for (t=0.; t<1.; t+=.01) { 
        int x = x0*(1.-t) + x1*t; 
        int y = y0*(1.-t) + y1*t; 

		if (x < yLimit && y < xLimit) {
			if (map[x][y] > 0){
				count++;
			}
		}
    } 
	return count;
}


uav_t * initUav(int ID, double xInit, double yInit, double xTarget, double yTarget, double stepSize, double velocity) {
    uav_t *uav = (uav_t *) malloc (sizeof (uav_t) * 1);
    uav->ID = ID;
    uav->position_coords[0] = xInit;
    uav->position_coords[1] = yInit;
    uav->target_coords[0] = xTarget;
    uav->target_coords[1] = yTarget;
    uav->distToTarget = euclideanDistance(xInit, yInit, xTarget, yTarget);
    uav->stepSize = stepSize;
    // This is the velocity of the uav
    uav->velocity = velocity;
    return uav;
}

env_t * initEnv(double xMin, double yMin, double xMax, double yMax, int ** map) {
    env_t *env = (env_t *) malloc (sizeof (env_t) * 1);
    env->mins[0] = xMin;
    env->mins[1] = yMin;
    env->maxs[0] = xMax;
    env->maxs[1] = yMax;
    env->map = map;
    return env;
}

void printUav(uav_t *uav) {
    printf ("UAV %d is at (%f, %f).\n", uav->ID, uav->position_coords[0], uav->position_coords[1]);
    printf ("The goal is at (%f, %f).\n", uav->target_coords[0], uav->target_coords[1]);
    printf ("It moves with velocity %f and with step size %f.\n", uav->velocity, uav->stepSize);
    printf ("It is currently at distance %f from the target.\n", uav->distToTarget);
}

void printEnv(env_t *env) {
    printf ("The world is of size %f x %f.\n", (env->maxs[0] - env->mins[0]), (env->maxs[1] - env->mins[1]));
    printf ("The bottom-left point is (%f, %f).\n", env->mins[0], env->mins[1]);
    printf ("The top-right point is (%f, %f).\n", env->maxs[0], env->maxs[1]);
}

int pso_path_countObstructions(double *vec, int dim, void *params) {
	pso_params_t * parameters;
    parameters = (pso_params_t *) params;
    
    // Ensure even-length vector
    if (dim % 2 != 0) {
        printf ("Inside 'pso_path_countObstruction': solution vector must be even-length!\n");
        exit (1);
    }

    int penaltyCount = 0;
    int i, count = 0;
    double xi = vec[0];
    double yi = vec[1];
    double xj = 0, yj = 0;

	penaltyCount = penaltyCount + line2 (parameters->start[0], parameters->start[1], 
            xi, yi, parameters->env->map, parameters->env->maxs[0], parameters->env->maxs[1]);

    count = 2;
    int iterations = dim / 2 - 1;
    for (i=0; i<iterations;i++) {
        xj = vec[count];
        count++;
        yj = vec[count];
        count++;
        penaltyCount = penaltyCount + line2 (xi, yi, xj, yj, parameters->env->map, parameters->env->maxs[0], parameters->env->maxs[1]);
        xi = xj;
        yi = yj;
    }

    penaltyCount = penaltyCount + line2 (xj, yj, parameters->stop[0], parameters->stop[1], 
            parameters->env->map, parameters->env->maxs[0], parameters->env->maxs[1]);
    return penaltyCount * 10;
}



double pso_path_penalty(double *vec, int dim, void *params) {
    // Penalizes the presence of obstacles in the path
    int countObstructions = pso_path_countObstructions (vec, dim, params);
    return pow (countObstructions, 1);
}


double pso_path(double *vec, int dim, void *params) {
    pso_params_t *parameters;
    parameters = (pso_params_t *) params;

    // Ensure even-length vector
    if (dim % 2 != 0) {
        printf ("Inside 'pso_path': solution vector must be even-length!\n");
        exit (1);
    }

    double distance = 0.0;
    int i, count = 0;

    double xi = vec[0];
    double yi = vec[1];
    double xj, yj = 0;

    distance = distance + euclideanDistance(parameters->start[0], parameters->start[1], xi, yi);

    count = 2;
    int iterations = dim / 2 - 1;
    for (i=0; i<iterations;i++) {
        xj = vec[count];
        count++;
        yj = vec[count];
        count++;
        distance = distance + euclideanDistance(xi, yi, xj, yj);
        xi = xj;
        yi = yj;
    }

    distance = distance + euclideanDistance(xj, yj, parameters->stop[0], parameters->stop[1]);
	distance = distance + pso_path_penalty(vec, dim, params);

    return distance;
}


void pso_set_path_settings(pso_settings_t *settings, pso_params_t *params, env_t *env, uav_t *uav, int waypoints) {
     
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


int getPSOParam_w_strategy(int code){
    if (code == 0)
        return PSO_W_CONST;
    if (code == 1)
        return PSO_W_LIN_DEC;
    return code;
}

int getPSOParam_nhood_topology(int code){
    if (code == 0)
        return PSO_NHOOD_GLOBAL;
    if (code == 1)
        return PSO_NHOOD_RING;
    if (code == 2)
        return PSO_NHOOD_RANDOM;
    return code;
}


