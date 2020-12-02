
#include "utils.h"

#ifndef PSO_H_
#define PSO_H_


// CONSTANTS
#define PSO_MAX_SIZE 100 // max swarm size
#define PSO_INERTIA 0.7298 // default value of w (see clerc02)


// === NEIGHBORHOOD SCHEMES ===

// global best topology
#define PSO_NHOOD_GLOBAL 0

// ring topology
#define PSO_NHOOD_RING 1

// Random neighborhood topology
// **see http://clerc.maurice.free.fr/pso/random_topology.pdf**
#define PSO_NHOOD_RANDOM 2



// === INERTIA WEIGHT UPDATE FUNCTIONS ===
#define PSO_W_CONST 0
#define PSO_W_LIN_DEC 1

// === DISCRETIZATION FLAGS ===
#define INTEGER 0
#define DECIMAL 1

// PSO SOLUTION -- Initialized by the user
typedef struct{

	double error;
	double *gbest; // should contain DIM elements!!

}pso_result_t;


// OBJECTIVE FUNCTION TYPE
typedef double (*pso_obj_fun_t)(double *, int, void *);

/*
 * pso_settings_t *pso_settings_new(int dim, double range_lo, double range_hi);
 * void pso_settings_free(pso_settings_t *settings);
 */

// set x value limits using two constants
double **pso_autofill_limits (double x_lo, double x_hi, int dim);
// print those limits
void pso_print_limits (double ** limits, int dim);

// return the swarm size based on dimensionality
int pso_calc_swarm_size(int dim);

// set the default PSO settings
void pso_set_default_settings(pso_settings_t *settings);

// minimize the provided obj_fun using PSO with the specified settings
// and store the result in *solution
void pso_solve(pso_obj_fun_t obj_fun, void *obj_fun_params, pso_result_t *solution, pso_settings_t *settings);

#endif // PSO_H_
