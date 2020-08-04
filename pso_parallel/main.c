
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h> // for printf
#include <getopt.h>
#include <sys/time.h>
#include "pso.h"
#include "mpi.h"

//==============================================================
//                  BENCHMARK FUNCTIONS
//==============================================================

double pso_sphere(double *vec, int dim, void *params) {

    double sum = 0;
    int i;
    for (i=0; i<dim; i++)
        sum += pow(vec[i], 2);

    return sum;
}



double pso_rosenbrock(double *vec, int dim, void *params) {

    double sum = 0;
    int i;
    for (i=0; i<dim-1; i++)
        sum += 100 * pow((vec[i+1] - pow(vec[i], 2)), 2) +	\
            pow((1 - vec[i]), 2);

    return sum;

}


double pso_griewank(double *vec, int dim, void *params) {

    double sum = 0.;
    double prod = 1.;
    int i;
    for (i=0; i<dim;i++) {
        sum += pow(vec[i], 2);
        prod *= cos(vec[i] / sqrt(i+1));
    }

    return sum / 4000 - prod + 1;

}



int main(int argc, char **argv) {

	pso_settings_t *settings = NULL;
        pso_obj_fun_t obj_fun = NULL;

        while((c = getopt(argc,argv,"n:m:t"))!=-1){
		switch(c){
			case 'n': N = atoi(optarg);
				  break;
	                case 'm': M = atoi(optarg);
			          break;
		        case 't': time = true;
			          break;
	                default:
				  printf("Incorrect options given\n");
			          return 1;
		}
	}
        if(argc != optind){
		printf("Too many options entered, exiting\n");
		return 1;
	}

    if (argc == 2) {
        if (strcmp(argv[1], "rosenbrock") == 0) {
            obj_fun = pso_rosenbrock;
            settings = pso_settings_new(1000, -2.048, 2.048);
            printf("Optimizing function: rosenbrock (dim=%d, swarm size=%d)\n",
                   settings->dim, settings->size);
        } else if (strcmp(argv[1], "griewank") == 0) {
            obj_fun = pso_griewank;
            settings = pso_settings_new(1000, -600, 600);
            printf("Optimizing function: griewank (dim=%d, swarm size=%d)\n",
                   settings->dim, settings->size);
        } else if (strcmp(argv[1], "sphere") == 0) {
            obj_fun = pso_sphere;
            settings = pso_settings_new(1000, -100, 100);
            printf("Optimizing function: sphere (dim=%d, swarm size=%d)\n",
                   settings->dim, settings->size);
        } else {
            printf("Unsupported objective function: %s", argv[1]);
            return 1;
        }
    } else if (argc > 2) {
        printf("Usage: demo [PROBLEM], where problem is optional with values [sphere|rosenbrock|griewank]\n ");
        return 1;
    }

    // handle the default case (no argument given)
    if (obj_fun == NULL || settings == NULL) {
        obj_fun = pso_sphere;
        settings = pso_settings_new(1000, -100, 100);
        printf("Optimizing function: sphere (dim=%d, swarm size=%d)\n",
                   settings->dim, settings->size);
    }

    // set some general PSO settings
    settings->goal = 1e-7;
    settings->size = 30;
    settings->nhood_strategy = PSO_NHOOD_RING;
    settings->nhood_size = 100;
    settings->w_strategy = PSO_W_LIN_DEC;

    // initialize GBEST solution
    pso_result_t solution;
    // allocate memory for the best position buffer
    solution.gbest = (double *)malloc(settings->dim * sizeof(double));

    // run optimization algorithm
    pso_solve(obj_fun, NULL, &solution, settings);

    // free the gbest buffer
    free(solution.gbest);

    // free the settings
    pso_settings_free(settings);

    return 0;

}

