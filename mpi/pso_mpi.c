#include <stdio.h>
#include <stdlib.h>
#include <time.h> // for time()
#include <math.h> // for cos(), pow(), sqrt() etc.
#include <float.h> // for DBL_MAX
#include <string.h> // for mem*
#include <stdint.h> //for uintptr_t
#include <gsl/gsl_rng.h>
#include <mpi.h>
#include <omp.h>

#include "utils.h"
#include "pso.h"


#define N 9
#define DIMS 100 
#define MAXN 2000

// function type for the different inform functions
typedef void (*inform_fun_t)(int *comm, double **pos_nb, double **pos_b, double *fit_b, double *gbest, int improved, pso_settings_t *settings);



//function type for the different inertia calculation functions
typedef double (*inertia_fun_t)(int step, pso_settings_t *settings);



double roundNum (double d) {
  double g = 0.0, e = 0.0, f = 0.0;
  g = floor (d);
  e = d - g;
  if (d > 0){
    if (e > 0.5){
      f = g + 1;
    } else {
      f = g;
    }
  } else {
      if (e > 0.5) {
          f = g - 1;
      } else {
          f = g;
      }
  }
  return f;
}


//==============================================================
// calulate swarm size based on dimensionality
int pso_calc_swarm_size(int dim) {
  int size = 10. + 2. * sqrt(dim);
  return (size > PSO_MAX_SIZE ? PSO_MAX_SIZE : size);
}



//==============================================================
//          INERTIA WEIGHT UPDATE STRATEGIES
//==============================================================
// calculate linearly decreasing inertia weight
double calc_inertia_lin_dec(int step, pso_settings_t *settings) {

  int dec_stage = 3 * settings->steps / 4;
  if (step <= dec_stage)
    return settings->w_min + (settings->w_max - settings->w_min) *	\
      (dec_stage - step) / dec_stage;
  else
    return settings->w_min;
}


/*
// calculate linearly decreasing inertia weight
double MPI_calc_inertia_lin_dec(int step, pso_settings_t *settings) {

	//Get rank and size of original comm
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);	
	
	
	//We determine the decreased weight based on the row and split the communicator on that basis.
	//We then use the original rank for ordering.
	//What happens then is that the value gets stored in val[i] and distributed 
	//among all procs (as defined above = 4) through MPI_Scatter. As it ran with MPI_Alltoall, the programme
	//accounted for a total of w*4 weights for a program run across 4 procs in parallel. This 
	//will need to be adjusted so that the program is only run once across 4 procs and MPI_Scatter should do this.
	
	
	int dec_stage = 3 * settings->steps / 4;

	//Construct new comm
	MPI_Comm row_comm;
	MPI_Comm_split(MPI_COMM_WORLD, dec_stage, rank, &row_comm);
	int row_rank, row_size;
	MPI_Comm_rank(row_comm, &row_rank);
	MPI_Comm_size(row_comm, &row_size);

	//printf("WORLD RANK/SIZE: %d%d --- PRIME RANK/SIZE: %d%d\n", world_rank, world_size, row_rank, row_size);
	//MPI_Comm_free(&row_comm);
	int i;
	double val[N], recv_buf[N];
	for(i = 0; i < N; i++){
		if (step <= dec_stage) {
			val[i] = settings->w_min + (settings->w_max - settings->w_min) *     \
			(dec_stage - step) / dec_stage;
		} else {
			val[i] = settings->w_min;
		}
	}
	//Was originally executed as an MPI_Alltoall routine. Problem: I don't want a process 0.	
	MPI_Scatter(&val, 1, MPI_DOUBLE, recv_buf, 1, MPI_DOUBLE, row_comm);
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	if(MPI_COMM_NULL != row_comm) {
		MPI_Comm_free(&row_comm);
	}
	return recv_buf[i];
}
*/


//==============================================================
//                     NEIGHBOURHOOD
//==============================================================

void inform_global(int *comm, double **pos_nb, double **pos_b, double *fit_b,
	    double *gbest, int improved, pso_settings_t *settings){
	
	int i, nparticles;
  	for (i=0; i<nparticles; i++)
    	memmove((void *)pos_nb[i], (void *)gbest,
            	sizeof(double) * settings->dim);
}


// ===============================================================
// general inform function :: according to the connectivity
// matrix COMM, it copies the best position (from pos_b) of the
// informers of each particle to the pos_nb matrix
void inform(int *comm, double **pos_nb, double **pos_b, double *fit_b,
	    int improved, pso_settings_t *settings)
{
  	int i, j, nparticles;
  	int b_n; // best neighbor in terms of fitness

  	// for each particle
  	for (j=0; j<nparticles; j++) {
    		b_n = j; // self is best
    		// who is the best informer??
    		for (i=0; i<nparticles; i++){
      		// the i^th particle informs the j^th particle
      		if (comm[i*nparticles + j] && fit_b[i] < fit_b[b_n])
        		// found a better informer for j^th particle
        		b_n = i;
    		// copy pos_b of b_n^th particle to pos_nb[j]
    		memmove((void *)pos_nb[j],
            		(void *)pos_b[b_n],
            		sizeof(double) * settings->dim);
		}
  	}
}



//==============================================================
//                     PSO SETTINGS
//==============================================================

// set x value limits (high and low) using two constants
double ** pso_autofill_limits (double x_lo, double x_hi, int dim){
    double ** limits = (double **) malloc (sizeof (double *) * 2);
    // limits[0]: lower bound
    limits[0] = (double *) malloc (sizeof (double) * dim);
    // limits[1]: higher bound
    limits[1] = (double *) malloc (sizeof (double) * dim);
    int i = 0;
    for (i = 0; i < dim; i++) {
        limits[0][i] = x_lo;
        limits[1][i] = x_hi;
    }
    return limits;
}

void pso_print_limits (double ** limits, int dim){
    int i = 0;  
    for (i = 0; i < dim; i++) {
        printf ("x%d: lower = %f, higher = %f\n", i+1, limits[0][i], limits[1][i]);
    }

}


//==============================================================
// return default pso settings
pso_settings_t *pso_settings_new(int dim, double r_lo, double r_hi) {

	pso_settings_t *settings = (pso_settings_t *)malloc(sizeof(pso_settings_t));
	if (settings == NULL) {return NULL;}
	
  	// set some default values
  	settings->dim = dim;
  	settings->goal = 1e-5;
			

	settings->r_lo = (double *)malloc(settings->dim * sizeof(double));
	if (settings->r_lo == NULL) {free(settings); return NULL;}

	settings->r_hi = (double *)malloc(settings->dim * sizeof(double));
	if (settings->r_hi == NULL) {free(settings); free(settings->r_lo); return NULL;}

	int i;
	for (i=0; i<settings->dim; i++) {
		settings->r_lo[i] = r_lo;
		settings->r_hi[i] = r_hi;
	}
	
  	settings->size = pso_calc_swarm_size(settings->dim);
  	settings->print_every = 100;
  	settings->steps = 100000;
  	settings->c1 = 1.496;
  	settings->c2 = 1.496;
  	settings->w_max = PSO_INERTIA;
  	settings->w_min = 0.3;

	settings->numset = DECIMAL;

  	settings->clamp_pos = 1;
 	settings->nhood_size = 5;
  	settings->w_strategy = PSO_W_LIN_DEC;

	settings->rng = NULL;
	settings->seed = time(0);  
 	return settings;

}


void pso_serial_settings(pso_settings_t *settings){
	
	settings->dim = 100;
	settings->x_hi = 20;
	settings->x_lo = -20;
	settings->goal = 1e-5;
	settings->limits = pso_autofill_limits(settings->x_lo, settings->x_hi, settings->dim);

  	settings->size = pso_calc_swarm_size(settings->dim);
  	settings->print_every = 10;
  	settings->steps = 10001;
  	settings->c1 = 1.496;
  	settings->c2 = 1.496;
  	settings->w_max = PSO_INERTIA;
  	settings->w_min = 0.3;

	settings->numset = DECIMAL;

  	settings->clamp_pos = 1;
 	settings->nhood_size = 5;
  	settings->w_strategy = PSO_W_LIN_DEC;

	settings->rng = NULL;
	settings->seed = time(0);

}


//destroy PSO settings
void pso_settings_free(pso_settings_t *settings) {
	free(settings->r_lo);
	free(settings->r_hi);
	free(settings);
}

//OLD MATRIX ALLOCATION FUNCTIONS
double **pso_matrix_new(int size, int dim){
	double **m = (double **)malloc(size *sizeof(double *));
	int i;
	for (i = 0; i<size; i++) {
		m[i] = (double *)malloc(dim * sizeof(double));
	}
	return m;
}

void pso_matrix_free(double **m, int size){ 
	int i;
	for(i=0; i<size; i++) {
		free(m[i]);
	}
	free(m);
}


//==============================================================
//                     PSO ALGORITHM
//==============================================================
void pso_solve(pso_obj_fun_t obj_fun, void *obj_fun_params, pso_result_t *solution, pso_settings_t *settings)
{


	//MPI Settings
	int nparticles;
	int rank;
	int nproc = N;	
	int *recvbuf = (int *)malloc((settings->dim+1) * nproc * sizeof(int));
	int *sendbuf = (int *)malloc((settings->dim+1) * sizeof(int));
	int min;

	//Get rank
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	// Split the number of particles across processes	
	if(rank == 0){
		nparticles = (int)settings->size / nproc;
		//printf("Number of particles is %d\n", nparticles);
	}		
	MPI_Bcast(&nparticles, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if(rank == 0){
		nparticles += (int)settings->size % nproc;
		//printf("Number of particles is %d\n", nparticles);
	}

	//RNG	
	int free_rng = 0;
  	// Particles
  	double **pos = pso_matrix_new(nparticles, settings->dim); // position matrix
  	double **vel = pso_matrix_new(nparticles, settings->dim) ; // velocity matrix
  	double **pos_b =  pso_matrix_new(nparticles, settings->dim); // best position matrix
  	double *fit = (double *)malloc(nparticles * sizeof(double)); // particle fitness vector
  	double *fit_b = (double *)malloc(nparticles * sizeof(double)) ; // best fitness vector
 	// Swarm
  	double **pos_nb = pso_matrix_new(nparticles, settings->dim); // what is the best informed
  	int *comm = (int *)malloc(nparticles * nparticles * sizeof(int)); // communications:who informs who
                                            // rows : those who inform
                                            // cols : those who are informed
	int improved = 0; // whether solution->error was improved during the last iteration
  	int i, d, step;
  	double a, b; // for matrix initialization
  	double rho1, rho2; // random numbers (coefficients)
  	double w = PSO_INERTIA; // current omega
  	inform_fun_t inform_fun = NULL; // neighborhood update function
  	inertia_fun_t calc_inertia_fun = NULL; // inertia weight update function
		

  	// CHECK RANDOM NUMBER GENERATOR
  	if (!settings->rng) {
    		// initialize random number generator
    		gsl_rng_env_setup();
    		// allocate the RNG
    		settings->rng = gsl_rng_alloc(gsl_rng_default);
    		// seed the generator
    		gsl_rng_set(settings->rng, settings->seed);
    		// remember to free the RNG - see free at the end
    		free_rng = 1;
  	}

	//Set inform function to global
	inform_fun = inform_global;	 


  	// SELECT APPROPRIATE INERTIA WEIGHT UPDATE FUNCTION
  	switch (settings->w_strategy) {
      	/* case PSO_W_CONST : */
     	/*     calc_inertia_fun = calc_inertia_const; */
      	/*     break; */
    		case PSO_W_LIN_DEC :
	 		//This is where the MPI calculation is done for the w_strategy.		
      			//calc_inertia_fun = MPI_calc_inertia_lin_dec;	
      			calc_inertia_fun = calc_inertia_lin_dec; 
			break;
    	}

  	// INITIALIZE SOLUTION
  	solution->error = DBL_MAX;
	

	/* START */



  	// SWARM INITIALIZATION
  	// for each particle
  	#pragma omp parallel for private(a,b) reduction(min:solution->error)
  	for (i=0; i<nparticles; i++) {
    		// for each dimension
    		for (d=0; d<settings->dim; d++) {
			// generate two numbers within the specified range
			//if (demo){
				a = settings->r_lo[d] + (settings->r_hi[d] - settings->r_lo[d])  *   \
				gsl_rng_uniform(settings->rng);
       				b = settings->r_lo[d] + (settings->r_hi[d] - settings->r_lo[d])  *      \
				gsl_rng_uniform(settings->rng);
       			//}
			//if (serial){
				//a = gsl_rng_uniform_int(settings->rng, settings->limits[1][i]);
		        	//b = gsl_rng_uniform_int(settings->rng, settings->limits[1][i]);
			//}
			/*
			a = settings->limits[0][i] + (settings->limits[1][i] - settings->limits[0][i]) 
   			 gsl_rng_uniform(settings->rng);
      			b = settings->limits[0][i] + (settings->limits[1][i] - settings->limits[0][i])
       			 gsl_rng_uniform(settings->rng);
			*/
      			// initialize position
      			pos[i][d] = a;
      			// best position is the same
      			pos_b[i][d] = a;
      			// initialize velocity
      			vel[i][d] = (a-b) / 2.;
    		}
    		
		
			

	
	
		// update particle fitness
    		fit[i] = obj_fun(pos[i], settings->dim, obj_fun_params);
   		fit_b[i] = fit[i]; // this is also the personal best
    




		// update gbest??
    		if (fit[i] < solution->error) {
			// update best fitness
     			solution->error = fit[i];
      	
			// copy particle pos to gbest vector
  			memmove((void *)solution->gbest, (void *)pos[i],
				sizeof(double) * settings->dim);
		}
    	}

	
	/*END*/


	
	// initialize omega using standard value
	w = PSO_INERTIA;



	
	// RUN ALGORITHM
  	for (step=0; step<settings->steps; step++) {
		#pragma omp parallel num_threads(4) shared(min)
		// update current step
    		settings->step = step;
    		// update inertia weight
    		// do not bother with calling a calc_w_const function
    		if (calc_inertia_fun != NULL) {
      			w = calc_inertia_fun(step, settings);
		}



		/*
    		// check optimization goal
    		if (solution->error <= settings->goal) {
     	 		// SOLVED!!
      			if (settings->print_every){
        			printf("Goal achieved @ step %d (error=%.3e) :-)\n", step, solution->error);
			}
     			break;
    		}	
		*/



    		//update pos_nb matrix (find best of neighborhood for all particles)
    		//inform_fun would go here
  		// copy the contents of gbest to pos_nb
  		//inform_fun(comm, (double **)pos_nb, (double **)pos_b, fit_b, solution->gbest,
		//			improved, settings);


    		//the value of improved was just used; reset it
    		improved = 0;




    		// update all particles
    		for (i=0; i<nparticles; i++) {
      			
			#pragma omp for private(a,b)
			// for each dimension
      			for (d=0; d<settings->dim; d++) {
        		// calculate stochastic coefficients
        			rho1 = settings->c1 * gsl_rng_uniform(settings->rng);
        			rho2 = settings->c2 * gsl_rng_uniform(settings->rng);


			// update velocity
				vel[i][d] = w * vel[i][d] +    \
          				rho1 * (pos_b[i][d] - pos[i][d]) +	\
          				rho2 * (pos_nb[i][d] - pos[i][d]);
        		// update position
        			pos[i][d] += vel[i][d];
				
				if (settings->numset == INTEGER){
					pos[i][d] = roundNum(pos[i][d]);
				}
				
				//if(demo){
        				// clamp position within bounds?
        				if (settings->clamp_pos) {
          					if (pos[i][d] < settings->r_lo[d]) {
            						pos[i][d] = settings->r_lo[d];
            						vel[i][d] = 0;
          					} else if (pos[i][d] > settings->r_hi[d]) {
           		 				pos[i][d] = settings->r_hi[d];
            						vel[i][d] = 0;
						}
        				} else {
          				// enforce periodic boundary conditions
          					if (pos[i][d] < settings->r_lo[d]) {
            						pos[i][d] = settings->r_hi[d] - fmod(settings->r_lo[d] - pos[i][d],
                                              			settings->r_hi[d] - settings->r_lo[d]);
            						vel[i][d] = 0;
          					} else if (pos[i][d] > settings->r_hi[d]) {
            						pos[i][d] = settings->r_lo[d] + fmod(pos[i][d] - settings->r_hi[d],
                                              			settings->r_hi[d] - settings->r_lo[d]);
            						vel[i][d] = 0;
          					}
        				}
				//}
				/*
				//if(serial) {
					if (settings->clamp_pos) {
          					if (pos[i][d] < settings->limits[0][i]) {
            						pos[i][d] = settings->limits[0][i];
            						vel[i][d] = 0;
          					} else if (pos[i][d] > settings->limits[1][i]) {
           		 				pos[i][d] = settings->limits[1][i];
            						vel[i][d] = 0;
						}
        				} else {
          				// enforce periodic boundary conditions
          					if (pos[i][d] < settings->limits[0][i]) {
            						pos[i][d] = settings->limits[1][i] - fmod(settings->limits[0][i] - pos[i][d],
                                              			settings->limits[1][i] - settings->limits[0][i]);
            						vel[i][d] = 0;
          					} else if (pos[i][d] > settings->limits[1][i]) {
            						pos[i][d] = settings->limits[0][i] + fmod(pos[i][d] - settings->limits[1][i],
                                              			settings->limits[1][i] - settings->limits[0][i]);
            						vel[i][d] = 0;
          					}
        				}
				//}
				*/
					
			}
                	// update particle fitness
                	fit[i] = obj_fun(pos[i], settings->dim, obj_fun_params);




                	// update personal best position?
      	        	if (fit[i] < fit_b[i]) {
        			fit_b[i] = fit[i];
        		
			// copy contents of pos[i] to pos_b[i]
        		memmove((void *)pos_b[i], (void *)pos[i],
                		sizeof(double) * settings->dim);
      			}

		

      			// update gbest?? - algorithm would end here in serial version,
      			// but additional code is required to perform the same update
      			// in each proc
      			
			#pragma omp for reduction(min:solution->error)
      			for(i=0; i<(int)nparticles; i++){		
				if (fit[i] < solution->error) {
        				//improved = 1;
        				// update best fitness
        				solution->error = fit[i];
        				// copy particle pos to gbest vector - removed because we have yet to find gbest
      				}
			}
			
			//This is further split for the additon of omp code to both parts
			#pragma omp for
			for (i=0; i<(int)nparticles; i++){
				if (solution->error == fit[i]){
					min = i;
				}	
			}		
		}
		//Update gbest with min position from personal best positions on each process - will then be gathered up on rank 0
		memmove((void *)solution->gbest, (void *)pos_b[min], sizeof(double) * settings->dim);			
		

		
		//Store global best position in the send buffer
		for (d=0; d<(int)settings->dim; d++){
			sendbuf[d] = (int)solution->gbest[d];
		}
			
		sendbuf[(int)settings->dim] = solution->error;
		
		
		//Gather data	
		MPI_Gather(&sendbuf, (settings->dim+1), MPI_INT, recvbuf, (settings->dim+1), MPI_INT, 0, MPI_COMM_WORLD);
 
 
		//Pass best position to recvbuf
		if (rank == 0){
			min = solution->error;
			int p = -1; //denote position
			int k;
			for (k = 0; k<nproc; k++){	
				if(min >= recvbuf[k*((int)settings->dim+1)+((int)settings->dim)]){
					min = recvbuf[k*((int)settings->dim+1)+((int)settings->dim)];
					p = k*((int)settings->dim+1);
				}
			}
			solution->error = min;
			k = 0;
			for(k = p; k<(int)settings->dim+p;k++){
				solution->gbest[k-p] = recvbuf[k];
			}
		}
	
		//Broadcast best position
		MPI_Bcast(solution->gbest, settings->dim, MPI_INT, 0, MPI_COMM_WORLD);

		//Print from rank 1
		if(rank == 1){		
			if (settings->print_every && (step % settings->print_every == 0)){
      				printf("Rank: %d, Step %d,    w=%.2f,    min_err=,    %.5e\n", rank, step, w, solution->error);	
			}
		}
	}

	 

	//free resources
	//if (serial){ 

		pso_matrix_free(pos, nparticles);
		pso_matrix_free(vel, nparticles);
		pso_matrix_free(pos_b, nparticles);
		pso_matrix_free(pos_nb, nparticles);
		free(comm);
		free(fit);
		free(fit_b);
	//}	
	
	if (free_rng)
		gsl_rng_free(settings->rng);
		
	if (rank == 0){
		free(sendbuf);
		free(recvbuf);
	}
	

}
