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
bool serial = true;
//bool parallel = 0;
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
    
    while ((c = getopt (argc, argv, "a:b:c:d:e:f:n:m:p:q:r:s:t:w:x:v:z")) != -1)
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
	return 0;
}


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


//For each Cartesian coord
void calculate_dims(int nproc, int* dims){
	int root = (int)sqrt(nproc);
	while(nproc % root != 0)
		root--;
	dims[0] = nproc/root;
	dims[1] = root;
}

//1d decomposition of n into size procs
int decomp1d(int n, int size, int rank, int *s, int *e){
	int nlocal, deficit;
	nlocal  = n / size;
	*s  = rank * nlocal + 1;
	deficit = n % size;
	*s  = *s + ((rank < deficit) ? rank : deficit);
	if (rank < deficit) nlocal++;
		*e = *s + nlocal - 1;

	if (*e > n || rank == size-1) *e = n;
		
	return MPI_SUCCESS;
}

// 2d decomposition of array into xprocs x yprocs 
int decomp2d(int nx, int ny, int xproc, int yproc, int* coords, int *xs, int *xe, int *ys, int *ye){

	decomp1d(nx, xproc, coords[0], xs, xe);
	decomp1d(ny, yproc, coords[1], ys, ye);
	return MPI_SUCCESS;
}


// intialises array ptr
void init_arr(int n, int m, double *x, double **x_ptr){
	int i;
	for(i=0;i<n;i++){
		x_ptr[i] = &x[i*m];
	}
}

//Clear the array
void clear_arr(int n, int m, double **x){
    	int i, j;
    	for(i=0;i<n;i++){
        	for(j=0;j<m;j++){
           		x[i][j] = 0.0;
		}
	}
}

//Initialise the mesh with boundary conditions
void init_range(double **unew, double **uold, double **f, int xs, int xe, int ys, int ye, int nx, int ny,
	double (*lbound)(int, int, int, int),
	double (*rbound)(int, int, int, int),
	double (*ubound)(int, int, int, int),
	double (*bbound)(int, int, int, int))
	{	
	int i;

	//lower boundary
	if (ys == 1){
		for (i=(xs-1);i<=(xe+1);i++){
			uold[i][0] = bbound(i,0,nx,ny);
	 		unew[i][0] = bbound(i,0,nx,ny);
	 	}
	}

	//Upper boundary
	if (ye == ny){
		for (i=(xs-1);i<=(xe+1);i++){
			uold[i][ny+1] = ubound(i,ny+1,nx,ny);
			unew[i][ny+1] = ubound(i,ny+1,nx,ny);
		}
	}

	//Left boundary
	if (xs == 1){
		for(i=ys;i<=ye;i++){
			uold[0][i] = lbound(0,i,nx,ny);
			unew[0][i] = lbound(0,i,nx,ny);
		}
	}

	//Right boundary
	if (xe == nx){
		for(i=ys;i<=ye;i++){
			uold[nx+1][i] = rbound(nx+1,i,nx,ny);
			unew[nx+1][i] = rbound(nx+1,i,nx,ny);
		}
	}
}

























