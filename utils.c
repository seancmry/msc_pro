#include <stdio.h>
#include <getopt.h>
#include <stdbool.h>
#include <time.h>
#include "utils.h"

/*Initial PSO settings */
//int popSize = 100;
//int maxIterations = 500; 

/* Serial and parallel option */
bool serial = true;
//bool parallel = 0;
bool demo = true; //for benchmark functions
int timing = 0;

/* Path options */
int inRoboID = 0;
double inStartX = 70.0;
double inStartY = 70.0;
double inEndX = 136.0;
double inEndY = 127.0;
//double inStepSize = 1;
//double inVelocity = 2;
//double inOriginX = 0;
//double inOriginY = 0;
double inHorizonX = 200;
double inHorizonY = 200;  // 70
//char inFileHandle[20] = "maps/sampleMap4.dat\0";
char inFileHandle[] = "sample_map_OpenRooms.txt";
int waypoints = 5;

/* PSO parameters */
//double pso_c1 = -1.0;
//double pso_c2 = -1.0;
//double pso_w_max = -1.0;
//double pso_w_min = -1.0;
//int pso_w_strategy_select = -1;
//int pso_nhood_size = -1;
//int pso_nhood_topology_select = -1;

//int pso_w_strategy = -1;
//int pso_nhood_topology = -1;

/* Option parsing */
int verbose = 0;
//char *inFileHandlePtr = NULL;


int parse_arguments(int argc, char **argv) {
    int c;
    
    while ((c = getopt (argc, argv, "v:z:a:b:h:i:d:e:f:n:m:s:t")) != -1)
        switch (c) {
            case 'v':
                verbose = 1;
                break;
            case 'z':
		serial = 1;
		break;
	    //case 'g':
		//parallel = 1;
		//break;
            case 'a':
                sscanf(optarg, "%lf", &inHorizonX);
                break;
            case 'b':
                sscanf(optarg, "%lf", &inHorizonY);
                break;
            case 'h':
                sscanf(optarg, "%lf", &inStartX);
                break;
	    case 'i':
		sscanf(optarg, "%lf", &inStartY);
		break;
            case 'd':
                demo = false;
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
            //case 'p': /* PSO c1 */
                //sscanf(optarg, "%lf", &pso_c1);
                //break;
            //case 'q': /* PSO c2 */
                //sscanf(optarg, "%lf", &pso_c2);
                //break;
            //case 'r': /* PSO w_max */
                //sscanf(optarg, "%lf", &pso_w_max);
                //break;
            case 's': /* PSO w_min */
                serial = false;
		//sscanf(optarg, "%lf", &pso_w_min);
                break;
            case 't': //timing  /* PSO w_strategy */
                //sscanf(optarg, "%d", &pso_w_strategy_select);
                timing = 1;
		break;
            //case 'w': /* PSO nhood_strategy */
                //sscanf(optarg, "%d", &pso_nhood_topology_select);
                //break;
	    //case 'x': /*PSO nhood_size */
		//sscanf(optarg, "%d", &pso_nhood_size);
		//break;
	    //case 'z': //PopSize
		//popSize = 100;
		//break;
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


int **readMap(char * fhandle, int height, int width) {

    FILE *file;
    file = fopen(fhandle, "r");
    size_t count; 
    char *line = (char *) malloc (sizeof (char) * width + 1);
    
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


/*
void options() {
    printf ("Dimension = (%f,%f), Start = (%f,%f), Target = (%f,%f)\n", 
            inHorizonX, inHorizonY, inStartX, inStartY, inEndX, inEndY);
    printf ("Map File = %s\n", inFileHandlePtr);
    printf ("PSO: c1 = %f, c2 = %f, weight strategy = %d, neighborhood topology = %d\n", 
           pso_c1, pso_c2, pso_w_strategy, pso_nhood_topology);
    if (pso_w_strategy == PSO_W_LIN_DEC)
        printf ("\tweight min = %f, weight max = %f\n", pso_w_min, pso_w_max);
    if (pso_nhood_topology_select == PSO_NHOOD_RANDOM)
        printf("\tneighborhood size = %d\n", pso_nhood_size);
}
*/


/*
void print_usage() {
return;

}
*/

