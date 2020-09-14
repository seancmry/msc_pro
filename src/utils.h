#ifndef _UTILS_H_
#define _UTILS_H_

int NUM_OBSTRUCT = 10;


//Function defs
int parse_args(int argc, char **argv) {
	
	int c;
	while ((c = getopt(argc,argv,"x:y:z"))!= -1) {
		switch(c){
			case 'x':
				do thing; break;
			case 'y':
				do thing; break;
			case 'z':
				do thing; break;
			default:
				fprintf(stderr,"Invalid option given\n");
				print_usage();
				return -1;
		}
	}
	return 0;

}


    int c;
    opterr = 0;

    while ((c = getopt (argc, argv, "a:b:c:d:e:f:n:m:p:q:r:s:t:w:x:v")) != -1)
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
            case 't': /* PSO w_strategy */
                sscanf(optarg, "%d", &pso_w_strategy_select);
                break;
            case 'w': /* PSO nhood_strategy */
                sscanf(optarg, "%d", &pso_nhood_topology_select);
		break;
	    case 'x': /* PSO nhood_size */
		sscanf(optarg, "%d", &pso_nhood_size);
		break;
	    default:
		abort();
	}
}

//write timings of different runs to file




//print usages
void print_usage() {




}
#endif
