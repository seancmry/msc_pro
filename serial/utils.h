#ifndef _UTILS_H_
#define _UTILS_H_

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



//write timings of different runs to file




//print usages
void print_usage() {




}
#endif
