#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "funcs.h"
//#include "mpi.h"

int main(){

	//Allocate buffer space for second receiver list
	list_a_t *first = NULL;
	list_b_t second;

	first = list_new(30);	
	
	//Space for solution buffer
	second.solution = (double *)malloc(first->dim * sizeof(double));

	//RUN ALG
	list(first, &second);	

	free(second.solution);
	
	//Free first struct	
	free(first);
	
	return 0;
}
