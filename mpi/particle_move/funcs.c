#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>

#include "mpi.h"
#include "funcs.h"


//Allocate matrices
double **matrix_new(int size, int dim){
	int i;
	double **m = (double **)malloc(size *sizeof(double *));
	for (i = 0; i<size; i++) {
		m[i] = (double *)malloc(dim * sizeof(double));
	}
	return m;
}

void matrix_free(double **m, int size){ 
	int i;
	for(i=0; i<size; i++) {
		free(m[i]);
	}
	free(m);
}




list_a_t *list_new(int dim) {

	list_a_t *first = (list_a_t *)malloc(sizeof(list_a_t));
	if (first == NULL) {return NULL;}
	
  	// set some default values
  	first->dim = dim;
	first->size = 20;
	return first;
}	




void list(list_a_t *first, list_b_t *second){

	int i,j;
	double **a = matrix_new(first->size, first->dim);
	
	srand(time(NULL));

	second->error = DBL_MAX;

	//Serial	
	for (i=0;i<first->size;i++){
		for(j=0;j<first->dim;j++){
			a[i][j] = rand() % 20;
		}
	
		//Ordinary memmove - from b to solution vector, the space for which is allocated in main:
		//memmove((void *)second->solution, (void *)a[i], sizeof(double) * first->dim);
		memmove_mpi((void *)second->solution, (void *)a[i], sizeof(double) * first->dim);
		printf("Solution: %f\n", second->solution[i]);	
	}
	
	matrix_free(a, first->size);
}


void *memmove_mpi(void* dest, const void* src, unsigned int n){	
	
	MPI_Status stat;
	int sendbuf, recvbuf;
	int rank = 0;
	int peer = (rank == 0) ? 1 : 0;
	char *pDest = (char *)dest;
	const char *pSrc =( const char*)src;
	unsigned char flag = 0; //Flag for copy requirement if overlap

	//Allocate memory for temp array
	char *tmp = (char *)malloc(sizeof(char) * n);
	if ((pSrc == NULL) && (pDest == NULL) && (NULL == tmp)){
		return NULL;
	}
	if((pSrc < pDest) && (pDest < pSrc + n)){
		for( pDest += n, pSrc += n; n--;){
			*--pDest = *--pSrc;
		}
	}
	else {
		while (n--){
			unsigned int a = 0;
			//copy src to tmp array
			for(a = 0;a < n; ++a){
				*(tmp + a) = *(pSrc+a);
				sendbuf = *(tmp + a);					
				//copy tmp to pDest using MPI_Sendrecv
				MPI_Sendrecv(&sendbuf, 0, MPI_INT, peer, MPI_ANY_TAG, &recvbuf, 1, MPI_INT, peer, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
				printf("MPI process %d RECEIVED value %d from MPI process %d \n", rank, recvbuf, peer);	
				(pDest + a) = (char)&recvbuf;
			}
			free(tmp); //Free allocated memory
		}
	}
	return dest;
}


void list_mpi(list_a_t *first, list_b_t *second){
	
	int i, j;
	double **b = matrix_new(first->size, first->dim);
	srand(time(NULL));
	second->error = DBL_MAX;

	//Parallel	
	for (i=0;i<first->size;i++){
		for(j=0;j<first->dim;j++){
			b[i][j] = rand() % 20;
		}
		 
		//printf("recvbuf: %f\n", second->solution[i]);				
	}
	matrix_free(b, first->size);
}


