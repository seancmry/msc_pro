#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <mpi.h>
 
/**
 *  * @brief Illustrates how to create a communicator representing a 2D torus
 *   * topology.
 *    **/
int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
 
        // Size of the default communicator
        int size;
        MPI_Comm_size(MPI_COMM_WORLD, &size);
	int dims[2] = {0, 0};
        MPI_Dims_create(size, 2, dims);
 
        // Make both dimensions periodic
        int periods[2] = {0, 0};
         
        // Let MPI assign arbitrary ranks if it deems it necessary
    	int reorder = true;
          
        // Create a communicator given the 2D torus topology.
        MPI_Comm new_communicator;
        MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &new_communicator);
                           
	//Declare neighbours
	enum DIRECTIONS {DOWN,UP,LEFT,RIGHT};
	char* neighbours_names[4] = {"down", "up", "left", "right"};
	int neighbours_ranks[4];
    
	//Let consider dims[0] = X, so the shift tells us our left and right neighbours
	MPI_Cart_shift(new_communicator, 0, 1, &neighbours_ranks[LEFT], &neighbours_ranks[RIGHT]);
	      
	// Let consider dims[1] = Y, so the shift tells us our up and down neighbours
	MPI_Cart_shift(new_communicator, 1, 1, &neighbours_ranks[DOWN], &neighbours_ranks[UP]);

        // My rank in the new communicator
        int my_rank;
        MPI_Comm_rank(new_communicator, &my_rank);
	
	//Get my coords in the new communicator
	int my_coords[2];
    	MPI_Cart_coords(new_communicator, my_rank, 2, my_coords);
 
	for(int i = 0; i < 4; i++)
    	{
        if(neighbours_ranks[i] == MPI_PROC_NULL){
            printf("[MPI process %d] I have no %s neighbour.\n", my_rank, neighbours_names[i]);
        }else{
            printf("[MPI process %d] I have a %s neighbour: process %d.\n", my_rank, neighbours_names[i], neighbours_ranks[i]);
		}
    	}

    	// Print my location in the 2D torus.
   	printf("[MPI process %d] I am located at (%d, %d).\n", my_rank, my_coords[0],my_coords[1]);
        
    	MPI_Finalize();
             
	return EXIT_SUCCESS;
}
