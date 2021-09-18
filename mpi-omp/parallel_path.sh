#!/bin/bash



function quit() 
{
	exit
}



for j in $(seq 25); do 
	mpirun -n 8 ./prog -m sample_map_OpenRooms.txt -a 100 -b 100 | tee -a n8.dat
done


quit


