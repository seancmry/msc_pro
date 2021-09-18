#!/bin/bash

function quit() 
{
	exit
}

for i in $(seq 1 10); do
	for j in $(seq 25); do 
		./prog -m sample_map_OpenRooms.txt -a 100 -b 100 -n $i | tee -a $i.dat
	done
done

quit

