#!/bin/bash

function quit() 
{
	exit
}

function F()
{
	mpirun -n 9 ./prog $1
}

for i in $(seq 25)
do
	F ackley | tee -a ackley.dat
done

for i in $(seq 25)
do 
	F griewank | tee -a griewank.dat
done

for i in $(seq 25)
do
	F rosenbrock | tee -a rosenbrock.dat
done

for i in $(seq 25)
do	
	F sphere | tee -a sphere.dat
done

quit

