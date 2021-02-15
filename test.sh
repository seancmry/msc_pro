#!/bin/bash

function quit() 
{
	exit
}

if [ -f output.dat ]
then
	rm output.dat
fi

function F()
{
	./prog $1 | tee -a output.dat
}

for i in $(seq 10)
do
	#F griewank
	#F ackley
	#F sphere
	F rosenbrock
done

quit
echo foo
