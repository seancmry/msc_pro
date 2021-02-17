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

for i in $(seq 25)
do
	F ackley
done

for i in $(seq 25)
do 
	F griewank
done

for i in $(seq 25)
do
	F rosenbrock
done

for i in $(seq 25)
do	
	F sphere
done

quit
echo foo
