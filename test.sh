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

F griewank

quit
echo foo
