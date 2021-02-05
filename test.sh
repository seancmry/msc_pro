#!/bin/bash

function quit() 
{
	exit
}

function F()
{
	./prog $1
}

F griewank
F sphere

quit
echo foo
