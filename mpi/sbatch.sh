#!/bin/bash

#SBATCH -n 16
#SBATCH -t 0-0:1:00 # 1 minutes
#SBATCH -p compute # partition name
#SBATCH -J test-job # sensible name for the job

# load modules
module load gcc/9.3.0 gsl/2.5 openmpi/3.1.6

# launch the code
mpirun -n 16 ./prog ackley

