#!/bin/sh
#SBATCH -n 1 # 1 CPU core
#SBATCH -t 0-00:01:00 # 1 min
#SBATCH -p compute # partition name
#SBATCH -J test-job # sensible name for the job

# load up the correct modules
module load gcc/9.3.0
module load gsl/2.5 
module load openmpi/3.1.6

# launch
mpirun -n 1 ./prog -m sample_map_OpenRooms.txt -a 100 -b 100 -n 2
