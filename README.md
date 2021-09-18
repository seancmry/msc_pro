# M.Sc. Project

Modules loaded for the running of the demo and serial path components prior to the cluster system updates in March 2021:
```
module load cports
module load gsl/2.2.1-gnu
```
Scripts now updated when loading the new version through the following:
```
module load gcc/9.3.0
module load gsl/2.5
module load openmpi/3.1.6
```
For the compilation, run
```
make
```
and for the execution of the demo functions

```
mpirun -n 4 ./prog [function name]
```
The list of available functions is given below

```
ackley
griewank
rosenbrock
sphere
```
For the purposes of running instrumented code, it is necessary to load the SCOREP module:
```
module load scorep/6.0
```
![alt text][img]

[img]: https://github.com/seancmry/msc_pro/blob/master/report/plots/screenoutput.PNG

