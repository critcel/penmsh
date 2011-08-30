#! /bin/bash


## set ifort to the correct path
module load intel_compilers/12.0.3.174 intel_mpi/4.0.1.007

name=pentranp.941s-EDF
rm *.mod
mpiifort  -O2 -mcmodel=medium ${name}.f  -o ${name}.o2
rm *.mod
mpiifort  -O3 -mcmodel=medium ${name}.f  -o ${name}.o3

##ifort -extend-source 72  -O2 ${name}.f -o ${name}.o2
##ifort -extend-source 72  -O3 ${name}.f -o ${name}.o3

