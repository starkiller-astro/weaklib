#!/bin/bash
## attempt to build weaklib and reproduce weaklib compile error with CCE on spock
## should be run from $MEMBERWORK or other GPFS directory (weaklib-tables is large)

## setup the compile/run environment
export WEAKLIB_MACHINE=spock_cce
module load PrgEnv-cray craype-accel-amd-gfx908 rocm cray-hdf5-parallel

## compile
make OPT_LEVEL=OPTIMIZE USE_GPU=TRUE USE_HIP=TRUE USE_OMP_OL=TRUE wlOpacityPerformanceTest