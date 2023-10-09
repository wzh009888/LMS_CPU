#!/bin/bash

export OMP_NUM_THREADS=32
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH -p node
#SBATCH --exclusive

mpirun ./demo2
#./demo2
