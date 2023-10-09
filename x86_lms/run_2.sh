#!/bin/bash

#SBATCH -N 1
#SBATCH -n 2
#SBATCH -p node
#SBATCH --exclusive

mpirun ./demo2
