#!/bin/bash

#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p node
#SBATCH --exclusive

mpirun ./demo2
