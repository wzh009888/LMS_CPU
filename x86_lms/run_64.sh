#!/bin/bash

#SBATCH -N 2
#SBATCH -n 64
#SBATCH -p node
#SBATCH --exclusive

mpirun ./demo2
