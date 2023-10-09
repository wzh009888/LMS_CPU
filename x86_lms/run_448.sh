#!/bin/bash

#SBATCH -N 8
#SBATCH -n 448
#SBATCH -p node
#SBATCH --exclusive

mpirun ./demo2

