#!/bin/bash

#SBATCH -N 4
#SBATCH -n 128
#SBATCH -p node
#SBATCH --exclusive

mpirun ./dp
