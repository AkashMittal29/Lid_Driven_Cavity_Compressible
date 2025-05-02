#!/bin/bash

#SBATCH -J 'Lid_driven'
#SBATCH -A mecfd_q
#SBATCH -t 00:05:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --cpus-per-task=4

module purge
module load intel openmpi

export OMP_NUM_THREADS=4
srun run.exe >output.log 2>error.log

echo "FINISHED"