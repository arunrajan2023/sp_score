#!/bin/bash
#SBATCH --job-name="py"
#SBATCH --output="envinfo.%j.%N.out"
#SBATCH --partition=large-shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --mem=400G	   #400GB requested for large-shared queue
#SBATCH --export=ALL
#SBATCH --account=ucn104
#SBATCH -t 01:30:00

module purge
module load slurm
module load cpu
module load gpu/0.15.4  intel/19.0.5.281  intel-mpi/2019.8.254

export PYTHONPATH=$PYTHONPATH:/home/arajan63/work/pscore/sp_score
export PYTHONUNBUFFERED=TRUE

#1 node 64 cores
/home/arajan63/miniconda3/bin/python tests/test_pscore.py
