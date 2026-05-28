#!/bin/bash
#----------------------------------------------------
# Example SLURM job script to run hybrid applications 
# (MPI/OpenMP or MPI/pthreads) on TACC's Stampede 
# system.
#----------------------------------------------------
#SBATCH -J mpi_job     # Job name
#SBATCH -o stdout.o%j # Name of stdout output file(%j expands to jobId) 
#SBATCH -e stderr.o%j # Name of stderr output file(%j expands to jobId)
#SBATCH -p normal         # Serial queue for serial and MPI jobs
#SBATCH -N 1              # Total number of nodes requested (16 cores/node)
#SBATCH -n 1              # Total number of mpi tasks requested
#SBATCH -t 2:00:00       # Run time (hh:mm:ss) - 1.5 hours
# The next line is required if the user has more than one project
#SBATCH -A TG-EAR160027  # <-- Allocation name to charge job against

# Set the number of threads per task(Default=1)
#export OMP_NUM_THREADS=16

# Run the MPI application
./command_comp.sh 380 380 50
./downsample.py 380
./txt2vtk `pwd` 380

