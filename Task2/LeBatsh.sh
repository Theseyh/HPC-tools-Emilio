#!/bin/bash
#SBATCH -o %x-%J.out                  # Output file
#SBATCH -e %x-%J.error                # Error file
#SBATCH --time=0-00:05:00             # Requested time to run the job
#SBATCH -n 64                          # Total number of processes
#SBATCH --mem-per-cpu=1G              # Memory per CPU
#SBATCH --exclusive                   # Request exclusive access to the node


module load intel impi

# Run the a.out executable 5 times with 4 MPI processes
for i in {1..5}; do
    echo "Run #$i"
    ./spmv
done
