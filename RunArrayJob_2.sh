#!/bin/bash
#SBATCH --job-name=Net_array_job
#SBATCH --array=1-98
#SBATCH --output=net_array_job_%A-%a.out
#SBATCH --partition normal
#SBATCH --ntasks 1
#SBATCH --cpus-per-task=10       # <-- use 10 cores
#SBATCH --mem-per-cpu=1024M
#SBATCH --time=0-48:00:00

n=$SLURM_ARRAY_TASK_ID
Dir="$( sed "${n}q;d" InputFiles.txt )"
FILES=$Dir"/Network*.csv"

echo "Running task $SLURM_ARRAY_TASK_ID on directory $Dir"
echo "Using $SLURM_CPUS_PER_TASK CPUs per task"

# Collect all realizations in this directory
REALIZATIONS=($FILES)
NUM_REAL=${#REALIZATIONS[@]}
BATCH_SIZE=$SLURM_CPUS_PER_TASK

# Loop over realizations in batches
for ((i=0; i<NUM_REAL; i+=BATCH_SIZE)); do
    for ((j=0; j<BATCH_SIZE && i+j<NUM_REAL; j++)); do
        f=${REALIZATIONS[i+j]}
        NetPath="$f"
        NetLabel=$(grep -o '_Net[0-9]*r[0-9]*' <<< "$f")
        FileLabel="2Strain${NetLabel}"

        echo "Running $FileLabel"

        ConfigPath="$PWD/ModelConfig_2Strains.txt"
        # Run one realization in background using 1 core
        srun -n1 -c1 MultiStrainSIRonNet.exe $ConfigPath $NetPath $FileLabel >> log_${FileLabel}.out &
    done
    wait  # wait for this batch to finish before next
done


