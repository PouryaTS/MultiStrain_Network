#!/bin/bash
#SBATCH --job-name=Net_array_job
#SBATCH --array=1-14
#SBATCH --output=net_array_job_%A-%a.out
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=1024M
#SBATCH --time=0-48:00:00

# Get directory for this array task
n=$SLURM_ARRAY_TASK_ID
Dir="$(sed "${n}q;d" InputFiles.txt)"

echo "Running task $SLURM_ARRAY_TASK_ID on directory $Dir"
echo "Using $SLURM_CPUS_PER_TASK CPUs"

# Collect all CSV files
REALIZATIONS=("$Dir"/Network*.csv)

# Config file
ConfigPath="$PWD/ModelConfig_2strains.txt"

# Parallel control
MAX_JOBS=$SLURM_CPUS_PER_TASK
running=0

# Loop over all files
for f in "${REALIZATIONS[@]}"; do

    # Extract labels
    NetLabel=$(grep -o '_Net[0-9]*r[0-9]*' <<< "$f")
    FileLabel="2Strain${NetLabel}"

    echo "Running $FileLabel"

    # Run simulation in background
    srun -n1 -c1 MultiStrainSIRonNet.exe $ConfigPath $NetPath $FileLabel >> log_${FileLabel}.out &
    ((running++))

    # If max parallel jobs reached → wait for one to finish
    if (( running >= MAX_JOBS )); then
        wait -n
        ((running--))
    fi

done

# Wait for all remaining jobs
wait

echo "All realizations finished"