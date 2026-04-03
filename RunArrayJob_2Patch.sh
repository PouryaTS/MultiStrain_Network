#!/bin/bash
#SBATCH --job-name=Net_array_job
#SBATCH --array=1-14
#SBATCH --output=net_array_job_%A-%a.out
#SBATCH --partition normal
#SBATCH --ntasks 1
#SBATCH --cpus-per-task=10       # <-- use 10 cores
#SBATCH --mem-per-cpu=1024M
#SBATCH --time=0-48:00:00

n=$SLURM_ARRAY_TASK_ID
Dir="$( sed "${n}q;d" InputFiles_MDSD.txt )"
FILES=$Dir"/Network*.csv"

echo "Running task $SLURM_ARRAY_TASK_ID on directory $Dir"
echo "Using $SLURM_CPUS_PER_TASK CPUs per task"

Set=$(basename "$Dir")

# Collect all realizations in this directory
REALIZATIONS=($FILES)
NUM_REAL=${#REALIZATIONS[@]}
BATCH_SIZE=$SLURM_CPUS_PER_TASK

# Loop over realizations in batches
for ((i=0; i<NUM_REAL; i+=BATCH_SIZE)); do
    for ((j=0; j<BATCH_SIZE && i+j<NUM_REAL; j++)); do
        f=${REALIZATIONS[i+j]}
        r=$(echo "$f" | sed -n 's/.*_Net[0-9]\+r\([0-9]\+\)_.*/\1/p')
        NetLabel1=$(grep -o '_Net[0-9]*r[0-9]*' <<< "$f")
        Net1="$f"

        # Build Net2 path (Net6)
        
        # Reverse index for Net2
        r2=$((39 - r))
        # Determine correct Set for Net2
        # Determine correct Set for Net2
        if [ "$r2" -le 19 ]; then
            Set2="Set1"
        else
            Set2="Set2"
        fi

	    pathnet2="/home/toranj/output/SHNB_Net/NetRealization/SHNBnet10/Net44/${Set2}"
        Net2="${pathnet2}/NetworkEdgelist_BothLayers_SHNB_Net44r${r2}_lambda=2.0_NBmean=4.0_NBstd=4.0_N=10000.csv"
        NetLabel2="_Net44r${r2}"
        FileLabel="2Strain${NetLabel1}${NetLabel2}"

        echo "Running $FileLabel : Net1=$Net1 , Net2=$Net2"

        # Each realization has its own input file
        netinputfile="network_input_${SLURM_JOB_ID}_${r}.txt"
        echo "$Net1" > $netinputfile
        echo "$Net2" >> $netinputfile

        ConfigPath="$PWD/ModelConfig_2Strains_2Patch.txt"

        # Run one realization in background using 1 core
        srun -n1 -c1 MultiStrainSIRonNet_2patches.exe $ConfigPath $netinputfile $FileLabel >> log_${FileLabel}.out &
    done
    wait  # wait for this batch to finish before next
done


