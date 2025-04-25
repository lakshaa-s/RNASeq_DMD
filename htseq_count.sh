#!/bin/bash
#$ -cwd           # Set the working directory for the job to the current directory
#$ -pe smp 1      # Request 12 core
#$ -l h_rt=72:0:0  # Request 72 hour runtime
#$ -l h_vmem=96G   # Request 1GB RAM

# Load Miniconda module
module load miniconda

# Activate the HTSeq environment
conda activate htseq

# Ensure the output directory exists
mkdir -p htseq

# Run HTSeq-count
htseq-count -f bam bam/${1}_RG_sorted_deduplicated.bam \
    /data/home/bt22615/hg38/genome.gtf \
    > htseq/${1}.txt

# Deactivate the conda environment
conda deactivate

echo "HTSeq-count completed for ${1}."
