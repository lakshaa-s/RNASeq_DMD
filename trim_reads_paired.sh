#!/bin/bash
#$ -cwd           # Set the working directory for the job to the current directory
#$ -pe smp 1      # Request 12 core
#$ -l h_rt=72:0:0  # Request 72 hour runtime
#$ -l h_vmem=96G   # Request 1GB RAM

module load trimgalore

# Output directory
OUTPUT_DIR="trimmed_files"
mkdir -p $OUTPUT_DIR

# Run Trim Galore
trim_galore --paired --fastqc -o $OUTPUT_DIR "$1_1.fastq" "$1_2.fastq" 

echo "Trimming complete. Files saved in $OUTPUT_DIR"
