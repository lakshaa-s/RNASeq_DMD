#!/bin/bash
#$ -cwd           # Set the working directory for the job to the current directory
#$ -pe smp 1      # Request 12 core
#$ -l h_rt=72:0:0  # Request 72 hour runtime
#$ -l h_vmem=96G   # Request 1GB RAM

# Load required modules
module load tophat


tophat2 -p 8 -o "$2" -G \
/data/home/bt22615/hg38/genome.gtf \
/data/home/bt22615/hg38/bowtie2/genome \
"trimmed_files/$1_1_val_1.fq" "trimmed_files/$1_2_val_2.fq"

# Index the BAM file
samtools index "$2"

echo "Alignment complete!"
