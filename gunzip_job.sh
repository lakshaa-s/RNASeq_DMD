#!/bin/bash
#$ -cwd           # Set the working directory for the job to the current directory
#$ -pe smp 1      # Request 12 core
#$ -l h_rt=72:0:0  # Request 72 hour runtime
#$ -l h_vmem=96G   # Request 1GB RAM

# Load any necessary modules (if required)
module load gzip

# Decompress all .gz files in the current directory
gunzip *.gz
