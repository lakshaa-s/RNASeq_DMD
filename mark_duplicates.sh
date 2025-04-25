#!/bin/bash
#$ -cwd           # Set the working directory for the job to the current directory
#$ -pe smp 1      # Request 12 core
#$ -l h_rt=72:0:0  # Request 72 hour runtime
#$ -l h_vmem=96G   # Request 1GB RAM

# Load required modules
module load samtools
module load java

# Ensure bam directory exists
mkdir -p bam

# Add read group information to the BAM file
samtools addreplacerg -r "@RG\tID:RG1\tSM:SampleName\tPL:Illumina\tLB:Library.fa" -o bam/${1}_RG_sorted.bam bam/${1}_sorted.bam

# Run Picard MarkDuplicates with the new command and filenames
java -Xmx32g -jar /data/scratch/bt22615/fastq_files/picard.jar MarkDuplicates \
    I="bam/${1}_RG_sorted.bam" \
    O="bam/${1}_RG_sorted_deduplicated.bam" \
    M="bam/${1}_RG_sorted_deduplicated_metrics.txt" \
    VALIDATION_STRINGENCY=LENIENT \
    REMOVE_DUPLICATES=true

echo "Deduplication completed for ${1}."
