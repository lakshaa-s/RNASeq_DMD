#!/bin/bash
# Navigate to the directory where you want to save FASTQ files
cd /data/scratch/bt22615/fastq_files

# Create output directory if it doesn't exist
mkdir -p raw_fastq

# List of SRA accession numbers (space-separated)
SRA_LIST="SRR12798466 SRR12798469 SRR12798472 SRR12798475 SRR12798478 SRR12798481 SRR16969732 SRR16969737 SRR16969742 SRR16969747 SRR16969752 SRR16969757 SRR16969762 SRR16969767"

# Load SRA Toolkit (if required on HPC)
module load sra-toolkit/2.11.0

# Loop through each accession and download FASTQ files
for SRA in $SRA_LIST
do
    echo "Downloading $SRA..."
    fastq-dump --split-files --gzip --readids --read-filter pass --outdir raw_fastq $SRA
done

echo "FASTQ download completed. Files saved in 'raw_fastq' directory."

