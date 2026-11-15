#!/bin/bash
# DIAMOND BLASTX
diamond blastx --threads 24 -d /mnt/program/db/VFDB/vfdb_reference.dmnd -q $1_1.fasta.gz -o $1_1_output.tsv --max-target-seqs 1

# Check if second file exists
if [ -f "$1_2.fasta.gz" ]; then
    # If it exists, run DIAMOND BLASTX on it
    diamond blastx --threads 24 -d /mnt/program/db/VFDB/vfdb_reference.dmnd -q $1_2.fasta.gz -o $1_2_output.tsv --max-target-seqs 1
    
    # Combine both output files
    cat $1_1_output.tsv $1_2_output.tsv > $1_output.tsv
    
    # Clean up intermediate files
    rm $1_1_output.tsv $1_2_output.tsv
else
    # If second file doesn't exist, just rename the first output
    mv $1_1_output.tsv $1_output.tsv
fi
