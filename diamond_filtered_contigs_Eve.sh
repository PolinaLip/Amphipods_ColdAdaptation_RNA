#!/bin/bash

#SBATCH --job-name=diamond_filtered_contigs_Eve
#SBATCH --time=3-20:00:00
#SBATCH --nodes=1
#SBATCH --nodelist=k72
#SBATCH --ntasks=10
#SBATCH --mem=100G
#SBATCH --output=diamond_filtered_contigs_Eve_job_output.txt

cd /scr/k70san2/lpolina/Unaligned/trimmed_and_cleaned_SSv2/clean/clean2/corrected/assembly_filtering/filtered_contig_annotation/

diamond blastx \
        -d /scr/k70san2/lpolina/databases/nr_taxon.dmnd \
        -q ../Eve_assembly.filtered.fasta \
        -o Eve_filtered_contigs.diamond.csv \
        --sensitive \
        --max-target-seqs 1 \
        --threads 10 \
        -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms skingdoms sphylums stitle salltitles
