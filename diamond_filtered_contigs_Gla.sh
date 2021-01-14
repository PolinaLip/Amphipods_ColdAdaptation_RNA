#!/bin/bash

#SBATCH --job-name=diamond_filtered_contigs_Gla
#SBATCH --time=3-20:00:00
#SBATCH --nodes=1
#SBATCH --nodelist=k74
#SBATCH --ntasks=10
#SBATCH --mem=100G
#SBATCH --output=diamond_filtered_contigs_Gla_job_output.txt

cd /scr/k70san2/lpolina/Unaligned/trimmed_and_cleaned_SSv2/clean/clean2/corrected/assembly_filtering/filtered_contig_annotation/

diamond blastx \
        -d /scr/k70san2/lpolina/databases/nr_taxon.dmnd \
        -q ../Gla_assembly.filtered.fasta \
        -o Gla_filtered_contigs.diamond.csv \
        --sensitive \
        --max-target-seqs 1 \
        --threads 10 \
        -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms skingdoms sphylums stitle salltitles

