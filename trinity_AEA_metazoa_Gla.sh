#!/bin/bash

#SBATCH --job-name=trinity_AEA_metazoa_Gla
#SBATCH --time=8:00:00
#SBATCH --nodelist=k74
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem=100G
#SBATCH --output=trinity_AEA_metazoa_Gla_job_output.txt

cd /scr/k70san2/lpolina/Unaligned/trimmed_and_cleaned_SSv2/clean/clean2/corrected/

/scr/k70san2/lpolina/apps/trinityrnaseq-v2.9.1/util/align_and_estimate_abundance.pl \
        --seqType fq \
        --transcripts assembly_filtering/Gla_assembly_metazoa.fasta \
        --samples_file assembly_filtering/Gla_samples.txt \
        --SS_lib_type FR \
        --est_method salmon \
        --output_dir assembly_filtering/AEA_metazoa_output_Gla \
        --thread_count 10 \
        --gene_trans_map GlSSv2_trinity/Trinity.fasta.gene_trans_map \
        --prep_reference

