#!/bin/bash

#SBATCH --job-name=trinity_Eve_assembly
#SBATCH --time=10-12:00:00                            # Runtime in D-HH:MM:SS
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=240G
#SBATCH --output=trinity_assembly_Eve_job_output.txt

cd /scr/k70san2/lpolina/Unaligned/trimmed_and_cleaned_SSv2/clean/clean2/corrected

/scr/k70san2/lpolina/apps/trinityrnaseq-v2.9.1/Trinity --seqType fq \
	--left 350_SSv2_Index-14_S25_clean2.fastq.1.cor.fq.gz,632_SSv2_Index-31_S29_clean2.fastq.1.cor.fq.gz,351_SSv2_Index-1_S26_clean2.fastq.1.cor.fq.gz,633_SSv2_Index-19_S30_clean2.fastq.1.cor.fq.gz,352_SSv2_Index-22_S27_clean2.fastq.1.cor.fq.gz,634_SSv2_Index-16_S31_clean2.fastq.1.cor.fq.gz,354_SSv2_Index-8_S28_clean2.fastq.1.cor.fq.gz,636_SSv2_Index-9_S32_clean2.fastq.1.cor.fq.gz\
	--right 350_SSv2_Index-14_S25_clean2.fastq.2.cor.fq.gz,632_SSv2_Index-31_S29_clean2.fastq.2.cor.fq.gz,351_SSv2_Index-1_S26_clean2.fastq.2.cor.fq.gz,633_SSv2_Index-19_S30_clean2.fastq.2.cor.fq.gz,352_SSv2_Index-22_S27_clean2.fastq.2.cor.fq.gz,634_SSv2_Index-16_S31_clean2.fastq.2.cor.fq.gz,354_SSv2_Index-8_S28_clean2.fastq.2.cor.fq.gz,636_SSv2_Index-9_S32_clean2.fastq.2.cor.fq.gz\
	--SS_lib_type FR \
	--max_memory 216G \
	--CPU 20 \
	--output EveSSv2_trinity 

