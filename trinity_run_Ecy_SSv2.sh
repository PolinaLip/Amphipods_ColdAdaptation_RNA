#!/bin/bash

#SBATCH --job-name=trinity_Ecy_assembly
#SBATCH --time=10-12:00:00                            # Runtime in D-HH:MM:SS
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=240G
#SBATCH --output=trinity_assembly_Ecy_job_output.txt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=polina.lipaeva@gmail.com

cd /scr/k70san2/lpolina/Unaligned/trimmed_and_cleaned_SSv2/clean/clean2/corrected

/scr/k70san2/lpolina/apps/trinityrnaseq-v2.9.1/Trinity --seqType fq \
        --left  129_SSv2_Index-46_S33_clean2.fastq.1.cor.fq.gz,131_SSv2_Index-17_S35_clean2.fastq.1.cor.fq.gz,130_SSv2_Index-3_S34_clean2.fastq.1.cor.fq.gz,133_SSv2_Index-13_S36_clean2.fastq.1.cor.fq.gz,560_SSv2_Index-23_S37_clean2.fastq.1.cor.fq.gz,562_SSv2_Index-6_S39_clean2.fastq.1.cor.fq.gz,561_SSv2_Index-20_S38_clean2.fastq.1.cor.fq.gz,563_SSv2_Index-30_S40_clean2.fastq.1.cor.fq.gz \
        --right 129_SSv2_Index-46_S33_clean2.fastq.2.cor.fq.gz,131_SSv2_Index-17_S35_clean2.fastq.2.cor.fq.gz,130_SSv2_Index-3_S34_clean2.fastq.2.cor.fq.gz,133_SSv2_Index-13_S36_clean2.fastq.2.cor.fq.gz,560_SSv2_Index-23_S37_clean2.fastq.2.cor.fq.gz,562_SSv2_Index-6_S39_clean2.fastq.2.cor.fq.gz,561_SSv2_Index-20_S38_clean2.fastq.2.cor.fq.gz,563_SSv2_Index-30_S40_clean2.fastq.2.cor.fq.gz \
        --SS_lib_type FR \
        --max_memory 216G \
        --CPU 20 \
        --output EcySSv2_trinity \
 
