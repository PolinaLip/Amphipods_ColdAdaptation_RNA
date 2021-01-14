#!/bin/bash

#SBATCH --job-name=trinity_Gl_assembly
#SBATCH --time=10-12:00:00                            # Runtime in D-HH:MM:SS
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=240G
#SBATCH --output=trinity_assembly_Gl_job_output.txt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=polina.lipaeva@gmail.com

cd /scr/k70san2/lpolina/Unaligned/trimmed_and_cleaned_SSv2/clean/clean2/corrected

/scr/k70san2/lpolina/apps/trinityrnaseq-v2.9.1/Trinity --seqType fq \
        --left 311_SSv2_Index-10_S41_clean2.fastq.1.cor.fq.gz,313_SSv2_Index-15_S43_clean2.fastq.1.cor.fq.gz,312_SSv2_Index-29_S42_clean2.fastq.1.cor.fq.gz,315_SSv2_Index-11_S44_clean2.fastq.1.cor.fq.gz,594_SSv2_Index-12_S45_clean2.fastq.1.cor.fq.gz,596_SSv2_Index-5_S47_clean2.fastq.1.cor.fq.gz,595_SSv2_Index-7_S46_clean2.fastq.1.cor.fq.gz,598_SSv2_Index-4_S48_clean2.fastq.1.cor.fq.gz \
        --right 311_SSv2_Index-10_S41_clean2.fastq.2.cor.fq.gz,313_SSv2_Index-15_S43_clean2.fastq.2.cor.fq.gz,312_SSv2_Index-29_S42_clean2.fastq.2.cor.fq.gz,315_SSv2_Index-11_S44_clean2.fastq.2.cor.fq.gz,594_SSv2_Index-12_S45_clean2.fastq.2.cor.fq.gz,596_SSv2_Index-5_S47_clean2.fastq.2.cor.fq.gz,595_SSv2_Index-7_S46_clean2.fastq.2.cor.fq.gz,598_SSv2_Index-4_S48_clean2.fastq.2.cor.fq.gz \
        --SS_lib_type FR \
        --max_memory 216G \
        --CPU 20 \
        --output GlSSv2_trinity \

