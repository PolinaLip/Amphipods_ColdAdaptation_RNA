# Amphipods_ColdAdaptation_RNA

Here are scripts which were used during work on transcriptomes analysis of Baikal amphipods (*Eulimnogammarus verrucosus*, *E. cyaneus*) and Holarctic species (*Gammarus lacustris*) exposed to average winter (1.5C) and average summer (12C) temperatures for one month. 

1. fix_NA_in_Annotatio.py: To fix problem with N/A in the Diamond annotation: remove the hypothetical proteins and recover taxonomy (taxid, superkingdom, kingdom and phylum) for the rest.
2. taxonomy_analysis.R: This script was used to extract only metazoan contigs for the further assemblies filtering (to avoid contigs coming from symbionts or commensales). The figure S2B is plotted by this script.
3. remove_contamination.py: This script allows to filter assemblies based on the output file from taxonomy_analysis.R.
4. Up_and_Down_plot.R: figure 2.
5. busco_stat.R: figure S2A.
6. DE_analysis.R: This script is for the differential expression analysis with DESeq2. Also: figure 3B, 3D, 4, 5A, 5B, 5D, 5F.
7. gsea_transcriptomics.R: To perform gsea analysis with fgsea package. Figures: 3A, 3C, 5C, 5E. gsea.R: initial scripts with different plots tryouts.
8. kegg.R: This script is for the Fisher’s exact test of the kegg pathway enrichment.
9. pathway_analysis.R: figure S4, S5, S6.  
