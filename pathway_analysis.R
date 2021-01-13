#BiocManager::install("pathview")
library(pathview)
library(dplyr)
data(gse16873.d) # test data

### Example from pathview package
pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id = "04110",
                   species = "hsa", out.suffix = "gse16873")

### Prepare data for analysis 
pathway_dir <- 'labeglo2/Transcriptomics/annotations/pathway_analisys/Gla/Gla_KEGG_annotation/kegg_annot_assembly/'
pathways_to_analyse <- read.table(file.path(pathway_dir, 'pathway_to_analyse2.txt'), # I manually chose the pathways I want to check
                                sep = '\t', header = F, colClasses = "character")
pathways_to_analyse <- subset(pathways_to_analyse, V1 != '04723') # for Gla

de_results <- res # res from DE_analysis.R
lfc_all <- ifelse(de_results$padj > 0.05 | is.na(de_results$padj), 
                  0, de_results$log2FoldChange)
names(lfc_all) <- row.names(res)

for (pathway in pathways_to_analyse$V1){
  pv.out <- pathview(gene.data = lfc_all, pathway.id = pathway, gene.idtype = 'kegg',
                     species = "ko", out.suffix = "gla_all_padj0.05")
}

pv.out <- pathview(gene.data = lfc_all, pathway.id = "04723", gene.idtype = 'kegg',
                   species = "ko", out.suffix = "test")

pv.out <- pathview(gene.data = lfc_all, pathway.id = "03008", gene.idtype = 'kegg',
                   species = "ko", out.suffix = "test6", limit = list(gene=2),
                   bins = list(gene=30))

### Pathway analysis for all three species together
lfc_all_gla <- lfc_all
lfc_all_ecy <- lfc_all
lfc_all_eve <- lfc_all

#save(lfc_all_gla, lfc_all_ecy, lfc_all_eve, file = 'lfc_all_species.RData')
load('labeglo2/Transcriptomics/annotations/pathway_analisys/lfc_all_species.RData')

lfc_all_species <- full_join(data.frame(x = lfc_all_gla, y = names(lfc_all_gla))[-1,], 
                             data.frame(x = lfc_all_ecy, y = names(lfc_all_ecy)), 
                             by = 'y')
lfc_all_species <- full_join(lfc_all_species, 
                             data.frame(x = lfc_all_eve, y = names(lfc_all_eve)),
                             by = 'y')

colnames(lfc_all_species) <- c('Gla', 'Kegg', 'Ecy', 'Eve')
row.names(lfc_all_species) <- lfc_all_species$Kegg
lfc_all_species <- lfc_all_species[,-2]

pv.out <- pathview(gene.data = lfc_all_species, 
                   pathway.id = "00472", gene.idtype = 'kegg',
                   species = "ko", out.suffix = "all_species", limit = list(gene=1),
                   bins = list(gene=20), multi.state = T, same.layer = T,
                   keys.align = "y")


