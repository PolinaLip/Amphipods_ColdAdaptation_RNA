library(fgsea)
library(Biobase)
library(ggplot2)
library(ggforce)
library(ggsci)
library(shadowtext)
library(viridis)

### Make paths to required files:
dir <- 'labeglo2/Transcriptomics/quantification/AEA_metazoa_Eve/'
dir_to_write <- 'labeglo2/Transcriptomics/quantification/AEA_metazoa_Eve/figures'
DE_sample <- 'condition_1.5C_vs_12C_0.05_FA.csv'
eggnog_FA <- 'eve_metazoa_eggnog2.csv'
eggnog_FA <- 'ecy_metazoa_eggnog3.csv'
eggnog_FA <- 'gla_metazoa_eggnog2.csv'
DE_sample <- 'Eve_FA_condition_1.5C_vs_12C.csv'

DE_gene_list <- read.csv(file.path(dir, DE_sample), sep = '\t')
go_terms <- read.csv(file.path(dir, eggnog_FA), sep = '\t', header = F)
go_info <- read.csv('labeglo2/Transcriptomics/annotations/functional_annotation/GO_info.csv', sep = '\t', header = T) # tsv file with three columns: GO id, GO description, and root (BC, MF, or CC)

### Files preprocessing:
go_terms <- subset(go_terms, V6 != '') # remove rows which do not contain gene name
GOs <- strsplit(as.character(go_terms$V7), ',') # extract GO terms (this is gene2GO list)
names(GOs) <- go_terms$V6 # assign gene names to it
DE_gene_stat <- DE_gene_list$stat # transcriptomics
names(DE_gene_stat) <- DE_gene_list$X
names(DE_gene_stat) <- DE_gene_list$gene
### Make GO2gene annotation list and put description of GO terms to GO2gene list:
GOs2gene <- reverseSplit(GOs) # make GO2gene from gene2Go
GOs2gene <- lapply(GOs2gene, unique) # remove duplicates of genes in GO2gene
tmp <- match(names(GOs2gene), go_info$GO_ID) # to select only required GOs from GO_info
names(GOs2gene) <- go_info[tmp,]$Description

### Split data on three categories:
GOs2gene_BP <- GOs2gene[go_info[tmp,]$Root == 'biological_process']
GOs2gene_MF <- GOs2gene[go_info[tmp,]$Root == 'molecular_function']
GOs2gene_CC <- GOs2gene[go_info[tmp,]$Root == 'cellular_component']

### Run GSEA
set.seed(42)
fgseaRes_BP <- fgsea(GOs2gene_BP, DE_gene_stat)
fgseaRes_MF <- fgsea(GOs2gene_MF, DE_gene_stat)
fgseaRes_CC <- fgsea(GOs2gene_CC, DE_gene_stat)
save(fgseaRes_BP, file = file.path(dir, 'gsea_BP_Eve1.5C_DE0.05.RData'))
save(GOs2gene_MF, file = file.path(dir, 'gsea_MF_EcyGla1.5C_DE0.01lfc2.RData'))
load('labeglo2/Transcriptomics/quantification/AEA_metazoa_Gla/gsea_BP_1.5vs12C.RData')
load('labeglo2/Transcriptomics/quantification/AEA_metazoa_Eve/gsea_BP_Eve1.5C_DE0.05.RData')

fgseaRes_sign_BP <- cbind(subset(fgseaRes_BP, pval < 0.05), 
                          rep('Biological process', 
                              nrow(subset(fgseaRes_BP, pval < 0.05))))
fgseaRes_sign_BP <- subset(fgseaRes_sign_BP, size > 2)
fgseaRes_sign_MF <- cbind(subset(fgseaRes_MF, pval < 0.05), 
                          rep('Molecular function', 
                              nrow(subset(fgseaRes_MF, pval < 0.05))))
fgseaRes_sign_MF <- subset(fgseaRes_sign_MF, size > 5)
fgseaRes_sign_CC <- cbind(subset(fgseaRes_CC, pval < 0.05), 
                          rep('Cellular compartment', 
                              nrow(subset(fgseaRes_CC, pval < 0.05))))
fgseaRes_sign_CC <- subset(fgseaRes_sign_CC, size > 2)

### Illustration of GSEA

fgseaRes_sign_BP$pathway <- reorder(fgseaRes_sign_BP$pathway, fgseaRes_sign_BP$NES)
fgseaRes_sign_BP$gene_ratio <- fgseaRes_sign_BP$size / length(DE_gene_stat) * 100

fgseaRes_sign_BP_up <- subset(fgseaRes_sign_BP, NES > 0)
fgseaRes_sign_BP_down <- subset(fgseaRes_sign_BP, NES < 0)

fgseaRes_sign_BP_up_merged <- fgseaRes_sign_BP_up
#---------------G.lacusrtis--BP-merging--UP----------------
fgseaRes_sign_BP_up_merged <- fgseaRes_sign_BP_up_merged[fgseaRes_sign_BP_up_merged$pathway != 'GO:0000398_mRNA splicing, via spliceosome']
fgseaRes_sign_BP_up_merged <- fgseaRes_sign_BP_up_merged[!grepl('GO:0000377', fgseaRes_sign_BP_up_merged$pathway)]
fgseaRes_sign_BP_up_merged <- fgseaRes_sign_BP_up_merged[!grepl('GO:0000375', fgseaRes_sign_BP_up_merged$pathway)]
fgseaRes_sign_BP_up_merged <- fgseaRes_sign_BP_up_merged[!grepl('GO:0022613|GO:0022618|GO:0006364|GO:0030490', fgseaRes_sign_BP_up_merged$pathway)] # ribosome biogenesis
fgseaRes_sign_BP_up_merged <- fgseaRes_sign_BP_up_merged[!grepl('GO:0051169', fgseaRes_sign_BP_up_merged$pathway)]
fgseaRes_sign_BP_up_merged <- fgseaRes_sign_BP_up_merged[!grepl('GO:0016072|GO:0006400|GO:0006399', fgseaRes_sign_BP_up_merged$pathway)] # ncRNA processing
fgseaRes_sign_BP_up_merged <- fgseaRes_sign_BP_up_merged[!grepl('GO:0016071', fgseaRes_sign_BP_up_merged$pathway)] # mRNA metabolic process to RNA metabolic process
fgseaRes_sign_BP_up_merged <- fgseaRes_sign_BP_up_merged[!grepl('GO:0034645', fgseaRes_sign_BP_up_merged$pathway)] # mRNA metabolic process to RNA metabolic process
fgseaRes_sign_BP_up_merged$pathway <- sub('.*_', '', fgseaRes_sign_BP_up_merged$pathway) # for figures
fgseaRes_sign_BP_up_merged$pathway <- reorder(fgseaRes_sign_BP_up_merged$pathway, fgseaRes_sign_BP_up_merged$NES)
#---------------E.verrucosus--BP-merging-UP--------------------------------------------
fgseaRes_sign_BP_up_merged <- fgseaRes_sign_BP_up_merged[!grepl('GO:0007338', fgseaRes_sign_BP_up_merged$pathway)] # fertilization
fgseaRes_sign_BP_up_merged <- fgseaRes_sign_BP_up_merged[!grepl('GO:0050000', fgseaRes_sign_BP_up_merged$pathway)] # establishment of chromosome localization
fgseaRes_sign_BP_up_merged <- fgseaRes_sign_BP_up_merged[!grepl('GO:0071826|GO:0022618', fgseaRes_sign_BP_up_merged$pathway)] # ribonucleoprotein complex biogenesis
fgseaRes_sign_BP_up_merged <- fgseaRes_sign_BP_up_merged[!grepl('GO:0090304|GO:0051259|GO:0034622|GO:0061077', fgseaRes_sign_BP_up_merged$pathway)]
fgseaRes_sign_BP_up_merged <- fgseaRes_sign_BP_up_merged[!grepl('GO:0009792', fgseaRes_sign_BP_up_merged$pathway)] # embryo development
fgseaRes_sign_BP_up_merged <- fgseaRes_sign_BP_up_merged[!grepl('GO:0065003', fgseaRes_sign_BP_up_merged$pathway)] # protein-containing complex subunit organization
#---------------E.cyaneus--BP-merging-UP-----------------------------------------
fgseaRes_sign_BP_up_merged <- fgseaRes_sign_BP_up_merged[!grepl('GO:0051716', fgseaRes_sign_BP_up_merged$pathway)] # response to stimulus
fgseaRes_sign_BP_up_merged <- fgseaRes_sign_BP_up_merged[!grepl('GO:0051603|GO:0044257', fgseaRes_sign_BP_up_merged$pathway)] # protein catabolic process
fgseaRes_sign_BP_up_merged <- fgseaRes_sign_BP_up_merged[!grepl('GO:0019941|GO:0006511|GO:0010498', fgseaRes_sign_BP_up_merged$pathway)] # modification-dependent macromolecule catabolic process
fgseaRes_sign_BP_up_merged <- fgseaRes_sign_BP_up_merged[!grepl('GO:0048878', fgseaRes_sign_BP_up_merged$pathway)] # homeostasis process
fgseaRes_sign_BP_up_merged <- fgseaRes_sign_BP_up_merged[!grepl('GO:0000377|GO:0000375|GO:0000398|GO:0006397', fgseaRes_sign_BP_up_merged$pathway)] # RNA splicing
fgseaRes_sign_BP_up_merged <- fgseaRes_sign_BP_up_merged[!grepl('GO:0003002', fgseaRes_sign_BP_up_merged$pathway)] # pattern specification process
fgseaRes_sign_BP_up_merged <- fgseaRes_sign_BP_up_merged[!grepl('GO:0044093|GO:0051347|GO:0033674', fgseaRes_sign_BP_up_merged$pathway)] # positive regulation of catalytic activity

fgseaRes_sign_BP_up_merged$pathway <- sapply(fgseaRes_sign_BP_up_merged$pathway, 
                                             function(s) paste(strwrap(s, width=45), 
                                                               collapse='\n'))
fgseaRes_sign_BP_up_merged$pathway <- factor(fgseaRes_sign_BP_up_merged$pathway)
fgseaRes_sign_BP_up_merged$pathway <- reorder(fgseaRes_sign_BP_up_merged$pathway, 
                                              fgseaRes_sign_BP_up_merged$gene_ratio)
fgseaRes_sign_BP_up_merged$pathway_short <- factor(sub('.*_', '', 
                                                fgseaRes_sign_BP_up_merged$pathway))
fgseaRes_sign_BP_up_merged$pathway_short <- reorder(fgseaRes_sign_BP_up_merged$pathway_short, 
                                              fgseaRes_sign_BP_up_merged$gene_ratio)

fgseaRes_sign_BP_up_merged <- subset(fgseaRes_sign_BP_up_merged, gene_ratio > 5)
ggplot(fgseaRes_sign_BP_up_merged, aes(gene_ratio, pathway_short)) +
  geom_point(aes(size = size, color = pval)) +
  ylab('') +
  xlab('Gene ratio, %') +
  scale_color_viridis('p-value', option = 'viridis', direction = -1, end = 0.6) + 
  scale_size_continuous('Number of genes', 
                        breaks = c(round(seq(min(fgseaRes_sign_BP_up_merged$size), 
                                             max(fgseaRes_sign_BP_up_merged$size), 
                                             length.out = 4)))) +
  theme_light() +
  theme(axis.text = element_text(size = 13), 
        axis.title.x = element_text(size = 13))

ggsave(filename = file.path(dir_to_write, 'dotplot_Eve_0.05_BP_up_2.png'),
       scale = 0.9, width = 8, height = 6.5)

fgseaRes_sign_BP_down_merged <- fgseaRes_sign_BP_down
#---
#-------------------G.lacusrtis--BP-merging--DOWN--------------------------------------------------
fgseaRes_sign_BP_down_merged <- fgseaRes_sign_BP_down_merged[
  !grepl('GO:1901293|GO:0009206|GO:0009260|GO:0009259|GO:0046390|GO:0046496|GO:0009124|GO:0009127|GO:0019693|GO:0006163|GO:0006164|GO:0009145|GO:0009142|GO:0009156|GO:0009152|GO:0009150|GO:0009161|GO:0009168|GO:0009165|GO:0006753|GO:0009201|GO:0019362', fgseaRes_sign_BP_down_merged$pathway)] # nucleotide metabolic process
fgseaRes_sign_BP_down_merged <- fgseaRes_sign_BP_down_merged[
  !grepl('GO:0002263|GO:0006955|GO:0002366|GO:0045321', fgseaRes_sign_BP_down_merged$pathway)] # immune system process
fgseaRes_sign_BP_down_merged <- fgseaRes_sign_BP_down_merged[
  !grepl('GO:0055074|GO:0072507|GO:0072503', fgseaRes_sign_BP_down_merged$pathway)] # cellular calcium ion homeostasis
fgseaRes_sign_BP_down_merged <- fgseaRes_sign_BP_down_merged[
  !grepl('GO:0036335', fgseaRes_sign_BP_down_merged$pathway)] # homeostasis of number of cells
fgseaRes_sign_BP_down_merged <- fgseaRes_sign_BP_down_merged[
  !grepl('GO:0016310', fgseaRes_sign_BP_down_merged$pathway)] # phosphate-containing compound metabolic process
fgseaRes_sign_BP_down_merged <- fgseaRes_sign_BP_down_merged[
  !grepl('GO:0072521|GO:0072524|GO:1901566|GO:0072522|GO:0006518', fgseaRes_sign_BP_down_merged$pathway)]  # organonitrogen compound metabolic process
fgseaRes_sign_BP_down_merged <- fgseaRes_sign_BP_down_merged[
  !grepl('GO:0055086|GO:0006796|GO:0006098|GO:0035336', fgseaRes_sign_BP_down_merged$pathway)] # nucleotide metabolic process
fgseaRes_sign_BP_down_merged <- fgseaRes_sign_BP_down_merged[
  !grepl('GO:0007018|GO:0097435|GO:0007010|GO:0070286|GO:0035082|GO:0043933|GO:0044085|GO:0070925|GO:0001578|GO:0044782|GO:0030031|GO:0034622|GO:0003341|GO:0007051|GO:0060271|GO:0065003|GO:0000226', 
         fgseaRes_sign_BP_down_merged$pathway)] # cellular component assembly
fgseaRes_sign_BP_down_merged <- fgseaRes_sign_BP_down_merged[
  !grepl('GO:0060285|GO:0019322|GO:0019321|GO:0033206|GO:0046843|GO:0006413|GO:0002181|GO:1901533', fgseaRes_sign_BP_down_merged$pathway)]
fgseaRes_sign_BP_down_merged <- fgseaRes_sign_BP_down_merged[
  !grepl('GO:0120036', fgseaRes_sign_BP_down_merged$pathway)]
#-----------------Eve--BP-merging--DOWN---------------
fgseaRes_sign_BP_down_merged <- fgseaRes_sign_BP_down_merged[
  !grepl('GO:0009583|GO:0009314', fgseaRes_sign_BP_down_merged$pathway)] # response to light stimulus
fgseaRes_sign_BP_down_merged <- fgseaRes_sign_BP_down_merged[
  !grepl('GO:1904950|GO:1903828|GO:0051224', fgseaRes_sign_BP_down_merged$pathway)] # negative regulation of transport
fgseaRes_sign_BP_down_merged <- fgseaRes_sign_BP_down_merged[
  !grepl('GO:0043067|GO:0010941|GO:0006915', fgseaRes_sign_BP_down_merged$pathway)] # regulation of apoptotic process
fgseaRes_sign_BP_down_merged <- fgseaRes_sign_BP_down_merged[
  !grepl('GO:0007165', fgseaRes_sign_BP_down_merged$pathway)] # signaling
fgseaRes_sign_BP_down_merged <- fgseaRes_sign_BP_down_merged[
  !grepl('GO:0003013', fgseaRes_sign_BP_down_merged$pathway)] # blood circulation
fgseaRes_sign_BP_down_merged <- fgseaRes_sign_BP_down_merged[
  !grepl('GO:0009581', fgseaRes_sign_BP_down_merged$pathway)] # detection of abiotic stimulus
fgseaRes_sign_BP_down_merged <- fgseaRes_sign_BP_down_merged[
  !grepl('GO:0043902|GO:0043436|GO:0006082|GO:0018209', fgseaRes_sign_BP_down_merged$pathway)]
#-------------------Ecy--BP-merging--DOWN-----------------------------------------------------------
fgseaRes_sign_BP_down_merged <- fgseaRes_sign_BP_down_merged[
  !grepl('GO:0032269|GO:0001933|GO:0042326|GO:0048523|GO:0010605|GO:0009890|GO:0045936|GO:0031324|GO:0031327|GO:0051248|GO:0031400|GO:0010558|GO:0010563|GO:0051172', 
         fgseaRes_sign_BP_down_merged$pathway)] # negative regulation of metabolic process
fgseaRes_sign_BP_down_merged <- fgseaRes_sign_BP_down_merged[
  !grepl('GO:0007268|GO:0098916|GO:0099537|GO:0099565|GO:0007267', fgseaRes_sign_BP_down_merged$pathway)] # synaptic signaling
fgseaRes_sign_BP_down_merged <- fgseaRes_sign_BP_down_merged[
  !grepl('GO:0010648|GO:0120031', fgseaRes_sign_BP_down_merged$pathway)] # negative regulation of signaling

fgseaRes_sign_BP_down_merged$pathway <- sub('.*_', '', fgseaRes_sign_BP_down_merged$pathway) # for figures

fgseaRes_sign_BP_down_merged$pathway <- sapply(fgseaRes_sign_BP_down_merged$pathway, 
                                               function(s) paste(strwrap(s, width=45), 
                                                                 collapse='\n'))
fgseaRes_sign_BP_down_merged$pathway <- factor(fgseaRes_sign_BP_down_merged$pathway)
fgseaRes_sign_BP_down_merged$pathway <- reorder(fgseaRes_sign_BP_down_merged$pathway, 
                                                fgseaRes_sign_BP_down_merged$gene_ratio)

ggplot(fgseaRes_sign_BP_down_merged, aes(gene_ratio, pathway)) +
  geom_point(aes(size = size, color = pval)) +
  ylab('') +
  xlab('Gene ratio, %') +
  scale_color_viridis('p-value', option = 'viridis', direction = -1, end = 0.6) + 
  scale_size_continuous('Number of genes', 
                        breaks = c(round(seq(min(fgseaRes_sign_BP_down_merged$size), 
                                             max(fgseaRes_sign_BP_down_merged$size), 
                                             length.out = 4)))) +
  theme_light() +
  theme(axis.text = element_text(size = 12), 
        axis.title.x = element_text(size = 12))

ggsave(filename = file.path(dir_to_write, 'dotplot_Ecy_0.05_BP_down.png'),
       scale = 0.7, width = 11, height = 7)





