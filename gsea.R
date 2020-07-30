library(fgsea)
library(Biobase)
library(ggplot2)
library(ggforce)
library(ggsci)
library(shadowtext)

### Make paths to required files:
dir <- 'labeglo2/Transcriptomics/quantification/AEA_metazoa_Ecy'
#dir <- 'labeglo2/MS_results/fold_change_calculation/second_exp/' # proteomics
#DE_sample <- 'Eve_BKvsPB_1.5C_filt.csv'  # proteomics
#DE_sample <- 'condition_1.5C_vs_12C_0.05_FAv3_1.5vs12_.csv' # the output from differential expression analysis (Note that it is better to take DE transcripts with low padj)
DE_sample <- 'condition_1.5C_vs_12C_0.05_FA.csv'
#DE_sample <- 'condition_Ecy_vs_Gla_0.05_1.5C.csv'
#DE_sample <- 'condition_Ecy_vs_Gla_0.01_lfc2_1.5C.csv'
eggnog_FA <- 'ecy_metazoa_eggnog3.csv' # the output after EggNOG mapper (functional annotation), the file should be without log info from Eggnog (the first 4 rows and last 3 rows)
#eggnog_FA <- 'gla_peptides_eggnog2.csv'
#eggnog_FA <- 'ecy_peptides_eggnog2.csv'
#eggnog_FA <- 'eve_peptides_eggnog2.csv'
eggnog_FA <- 'eve_metazoa_eggnog2.csv'
eggnog_FA <- 'gla_metazoa_eggnog2.csv'
#eggnog_FA <- 'ecy_gla_eve_peptides_eggnog.csv' # proteomics
DE_gene_list <- read.csv(file.path(dir, DE_sample))
#DE_gene_list <- read.csv(file.path(dir, DE_sample), sep = '\t')
go_terms <- read.csv(file.path(dir, eggnog_FA), sep = '\t', header = F)
go_info <- read.csv('labeglo2/Transcriptomics/annotations/functional_annotation/GO_info.csv', sep = '\t', header = T) # tsv file with three columns: GO id, GO description, and root (BC, MF, or CC)

### Files preprocessing:

go_terms <- subset(go_terms, V6 != '') # remove rows which do not contain gene name
GOs <- strsplit(as.character(go_terms$V7), ',') # extract GO terms (this is gene2GO list)
names(GOs) <- go_terms$V6 # assign gene names to it
#GOs <- GOs[unique(names(GOs))] # proteomics
### Take genes and stat from DE analysis. Stat is a z-statistics which is LFC/stderr (Wald test):
DE_gene_stat <- DE_gene_list$stat # transcriptomics
#DE_gene_stat <- DE_gene_list$t.mod # proteomics
names(DE_gene_stat) <- DE_gene_list$X
names(DE_gene_stat) <- rownames(DE_gene_list) # proteomics

DE_gene_lfc <- DE_gene_list$log2FoldChange
names(DE_gene_lfc) <- DE_gene_list$X

### Make GO2gene annotation list and put description of GO terms to GO2gene list:
GOs2gene <- reverseSplit(GOs) # make GO2gene from gene2Go
GOs2gene <- lapply(GOs2gene, unique) # remove duplicates of genes in GO2gene
tmp <- match(names(GOs2gene), go_info$GO_ID) # to select only required GOs from GO_info
names(GOs2gene) <- go_info[tmp,]$Description

### Split data on three categories:
GOs2gene_BP <- GOs2gene[go_info[tmp,]$Root == 'biological_process']
GOs2gene_MF <- GOs2gene[go_info[tmp,]$Root == 'molecular_function']
GOs2gene_CC <- GOs2gene[go_info[tmp,]$Root == 'cellular_component']

#tmp <- sapply(GOs2gene, function(x) length(x) > 3)
#GOs2gene_less <- GOs2gene[tmp]

### Run GSEA
fgseaRes_BP <- fgsea(GOs2gene_BP, DE_gene_stat, nperm = 10000)
fgseaRes_MF <- fgsea(GOs2gene_MF, DE_gene_stat, nperm = 10000)
fgseaRes_CC <- fgsea(GOs2gene_CC, DE_gene_stat, nperm = 10000)
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
fgseaRes_sign_CC <- cbind(subset(fgseaRes_CC, pval < 0.01), 
                          rep('Cellular compartment', 
                          nrow(subset(fgseaRes_CC, pval < 0.01))))
fgseaRes_sign_CC <- subset(fgseaRes_sign_CC, size > 2)

fgseaRes_sign <- rbind(fgseaRes_sign_BP, fgseaRes_sign_MF)
fgseaRes_sign <- rbind(fgseaRes_sign, fgseaRes_sign_CC)
colnames(fgseaRes_sign) <- c(colnames(fgseaRes_sign)[1:8], 'Group')
fgseaRes_sign$Group <- as.factor(fgseaRes_sign$Group)

### Illustration of GSEA

fgseaRes_sign$pathway <- reorder(fgseaRes_sign$pathway, fgseaRes_sign$NES)

ggplot(fgseaRes_sign, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=Group)) +
  coord_flip() +
  scale_fill_d3() +
  labs(x="Pathway", y="Normalized Enrichment Score") + 
#  theme(axis.text.y = element_text(color = fgseaRes_sign$Group))
  theme_light()

ggsave(file.path(dir, 'bars_gsea_1.5vs6C.png'), scale = 1.5)

ggplot(fgseaRes_sign, aes(pathway, NES)) +
  geom_segment(aes(xend = pathway, y = 0, yend = NES, color = NES > 0), size=1) +
  geom_point(aes(color = NES > 0), size = 5, shape=21, fill='white', stroke=1) +
  geom_text(aes(label=size), alpha=0.5, size = 3) +
  geom_text(aes(y=0, label=sprintf(' %s ', pathway)),
            hjust = ifelse(fgseaRes_sign$NES > 0, 1, 0), alpha=0.8, size = 3) +
  coord_flip() +
  scale_x_discrete('', labels = NULL) +
  scale_y_continuous('Normalized Enrichment Score') +
  scale_color_manual('Regulation', values = c('dodgerblue1', 'firebrick1'), labels = c('DOWN','UP')) +
  ggforce::facet_col(fgseaRes_sign$Group, scales = "free_y", space = "free") +
  theme_minimal()

ggsave(file.path(dir, 'lollipop_gsea_facet_0.02.png'), scale = 2.2)

### For big outputs it is better to draw every category in its own plot

# Biological process

fgseaRes_sign_BP$pathway <- reorder(fgseaRes_sign_BP$pathway, fgseaRes_sign_BP$NES)
fgseaRes_sign_BP <- subset(fgseaRes_sign_BP, size >= 2)

ggplot(fgseaRes_sign_BP, aes(pathway, NES)) +
  geom_segment(aes(xend = pathway, y = 0, yend = NES, color = NES > 0), size=1) +
  geom_point(aes(color = NES > 0), size = 5, shape=21, fill='white', stroke=1) +
  geom_text(aes(label=size), alpha=0.5, size = 2) +
  geom_text(aes(y=0, label=sprintf(' %s ', pathway)),
            hjust = ifelse(fgseaRes_sign_BP$NES > 0, 1, 0), alpha=0.8, size = 3) +
  coord_flip() +
  scale_x_discrete('', labels = NULL) +
  scale_y_continuous('Normalized Enrichment Score') +
  scale_color_manual('Regulation', values = c('dodgerblue1', 'firebrick1'), labels = c('DOWN','UP')) +
  theme_minimal()

ggsave(file.path(dir, 'lollipop_gsea_BP_0.002.png'), scale = 3.2)

fgseaRes_sign_BP_up <- subset(fgseaRes_sign_BP, NES > 0)
ggplot(subset(fgseaRes_sign_BP_up, NES > 0), aes(pathway, NES)) +
  geom_segment(aes(xend = pathway, y = 0, yend = NES), size=1, color='firebrick1') +
  geom_point(size = 5, shape=21, fill='firebrick1', stroke=1, color='firebrick1') +
  geom_text(aes(label=size), alpha=1, size = 2) +
  coord_flip() +
  scale_x_discrete('') +
  scale_y_continuous('Normalized Enrichment Score') +
  scale_color_manual('Regulation', values = c('dodgerblue1', 'firebrick1'),
                     labels = c('DOWN', 'UP')) +
  theme_minimal()

ggsave(file.path(dir, 'lollipop_gsea_BP_up_0.002.png'), scale = 1.9)

ggplot(subset(fgseaRes_sign_BP_up, NES > 0), aes(pathway, NES)) +
  geom_bar(stat='identity', fill='#FF0033', alpha=0.7) +
  geom_text(aes(y = -0.2, label=sapply(leadingEdge, length)), alpha=0.7, 
            size = 3.5, hjust=1) +
  annotate('text', x = nrow(fgseaRes_sign_BP_up)+2.2, y = -0.31, 
           label='Leading\nedge', size=3.2, alpha=0.7) +
  coord_flip() +
  scale_x_discrete('', expand=expansion(add=c(0, 3.3))) +
  scale_y_continuous('Normalized Enrichment Score', 
                     limits=c(-0.4, max(fgseaRes_sign_BP_up$NES))) +
  theme_minimal()+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.y = element_text(colour = 'black'))

ggsave(file.path(dir, 'EcyGla0.01lfc2_BP_up_0.001.png'), scale = 1.3, height = 7, width = 7)

# BP with merging

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
#---------------E.cyaneus--BP-merging-UP-----------------------------------------
fgseaRes_sign_BP_up_merged <- fgseaRes_sign_BP_up_merged[!grepl('GO:0051716', fgseaRes_sign_BP_up_merged$pathway)] # response to stimulus
fgseaRes_sign_BP_up_merged <- fgseaRes_sign_BP_up_merged[!grepl('GO:0051603|GO:0044257', fgseaRes_sign_BP_up_merged$pathway)] # protein catabolic process
fgseaRes_sign_BP_up_merged <- fgseaRes_sign_BP_up_merged[!grepl('GO:0019941|GO:0006511|GO:0010498', fgseaRes_sign_BP_up_merged$pathway)] # modification-dependent macromolecule catabolic process
fgseaRes_sign_BP_up_merged <- fgseaRes_sign_BP_up_merged[!grepl('GO:0048878', fgseaRes_sign_BP_up_merged$pathway)] # homeostasis process
fgseaRes_sign_BP_up_merged <- fgseaRes_sign_BP_up_merged[!grepl('GO:0000377|GO:0000375|GO:0000398|GO:0006397', fgseaRes_sign_BP_up_merged$pathway)] # RNA splicing
fgseaRes_sign_BP_up_merged <- fgseaRes_sign_BP_up_merged[!grepl('GO:0003002', fgseaRes_sign_BP_up_merged$pathway)] # pattern specification process
fgseaRes_sign_BP_up_merged <- fgseaRes_sign_BP_up_merged[!grepl('GO:0044093|GO:0051347|GO:0033674', fgseaRes_sign_BP_up_merged$pathway)] # positive regulation of catalytic activity
#---------------E.verrucosus--BP-merging-UP--------------------------------------------
fgseaRes_sign_BP_up_merged <- fgseaRes_sign_BP_up_merged[!grepl('GO:0007338', fgseaRes_sign_BP_up_merged$pathway)] # fertilization
fgseaRes_sign_BP_up_merged <- fgseaRes_sign_BP_up_merged[!grepl('GO:0050000', fgseaRes_sign_BP_up_merged$pathway)] # establishment of chromosome localization
fgseaRes_sign_BP_up_merged <- fgseaRes_sign_BP_up_merged[!grepl('GO:0071826|GO:0022618', fgseaRes_sign_BP_up_merged$pathway)] # ribonucleoprotein complex biogenesis
fgseaRes_sign_BP_up_merged <- fgseaRes_sign_BP_up_merged[!grepl('GO:0090304|GO:0051259|GO:0034622|GO:0061077', fgseaRes_sign_BP_up_merged$pathway)]
fgseaRes_sign_BP_up_merged <- fgseaRes_sign_BP_up_merged[!grepl('GO:0009792', fgseaRes_sign_BP_up_merged$pathway)] # embryo development
fgseaRes_sign_BP_up_merged <- fgseaRes_sign_BP_up_merged[!grepl('GO:0065003', fgseaRes_sign_BP_up_merged$pathway)] # protein-containing complex subunit organization
#---------------Gla-BP-merging-1.5vs6C---UP
fgseaRes_sign_BP_up_merged <- fgseaRes_sign_BP_up_merged[
  !grepl('GO:0000377|GO:0000375', fgseaRes_sign_BP_up_merged$pathway)] # RNA splicing
###
fgseaRes_sign_BP_up_merged$pathway <- sub('.*_', '', 
                                          fgseaRes_sign_BP_up_merged$pathway) # for figures
fgseaRes_sign_BP_up_merged$pathway <- reorder(fgseaRes_sign_BP_up_merged$pathway, 
                                              fgseaRes_sign_BP_up_merged$NES)

fgseaRes_sign_BP_up_merged$pathway <- sapply(fgseaRes_sign_BP_up_merged$pathway, 
                                             function(s) paste(strwrap(s, width=45), 
                                                               collapse='\n'))
fgseaRes_sign_BP_up_merged$pathway <- factor(fgseaRes_sign_BP_up_merged$pathway)
fgseaRes_sign_BP_up_merged$pathway <- reorder(fgseaRes_sign_BP_up_merged$pathway, 
                                              fgseaRes_sign_BP_up_merged$NES)

ggplot(subset(fgseaRes_sign_BP_up_merged, NES > 0), aes(pathway, NES)) +
  geom_bar(stat='identity', fill='#FF0033', alpha=0.7) +
  geom_text(aes(y = -0.2, label=sapply(leadingEdge, length)), 
            alpha=0.7, size = 3, hjust=1) +
  annotate('text', x = nrow(fgseaRes_sign_BP_up_merged)+1.4, y = -0.31, 
           label='Leading\nedge', size=3.2, alpha=0.7) +
  coord_flip() +
  scale_x_discrete('', expand=expansion(add=c(0, 3.3))) +
  scale_y_continuous('NES', 
                     limits=c(-0.5, max(fgseaRes_sign_BP_up_merged$NES))) +
  theme_minimal()+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.y = element_text(colour = 'black'))

ggsave(file.path(dir, 'bars_gsea_1.5vs12C_BP_up_0.05_merged_readable.png'), 
       scale = 0.75, width = 6.5, height = 6)

chosen_contigs_to_heatmap <- DE_gene_lfc[
  unlist(
    fgseaRes_sign_BP_up_merged[grepl('homeostatic process', 
                                     fgseaRes_sign_BP_up_merged$pathway)]$leadingEdge
    )
  ]
chosen_contigs_to_heatmap_stat <- DE_gene_stat[
  unlist(
    fgseaRes_sign_BP_up_merged[grepl('^(RNA processing|mRNA metabolic process)$', 
                                     fgseaRes_sign_BP_up_merged$pathway)]$leadingEdge
  )
  ]

write.csv(names(chosen_contigs_to_heatmap_stat), 'GO:0010467_gene_list.txt', 
          quote = F, row.names = F)
## DOWN 
fgseaRes_sign_BP_down <- subset(fgseaRes_sign_BP, NES < 0)
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
#-------------------Ecy--BP-merging--DOWN-----------------------------------------------------------
fgseaRes_sign_BP_down_merged <- fgseaRes_sign_BP_down_merged[
  !grepl('GO:0032269|GO:0001933|GO:0042326|GO:0048523|GO:0010605|GO:0009890|GO:0045936|GO:0031324|GO:0031327|GO:0051248|GO:0031400|GO:0010558|GO:0010563|GO:0051172', 
         fgseaRes_sign_BP_down_merged$pathway)] # negative regulation of metabolic process
fgseaRes_sign_BP_down_merged <- fgseaRes_sign_BP_down_merged[
  !grepl('GO:0007268|GO:0098916|GO:0099537|GO:0099565|GO:0007267', fgseaRes_sign_BP_down_merged$pathway)] # synaptic signaling
fgseaRes_sign_BP_down_merged <- fgseaRes_sign_BP_down_merged[
  !grepl('GO:0010648|GO:0120031', fgseaRes_sign_BP_down_merged$pathway)] # negative regulation of signaling
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
#--------------------
#--------G.la--BP-merging-15vs6C---
fgseaRes_sign_BP_down_merged <- fgseaRes_sign_BP_down_merged[
  !grepl('GO:0000956|GO:0006402|GO:0000184|GO:0034655', 
         fgseaRes_sign_BP_down_merged$pathway)] # RNA catabolic process
fgseaRes_sign_BP_down_merged <- fgseaRes_sign_BP_down_merged[
  !grepl('GO:0019439|GO:0044270|GO:0046700', 
         fgseaRes_sign_BP_down_merged$pathway)] # organic cyclic compound catabolic process
fgseaRes_sign_BP_down_merged <- fgseaRes_sign_BP_down_merged[
  !grepl('GO:1901293|GO:0009260|GO:0006164|GO:0046390|GO:0009152|GO:0072522', 
         fgseaRes_sign_BP_down_merged$pathway)] # nucleotide biosynthetic process
fgseaRes_sign_BP_down_merged <- fgseaRes_sign_BP_down_merged[
  !grepl('GO:0006518|GO:0043604|GO:0043043', 
         fgseaRes_sign_BP_down_merged$pathway)] # translation
fgseaRes_sign_BP_down_merged <- fgseaRes_sign_BP_down_merged[
  !grepl('GO:0006612|GO:0006613|GO:0006614|GO:0045047', 
         fgseaRes_sign_BP_down_merged$pathway)] # protein targeting
fgseaRes_sign_BP_down_merged <- fgseaRes_sign_BP_down_merged[
  !grepl('GO:0072599', 
         fgseaRes_sign_BP_down_merged$pathway)] # localization to EPR
fgseaRes_sign_BP_down_merged <- fgseaRes_sign_BP_down_merged[
  !grepl('GO:0098916|GO:0099536|GO:0007268', 
         fgseaRes_sign_BP_down_merged$pathway)] # synaptic sygnaling
fgseaRes_sign_BP_down_merged <- fgseaRes_sign_BP_down_merged[
  !grepl('GO:0070588|GO:0098655|GO:0098660|GO:0015672|GO:1902600', 
         fgseaRes_sign_BP_down_merged$pathway)] # inorganic cation transmembrane transport
fgseaRes_sign_BP_down_merged <- fgseaRes_sign_BP_down_merged[
  !grepl('GO:0090150|GO:0002181|GO:0000028|GO:0060271|GO:0030031', fgseaRes_sign_BP_down_merged$pathway)] 
#----------------------------------------------------
fgseaRes_sign_BP_down_merged$pathway <- sub('.*_', '', fgseaRes_sign_BP_down_merged$pathway) # for figures
fgseaRes_sign_BP_down_merged$pathway <- reorder(fgseaRes_sign_BP_down_merged$pathway, fgseaRes_sign_BP_down_merged$NES)
fgseaRes_sign_BP_down_merged$pathway <- sapply(fgseaRes_sign_BP_down_merged$pathway, 
                                             function(s) paste(strwrap(s, width=45), 
                                                               collapse='\n'))
fgseaRes_sign_BP_down_merged$pathway <- factor(fgseaRes_sign_BP_down_merged$pathway)
fgseaRes_sign_BP_down_merged$pathway <- reorder(fgseaRes_sign_BP_down_merged$pathway, 
                                              fgseaRes_sign_BP_down_merged$NES)
ggplot(fgseaRes_sign_BP_down_merged, aes(pathway, NES)) +
  geom_bar(stat='identity', fill='blue2', alpha=0.7) +
  geom_text(aes(y = 0.4, label=sapply(leadingEdge, length)), 
            alpha=0.7, size = 3.5, hjust=1) +
  annotate('text', x = nrow(fgseaRes_sign_BP_down_merged)+1.5, y = 0.25, 
           label='Leading\nedge', size=3.8, alpha=0.7) +
  coord_flip() +
  scale_x_discrete('', expand=expansion(add=c(0, 3.3))) +
  scale_y_continuous('NES', 
                     limits=c(min(fgseaRes_sign_BP_down_merged$NES), 0.5)) +
  theme_minimal()+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.y = element_text(colour = 'black', size = 10))

ggsave(file.path(dir, 'bars_gsea_1.5vs12C_BP_down_0.05_merged_readable3.png'), 
       scale = 0.7, width = 8, height = 13)

chosen_contigs_to_heatmap_stat <- DE_gene_stat[
  unlist(
    fgseaRes_sign_BP_down_merged[grepl('^(sperm motility)$', 
                                     fgseaRes_sign_BP_down_merged$pathway)]$leadingEdge
  )
  ]
# Molecular function  

fgseaRes_sign_MF$pathway <- reorder(fgseaRes_sign_MF$pathway, fgseaRes_sign_MF$NES)

ggplot(fgseaRes_sign_MF, aes(pathway, NES)) +
  geom_segment(aes(xend = pathway, y = 0, yend = NES, color = NES > 0), size=1) +
  geom_point(aes(color = NES > 0), size = 5, shape=21, fill='white', stroke=1) +
  geom_text(aes(label=size), alpha=0.5, size = 2.5) +
  geom_text(aes(y=0, label=sprintf(' %s ', pathway)),
            hjust = ifelse(fgseaRes_sign_MF$NES > 0, 1, 0), alpha=0.8, size = 2.5) +
  coord_flip() +
  scale_x_discrete('', labels = NULL) +
  scale_y_continuous('Normalized Enrichment Score') +
  scale_color_manual('Regulation', values = c('dodgerblue1', 'firebrick1'), 
                     labels = c('DOWN','UP')) +
  theme_minimal()

ggsave(file.path(dir, 'lollipop_gsea_MF.png'), scale = 1.5)

fgseaRes_sign_MF <- fgseaRes_sign_MF[
  !grepl('GO:0008186|GO:0016411|GO:0005528|GO:0004553|GO:0051959', fgseaRes_sign_MF$pathway)] 

fgseaRes_sign_MF$pathway <- sapply(fgseaRes_sign_MF$pathway, function(s)
  paste(strwrap(s, width=75), collapse='\n'))

fgseaRes_sign_MF$pathway <- factor(fgseaRes_sign_MF$pathway)
fgseaRes_sign_MF$pathway <- reorder(fgseaRes_sign_MF$pathway, fgseaRes_sign_MF$NES)

ggplot(fgseaRes_sign_MF, aes(pathway, NES)) +
  geom_bar(aes(fill=NES>0), stat='identity') +
  geom_text(aes(y = ifelse(NES>0, -0.3, 0.3), label=sapply(leadingEdge, length)), 
            alpha=0.7, size = 3.5, hjust=0.5) +
  annotate('text', x = nrow(fgseaRes_sign_MF)+1.5, y = 0, 
           label='Leading\nedge', size=3.2, alpha=0.7) +
  coord_flip() +
  scale_x_discrete('', expand=expansion(add=c(0, 3.3))) +
  scale_y_continuous('Normalized Enrichment Score', 
                     limits=c(min(fgseaRes_sign_MF$NES)-0.3, max(fgseaRes_sign_MF$NES))) +
  scale_fill_manual('Regulation', values = alpha(c('blue2', '#FF0033'), 0.5), 
                     labels = c('DOWN','UP')) +
  theme_minimal()+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.y = element_text(colour = 'black'))

ggsave(file.path(dir, 'EcyGlaDE0.01lfc2_MF_0.05.png'), scale = 1.35, height = 6, width = 7)

chosen_contigs_to_heatmap <- DE_gene_stat[
  unlist(
    fgseaRes_sign_MF[grepl('calmodulin binding', 
                           fgseaRes_sign_MF$pathway)]$leadingEdge
  )
  ]
#----------Gla-MF-Merging-----------------------------------------------------------------
fgseaRes_sign_MF <- fgseaRes_sign_MF[
  !grepl('GO:0008569|GO:0003777|GO:1990939', fgseaRes_sign_MF$pathway)] #motor activity
fgseaRes_sign_MF <- fgseaRes_sign_MF[
  !grepl('GO:0005254|GO:1901363|GO:0097159', fgseaRes_sign_MF$pathway)]
#---------------------------------------------------------------------------------
# Cellular component

fgseaRes_sign_CC$pathway <- reorder(fgseaRes_sign_CC$pathway, fgseaRes_sign_CC$NES)

ggplot(fgseaRes_sign_CC, aes(pathway, NES)) +
  geom_segment(aes(xend = pathway, y = 0, yend = NES, color = NES > 0), size=1) +
  geom_point(aes(color = NES > 0), size = 5, shape=21, fill='white', stroke=1) +
  geom_text(aes(label=size), alpha=0.5, size = 2.5) +
  geom_text(aes(y=0, label=sprintf(' %s ', pathway)),
            hjust = ifelse(fgseaRes_sign_CC$NES > 0, 1, 0), alpha=0.8, size = 2.5) +
  coord_flip() +
  scale_x_discrete('', labels = NULL) +
  scale_y_continuous('Normalized Enrichment Score') +
  scale_color_manual('Regulation', values = c('dodgerblue1', 'firebrick1'), labels = c('DOWN','UP')) +
  theme_minimal()

ggsave(file.path(dir, 'lollipop_gsea_CC.png'), scale = 2)

## Prepare table for Revigo

test <- fgseaRes_sign[,1:2]
test$pathway <- sub('_.*', '', test$pathway)
write.table(test, 
          'labeglo2/Transcriptomics/quantification/AEA_metazoa_Ecy/Ecy_for_Revigo.csv',
          sep = '\t', quote = F, row.names = F)

test <- fgseaRes_sign_BP_down_merged[,1:2]
test$pathway <- sub('_.*', '', test$pathway)
write.table(test, file.path(dir, 'Gla_BPdown_1.5vs6C_for_Revigo.csv'),
            sep = '\t', quote = F, row.names = F)


