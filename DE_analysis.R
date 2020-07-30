library(tximport)
library(DESeq2)
library(apeglm)
library(RColorBrewer)
library(pheatmap)
library(viridis)
library(ggplot2)
library(ggrepel)
library(gplots)

### Make paths to the samples

path_to_files <- function(dir_to_samples, sample_dir_names){
  f <- file.path(dir_to_samples, sample_dir_names, 'quant.sf')
  names(f) <- paste0('sample', 1:length(sample_dir_names))
  return(f)
}

dir <- 'labeglo2/Transcriptomics/quantification/AEA_metazoa_Gla' # specify path to your samples
dir <- 'labeglo2/Transcriptomics/quantification/AEA_GlaReads_EcyAssembly/GlaReads_EcyAssembly'
dir <- 'labeglo2/Transcriptomics/quantification/GlaReads_EcyAssembly_all_samples/'

#files_names <- c('350', '351',  '352', '354', '632', '633',  '634',  '636',
#                 'Eve_6C_rep1', 'Eve_6C_rep3', 'Eve_6C_rep4',
#                 'Eve_24C_rep1','Eve_24C_rep2','Eve_24C_rep3', 'Eve_24C_rep4') # specify the names of folders with quant.sf data 

files_names <- c('311', '312', '313', '315', '594', '595', '596', '598') # specify the names of folders with quant.sf data, Gla
files_names <- c('129', '130', '131', '133', '560', '561', '562', '563') # Ecy
files_names <- c('350', '351',  '352', '354', '632', '633',  '634',  '636') # Eve
files_names <- c('594_glaA', '595_glaA', '596_glaA', '598_glaA', 
                 '560_glaA', '561_glaA', '562_glaA', '563_glaA') # Gla 1.5 and Ecy 1.5
files_names <- c('594_ecyA', '595_ecyA', '596_ecyA', '598_ecyA', 
                 '560_ecyA', '561_ecyA', '562_ecyA', '563_ecyA')
files <- path_to_files(dir, files_names)
all(file.exists(files))

annot_file <- file.path(dir, 'contigs_functional_annot.txt')
annot_file <- file.path(dir, 'contigs_DAnnotation.txt')

### Read table with contigs ids and genes they assigned 

txi_creation <- function(dir_to_samples, 
                         file_with_annotation, 
                         path_to_files_, 
                         type_of_quant="salmon"){
  c2g <- read.csv(file.path(dir_to_samples, file_with_annotation),
                  sep = '\t', header = F, col.names = c('ID', 'Gene', 'Description'))
  txi_ <- tximport(path_to_files_, type = type_of_quant, tx2gene = c2g) # tximport sums up TPMs of the same genes if it meets the same genes name; txi$abundunce - TPMs, txi$counts - number of reads, txi$length - effective length
  return(txi_)
}

txi <- txi_creation(dir, 'contigs_functional_annot.txt', files)
txi <- txi_creation(dir, 'contigs_DAnnotation.txt', files)
txi <- txi_creation(dir, 'Gla_kegg_pred_cds_annotation_wo_p_and_dupl.txt', files)
txi <- txi_creation(dir, 'kegg_annot_Gla_assembly.txt', files)
txi <- txi_creation(dir, 'pred_peptides_functional_annot_wo_p_and_Dupl.txt', files)
txi$counts <- txi$counts[-1,]

#write.table(txi$counts, 'labeglo2/Transcriptomics/quantification/AEA_metazoa_Ecy/txi_count_1.5_12.tsv', 
#            sep='\t', quote = F)
### Prepare table with information about samples (it needs for DESeq running)

table_with_sample_descr <- function(samples_number, conditions){ # conditions is a vector of conditions in the order: untreated, treated
  d <- data.frame(row.names = paste0('sample', 1:samples_number), 
    condition = rep(conditions, each = samples_number/length(conditions)),
    replicate = rep(paste0('replicate_', 1:(samples_number/length(conditions))), 
                    length(conditions)))
  d$condition <- factor(d$condition, levels = conditions)
  return(d)
}

sample_description <- table_with_sample_descr(8, c('12C', '1.5C'))
sample_description <- table_with_sample_descr(8, c('Gla', 'Ecy'))
sample_description$condition <- factor(sample_description$condition, 
                                       levels = c('Ecy', 'Gla'))
sample_description$condition <- factor(sample_description$condition, 
                                       levels = c('1.5C', '12C', '6C', '24C'))
sample_description$condition <- ifelse(sample_description$condition == '1.5C', 
                                       '1.5C', '6or12C')
sample_description$replicate <- c(paste0('12C_', 1:4), paste0('1.5C_', 1:4), 
                                  paste0('6C_', 1:3), paste0('24C_', 1:4))
sample_description$condition <- factor(sample_description$condition, 
                                       levels = c('6or12C','1.5C'))
sample_description <- sample_description[!rownames(sample_description) == 'sample10',]
row.names(sample_description) <- paste0('sample', 1:7)

### Run differential expression analysis based on the negative binominal distribution

deseq <- DESeqDataSetFromTximport(txi, colData = sample_description, 
                                      design = ~ condition)
deseq$condition

deseq <- DESeq(deseq)
res <- results(deseq, name = resultsNames(deseq)[2]) # results() extracts a results table with log2 fold changes, p values and adjusted p values  
#res[grepl('XP_015433778.1|RNA01364.1|RNA33404.1', row.names(res)),] # select set of genes
res[grepl('OXCT1|OXCT|SCOT|CG1140|catI', row.names(res), ignore.case = T),]
res[grepl('CYB5R3|DIA1|Cyb5r3|MCR1|YKL150W|YKL605|CBR1|CBR|CYB5R', row.names(res), ignore.case = T),] # NADH-cytochrome b5 reductase -> one of the component of aerobic desaturation (PUFA synthesis)
res[grepl('CYB5A|CYB5B|CYB5m|Omb5|CB5-A|N1949|YNL111C|Cyt-b5|CG2140', row.names(res), ignore.case = T),] # cytochrome b5 -> one of the component of aerobic desaturation (PUFA synthesis)
res[grepl('FADS|SCD|des6|Delta6FAD|FAD2|desat|OLE1|fat-6', row.names(res), ignore.case = T),] # acyl-coa desaturases -> one of the component of aerobic desaturation (PUFA synthesis)
res[grepl('GPT2|AGXT2|GFPT|OAT', row.names(res), ignore.case = T),]
res[grepl('Vg', row.names(res), ignore.case = F),]
res[grepl('GLUD|Gdh', row.names(res), ignore.case = T),]

### Pathways

kegg_names <- read.csv('labeglo2/Transcriptomics/annotations/pathway_analisys/Gla/endocytosis_kegg_names.txt', header = F) # file with list of kegg names from the pathway of interest
# Endocytosis
pathway_res <- res[grepl('AP2A2|arf-1.2|ARF4|ARFGEF1|ARPC2|ARPC4|ARPC5|CAPZB|CCDC53|CHMP2A|CHMP2B|CHMP6|CHMP7|CLINT1|DAB2|DNM2|dyn-1|efa-6|EHD4|EPN2|FAM21A|fgr|GBF1|HGS|Hsp68|HSPA1A|IST1|KIAA0196|PARD3|PDCD6IP|PIP5K1B|RAB10|RAB11A|RAB11FIP2|Rab3|Rab7|rab8a|SH3KBP1|SMAD2|SNX4|Sop2|TSG101|WASHC1|ZFYVE20', row.names(res)),]
pathway_res <- res[grepl('ACD88959.1|ACI90359.1|ACO11926.1|ELU01524.1|ELU18512.1|GBL84746.1|GBM03826.1|GBM49310.1|KAA0187646.1|KAE9436284.1|KAE9529754.1|KAF0313898.1|KOF99293.1|KPM10047.1|KPP61432.1|KRZ20197.1|NP_001079320.1|NP_001171767.1|ODN03957.1|OQR75673.1|P21613.1|PRD35706.1|PVD38106.1|QGC85414.1|RMZ98648.1|RNA06970.1|RNA10777.1|RNA12143.1|RNA23416.1|RNA30539.1|RNA36530.1|RNA40271.1|RNA40936.1|RNA42646.1|RNA43200.1|RNA45177.1|RUS82360.1|RWS17700.1|TDH14235.1|XP_001634219.1|XP_001947218.2|XP_002108937.1|XP_002608872.1|XP_002734355.1|XP_003387703.1|XP_005092492.1|XP_005092554.1|XP_007260027.2|XP_007552695.1|XP_008150111.2|XP_009016976.1|XP_009051123.1|XP_009060803.1|XP_009063300.1|XP_011423452.1|XP_011435429.1|XP_011436985.1|XP_011448123.1|XP_011448246.1|XP_011682840.2|XP_013089460.1|XP_013091222.1|XP_013145319.1|XP_013387153.1|XP_013403219.1|XP_013410030.1|XP_013413854.1|XP_013414079.1|XP_014788045.1|XP_015723034.1|XP_015753214.1|XP_015784105.1|XP_015833074.1|XP_015911282.1|XP_015913996.1|XP_015923451.1|XP_015929739.1|XP_015930111.1|XP_017280720.1|XP_018011276.1|XP_018012068.1|XP_018015153.1|XP_018015160.1|XP_018015331.1|XP_018020317.1|XP_018020389.1|XP_018022226.1|XP_018024609.1|XP_018025977.1|XP_018323433.1|XP_018595174.1|XP_018902622.1|XP_019622749.1|XP_019635230.1|XP_019926929.1|XP_020389187.1|XP_020470088.1|XP_020615495.1|XP_020775282.1|XP_020897630.1|XP_021350910.1|XP_021371058.1|XP_021376808.1|XP_021938486.1|XP_022109833.1|XP_022120296.1|XP_022293720.1|XP_022299881.1|XP_022314593.1|XP_022344077.1|XP_022666255.1|XP_022700037.1|XP_023165599.2|XP_023220284.1|XP_023242330.1|XP_023706140.1|XP_023933170.1|XP_024920344.1|XP_025076692.1|XP_025081538.1|XP_025089563.1|XP_026470738.1|XP_026544589.1|XP_026809367.1|XP_027002072.1|XP_027049103.1|XP_027217697.1|XP_027232977.1|XP_028916161.1|XP_029467529.1|XP_029646664.1|XP_029939337.1|XP_029981088.1|XP_030847587.1|XP_030848930.1', row.names(res)),] 
#write.table(pathway_res, file = file.path(dir, 'endocytosis_DE_genes_1.5vs12.csv'), quote = F, sep = '\t')
pathway_res <- pathway_res[row.names(pathway_res) != 'OQR75673.1',]
colnames(pathway_res)
pathway_res_with_kegg <- read.csv('labeglo2/Transcriptomics/annotations/pathway_analisys/Gla/endocytosis_DE_genes_1.5vs12_KEGG_NAMES.csv', 
                                  sep = '\t')
# mTOR
pathway_res <- res[grepl('AAB51350.1|ACD37540.1|ACD37572.1|ARU12812.1|CRL08365.1|NP_001164247.1|NP_001191427.1|OPJ87379.1|OQV24379.1|PSN51071.1|QBB01716.1|RNA11603.1|RNA22546.1|RNA33212.1|RUS76505.1|RWS16001.1|RWS29667.1|RXG67127.1|VVC28224.1|XP_005996387.1|XP_006132953.1|XP_006788409.1|XP_009067006.1|XP_012522263.1|XP_013387556.1|XP_013396754.1|XP_013407885.1|XP_013416422.1|XP_013416707.1|XP_013774923.1|XP_014673959.1|XP_014779040.1|XP_014783093.1|XP_015919227.1|XP_018019914.1|XP_018022048.1|XP_018022624.1|XP_018027170.1|XP_018319431.1|XP_018565875.1|XP_019644561.1|XP_019875529.1|XP_020942410.1|XP_021185112.1|XP_021343299.1|XP_021378961.1|XP_022249574.1|XP_022338000.1|XP_022341054.1|XP_022912893.1|XP_023705414.1|XP_023933170.1|XP_024892256.1|XP_025076692.1|XP_025078135.1|XP_025081196.1|XP_025098693.1|XP_025114172.1|XP_025711122.1|XP_026293028.1|XP_027018255.1|XP_027030209.1|XP_028147582.1|XP_028319273.1|XP_028404561.1|XP_029022103.1|XP_029634367.1|XP_029635867.1|XP_029838482.1|XP_030068144.1|XP_031350878.1', row.names(res)),]
write.table(pathway_res, file = file.path(dir, 'mTOR_DE_genes_1.5vs12.csv'), quote = F, sep = '\t')

# Ras
pathway_res <- res[grepl('ACI90359.1|ACT88125.1|AEN94436.1|EFX88843.1|ELT88415.1|ELT92419.1|GBL91632.1|KHJ40867.1|OQR78650.1|QBB01734.1|RNA01031.1|RNA11603.1|RWS16372.1|TMW48634.1|XP_002425376.1|XP_002600468.1|XP_003387703.1|XP_006499573.1|XP_013089460.1|XP_013404807.1|XP_013407885.1|XP_013409440.1|XP_013410607.1|XP_014662426.1|XP_015247035.1|XP_015929739.1|XP_017893513.1|XP_018011064.1|XP_018015054.1|XP_018015883.1|XP_019622515.1|XP_019875529.1|XP_021341911.1|XP_023320056.1|XP_023706140.1|XP_023933170.1|XP_025076692.1|XP_025078135.1|XP_025112658.1|XP_025114649.1|XP_026084923.1|XP_026101014.1|XP_027018255.1|XP_027048263.1|XP_027210598.1|XP_028147582.1|XP_029634367.1|XP_029791923.1|XP_031350878.1', row.names(res)),]
write.table(pathway_res, file = file.path(dir, 'ras_DE_genes_1.5vs12.csv'), quote = F, sep = '\t')

# mRNA surveillance pathway
pathway_res <- res[grepl('AEC22818.1|AXF35703.1|AXF35717.1|CAX73342.1|ELU14652.1|ETN63045.1|GBM57768.1|KAB7506016.1|KAE9424857.1|KAF0287225.1|KFM57642.1|KFM79198.1|KFP38056.1|RMZ95943.1|RNA03132.1|RNA16196.1|RNA27667.1|RNA37021.1|RNA37638.1|RUS80738.1|SVE91448.1|VDO99168.1|XP_002414033.1|XP_003402681.1|XP_003995043.1|XP_004549214.1|XP_006824033.1|XP_008944423.1|XP_009062899.1|XP_011415665.1|XP_013380004.1|XP_013384817.1|XP_013389129.1|XP_013410745.1|XP_013413581.1|XP_015930587.1|XP_016394645.1|XP_018006606.1|XP_018006695.1|XP_018010441.1|XP_018011000.1|XP_018013254.1|XP_018015639.1|XP_018016380.1|XP_018019677.1|XP_018021902.1|XP_018023050.1|XP_018024083.1|XP_018027579.1|XP_019629205.1|XP_019713648.1|XP_020455096.1|XP_021197009.1|XP_022256823.1|XP_022785339.1|XP_022799141.1|XP_023227579.1|XP_023231017.1|XP_023931847.1|XP_026287481.1|XP_026495179.1|XP_027036616.1|XP_027217163.1|XP_028809182.1|XP_029201949.1|XP_029351959.1|XP_029658526.1|XP_029843320.1|XP_029948262.1|XP_030285256.1', row.names(res)),]
write.table(pathway_res, file = file.path(dir, 'RNAsurv_DE_genes_1.5vs12.csv'), quote = F, sep = '\t')

### Paint Pathway

res_sign <- subset(res, padj < 0.05)
res_ordered <- res_sign[order(res_sign$log2FoldChange),]
res_sign_lfc <- subset(res_sign, log2FoldChange > 2)
rownames(res_sign_lfc)
#res_sign_lfc_up <- subset(res_sign_lfc, log2FoldChange > 0)
#res_sign_lfc_down <- subset(res_sign_lfc, log2FoldChange < 0)
#write.table(res_sign_lfc_up, file = file.path(dir, 'FA_UP_0.001_lfc2_DE_genes_1.5vs12.csv'), quote = F, sep = '\t')
#write.table(res_sign_lfc_down, file = file.path(dir, 'FA_DOWN_0.001_lfc2_DE_genes_1.5vs12.csv'), quote = F, sep = '\t')

### Calculate sample-to-sample distances and draw a heatmap of these distances

tr_counts <- vst(deseq) # to perform variance-stabilizing transformation (and get homoskedastic values, i.e. constant variance along the range of mean values)  

ss_dist_heatmap <- function(vst_output, dir_to_samples, name=NULL){
  vst_m <- assay(vst_output) # to get matrix with transformed counts from vst() object
  sample_dist <- dist(t(vst_m), method = "euclidean") # t - transform data to samples-genes matrix (from gene-samples)
  sampleDistMatrix <- as.matrix(sample_dist)
  rownames(sampleDistMatrix) <- paste(vst_output$condition, 
                                      vst_output$replicate, sep = '_')
  colnames(sampleDistMatrix) <- NULL
  condition_name <- paste(levels(vst_output$condition)[2], 
                          levels(vst_output$condition)[1], sep = '_vs_')
  rownames(sampleDistMatrix) <- sub('C.*', '°C', rownames(sampleDistMatrix))
  
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sample_dist,
           clustering_distance_cols = sample_dist,
           color = plasma(10), 
           fontsize_row = 16,
           filename = file.path(dir_to_samples, paste0(name, '_', condition_name, '_heatmap_SampleDIST.png')),
           width = 7,
           height = 5.5)
  if (dev.cur() != 1) {dev.off()}
}
ss_dist_heatmap(tr_counts, dir, name = 'EcyGla2_1.5') # !!! put specific name to not overwrite previous plots

### Draw PCA of samples

PCA_samples <- function(vst_output, dir_to_samples, name=NULL){
  pcaData <- plotPCA(vst_output, intgroup = 'condition', returnData = TRUE)
  pcaData$truename <- c(paste(vst_output$condition, vst_output$replicate, sep = '_'))
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  condition_name <- paste(levels(vst_output$condition)[2], 
                          levels(vst_output$condition)[1], sep = '_vs_')
  g <- ggplot(pcaData, aes(PC1, PC2, color = condition)) +
    geom_point(size = 3) +
    geom_text_repel(aes(y = PC2, label = truename), point.padding = 0.1,
                    segment.alpha = 0) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed() +
    scale_x_continuous(limits = c(min(pcaData$PC1)*1.2, max(pcaData$PC1)*1.2)) +
    scale_color_manual(values = c('#3A9944', '#7C5CCC')) +
    theme_light()
  ggsave(file.path(dir_to_samples, paste0(name, '_', condition_name, 
                                          '_PCA_samples.png')), g, 
         width = 7, height = 6, scale = 1)
}
PCA_samples(tr_counts, dir, name = 'EcyGla_1.5')

### Check for Outliers (visualisation of Cook's distances)

cooks_dist <- function(dds, vst_output, dir_to_samples, name=NULL){
  cooks_log <- as.data.frame(log10(assays(dds)[["cooks"]]))
  cooks_log$gene <- row.names(cooks_log)
  cl <- reshape(cooks_log,
                direction = "long",
                varying = list(names(cooks_log)[1:(length(names(cooks_log))-1)]),
                v.names = 'Cook_values',
                timevar = 'Sample',
                times = paste(vst_output$condition, vst_output$replicate, sep = '_'),
                idvar = 'gene')
  condition_name <- paste(levels(vst_output$condition)[2], 
                          levels(vst_output$condition)[1], sep = '_vs_')
  g <- ggplot(cl) +
    geom_boxplot(aes(x = Sample, y = Cook_values)) +
    theme_light() +
    ylab("Cook's distance, log10")
  ggsave(file.path(dir_to_samples, paste0(name, '_', condition_name,
                                          '_CookDist_boxplot.png')), g, 
         width = 8,
         height = 5,
         scale = 1) 
}

cooks_dist(deseq, tr_counts, dir, name = 'EcyGla_1.5')

### Get some summary from the results:

get_summary_DE <- function(DE_results){
  s <- data.frame()
  for(pv in c(0.1, 0.05, 0.01, 0.005, 0.001)) {
    temp_res <- subset(DE_results, padj < pv)
    Up <- sum(temp_res$log2FoldChange > 0)
    Down <- sum(temp_res$log2FoldChange < 0)
    Sum <- sum(DE_results$padj < pv, na.rm = TRUE)
    s <- rbind(s, c(paste0('< ', pv), 'UP', Up, round(Up/Sum * 100, 2)), stringsAsFactors=F) 
    s <- rbind(s, c(paste0('< ', pv), 'DOWN', Down, round(Down/Sum * 100, 2)), stringsAsFactors=F)
  }
  colnames(s) <- c('Threshold', 'Regulation', 'Transcripts_number', 'Perc')
  s$Transcripts_number <- as.numeric(s$Transcripts_number)
  s$Perc <- as.numeric(s$Perc)
  return(s)
} # it gives number of Up- and Down- regulated genes for padj < 0.1, 0.05, 0.01, 0.005, 0.001
res_sum <- get_summary_DE(res)

summary_visual <- function(output_get_summary_DE, dir_to_samples, name=NULL){
  colors <- c('#a05195','#ffa600')
  g <- ggplot(output_get_summary_DE, aes(Threshold, Transcripts_number)) +
    geom_col(aes(fill = Regulation), position = position_stack(reverse = TRUE),
             color='grey10') +
    coord_flip() +
    scale_fill_manual(values = colors) +
    xlab("Adjusted P value") +
    ylab("Number of transcripts") +
    theme_light()
  ggsave(file.path(dir_to_samples, paste0(name, '_results_summary.png')), g, 
         width = 9, height = 4, scale = 0.8)
  }

summary_visual(res_sum, dir, name = 'EcyGla_1.5C')

### log2 fold change visualisation (MA-plot)

get_ma_plot <- function(dds, dir_to_samples, name=NULL){
  resLFC <- lfcShrink(dds, coef = resultsNames(dds)[2], type = "apeglm") # to remove low count genes (and get rid of the noise)
  png(file.path(dir_to_samples, paste0(name, '_MA-plot.png')), width = 720, height = 720, pointsize = 24)
  plotMA(resLFC, ylim = c(-2,2), alpha = 0.005) # MA-plot for the shrunken log2 fold changes, wo low count genes
  dev.off()
}
get_ma_plot(deseq, dir, name='FA_6_vs_12C')

### Order result data by pvalue and save the table

results_sign <- function(dds_results, 
                         dir_to_samples, 
                         dds, 
                         adj_pv_thr=0.05,
                         lfc = 2,
                         special_name=NULL){
  resOrdered <- dds_results[order(dds_results$pvalue),]
  res_signif <- subset(resOrdered, padj < adj_pv_thr)
  res_signif <- subset(res_signif, abs(log2FoldChange) > lfc)
  write.table(as.data.frame(res_signif), 
            file = file.path(dir_to_samples, paste0(resultsNames(dds)[2], '_',
                                                    adj_pv_thr, '_lfc', lfc, '_', special_name, 
                                                    '.csv')), quote = F, sep = '\t')
}
results_sign(res, dir, deseq, special_name = '1.5C', adj_pv_thr = 0.01, lfc = 2)

### Make a heatmap for genes

heatmap_genes <- function(txi_table, 
                          dds_results, 
                          vst_output, 
                          adj_pv_thr=0.05,
                          dir_to_samples,
                          name=NULL,
                          names_for_rownames=NULL){
  tpm <- txi_table$abundance # take TPMs values
  stpm <- apply(tpm, 1, scale) # standardize data by each row
  stpm <- t(stpm) 
  
  resOrdered <- dds_results[order(dds_results$pvalue),]
  res_signif <- subset(resOrdered, padj < adj_pv_thr)
  
  stpm_sign <- stpm[row.names(stpm) %in% row.names(res_signif),] # take only genes with padj < 0.05 (of log2 fold change)
  stpm_sign <- as.matrix(stpm_sign)
  colnames(stpm_sign) <- paste(vst_output$condition, vst_output$replicate, sep = '_')
  condition_name <- paste(levels(vst_output$condition)[2], 
                          levels(vst_output$condition)[1], sep = '_vs_')
  
  if (!is.null(names_for_rownames)){
  row.names(stpm_sign) <- names_for_rownames
  }
  
  pheatmap(stpm_sign,
           color = plasma(10),
           show_rownames = F,
           border_color = T,
           cluster_cols = T, 
           height = 16,
           width = 7,
           clustering_distance_cols = "correlation",
           clustering_distance_rows = "correlation",
           filename = file.path(dir_to_samples, 
                                paste0(name, '_', condition_name, 
                                       '_heatmap_genes_', adj_pv_thr, '.png')))
  if (dev.cur() != 1) {dev.off()}
}

#toheatmap <- names(chosen_contigs_to_heatmap)
heatmap_genes(txi, res, tr_counts, dir_to_samples = dir, adj_pv_thr=0.05, 
              name='Ecy_1.5Cvs12C') 
#              names_for_rownames = pathway_res_with_kegg$Name)

heatmap_genes_regulateLFC <- function(dds_results,
                                      vst_output,
                                      txi_table,
                                      adj_pv_thr=0.05,
                                      lfc_wanted=3,
                                      dir_to_samples,
                                      name=NULL,
                                      gene_names=T){
  tpm <- txi_table$abundance # take TPMs values
  stpm <- apply(tpm, 1, scale) # standardize data by each row
  stpm <- t(stpm) 
  
  resOrdered <- dds_results[order(dds_results$pvalue),]
  res_signif <- subset(resOrdered, padj < adj_pv_thr)
  res_signif_LFCreg <- subset(res_signif, abs(log2FoldChange) > lfc_wanted)
  stmp_sign_LFCreg <- stpm[row.names(stpm) %in% row.names(res_signif_LFCreg),] 
  stmp_sign_LFCreg <- as.matrix(stmp_sign_LFCreg)
  colnames(stmp_sign_LFCreg) <- rep(c('12°C', '1.5°C'), each = 4)
#  colnames(stmp_sign_LFCreg) <- paste(vst_output$condition, vst_output$replicate, sep = '_')
  condition_name <- paste(levels(vst_output$condition)[2], 
                          levels(vst_output$condition)[1], sep = '_vs_')
  
  pheatmap(stmp_sign_LFCreg,
           color = plasma(10),
           filename = file.path(dir_to_samples, 
                                paste0(name, '_', condition_name, 
                                       '_heatmap_', adj_pv_thr, '_LFC', 
                                       lfc_wanted, '.png')),
           border_color = "black",
           cluster_cols = F,
           clustering_distance_cols = "euclidean",
           clustering_distance_rows = "correlation",
           angle_col = 0,
           fontsize_col = 11,
           fontsize_row = 14,
           show_rownames = gene_names,
           width = 6,
           height = 5)
  if (dev.cur() != 1) {dev.off()}
}

heatmap_genes_regulateLFC(res, tr_counts, txi, dir_to_samples=dir, 
                          adj_pv_thr=0.05, name='Ecy2', lfc_wanted = 2,
                          gene_names = T)

get_summary_DE_LFC <- function(DE_results, vst_output, pvalue_adj=0.05, 
                               lfc_range=1:10){

  resOrdered <- DE_results[order(DE_results$pvalue),]
  res_signif <- subset(resOrdered, padj < pvalue_adj)
  
  s <- data.frame()
  for(lfc in lfc_range) {
    res_lfc <- subset(res_signif, abs(log2FoldChange) > lfc)
    
    Up <- sum(res_lfc$log2FoldChange > lfc)
    Down <- sum(res_lfc$log2FoldChange < -lfc)
    Sum <- sum(Up + Down)
    s <- rbind(s, c(lfc, paste('Up-regulated at ', levels(vst_output$condition)[2]), 
                    Up, round(Up/Sum * 100, 2)), stringsAsFactors=F) 
    s <- rbind(s, c(lfc, paste('Up-regulated at ', levels(vst_output$condition)[1]), 
                    Down, round(Down/Sum * 100, 2)), stringsAsFactors=F)
  }
  colnames(s) <- c('LFC_threshold', 'Regulation', 'Transcripts_number', 'Perc')
  s$Transcripts_number <- as.numeric(s$Transcripts_number)
  s$Perc <- as.numeric(s$Perc)
  s$LFC_threshold <- as.numeric(s$LFC_threshold)
  return(s)
} 

summ_lfc <- get_summary_DE_LFC(res, tr_counts, lfc_range = 1:7)

ggplot(summ_lfc, aes(LFC_threshold, Transcripts_number)) +
  geom_point(aes(color = Regulation)) +
  geom_line(aes(color = Regulation)) +
  scale_color_manual(values = c('#ffa600', '#665191')) +
  xlab('log2 fold change threshold, > |LFC|') +
  ylab('Number of transcripts') +
  theme_light()

ggsave(filename = file.path(dir, 'DA_lfc_range_transcriptsNumber.png'))
ggplot(summ_lfc, aes(LFC_threshold, Perc)) +
  geom_point(aes(color = Regulation)) +
  geom_line(aes(color = Regulation)) + 
  scale_color_manual(values = c('#ffa600', '#665191')) +
  xlab('log2 fold change threshold, > |LFC|') +
  ylab('% of transcripts') +
  theme_light()

ggsave(filename = file.path(dir, 'DA_lfc_range_perc.png'))

### Draw heatmap with gene description 

heatmap_genes_descr <- function(dds_results,
                                      vst_output,
                                      txi_table,
                                      annotation_file,
                                      adj_pv_thr=0.05,
                                      lfc_wanted_down_limit=3,
                                      lfc_wanted_up_limit=Inf,
                                      dir_to_samples,
                                      name=NULL,
                                      gene_names=T, 
                                      row_descr='Description'){
  tpm <- txi_table$abundance # take TPMs values
  stpm <- apply(tpm, 1, scale) # standardize data by each row
  stpm <- t(stpm) 
  
  resOrdered <- dds_results[order(dds_results$pvalue),]
  res_signif <- subset(resOrdered, padj < adj_pv_thr)
  res_signif_LFCreg <- subset(res_signif, abs(log2FoldChange) > lfc_wanted_down_limit
                              & abs(log2FoldChange) <= lfc_wanted_up_limit)
  stmp_sign_LFCreg <- stpm[row.names(stpm) %in% row.names(res_signif_LFCreg),] 
  stmp_sign_LFCreg <- as.matrix(stmp_sign_LFCreg)
  colnames(stmp_sign_LFCreg) <- paste(vst_output$condition, vst_output$replicate, sep = '_')
  
  annotat_file <- read.csv(annotation_file, sep = '\t', 
                              header = F, 
                              col.names = c('ID', 'Gene', 'Description', 
                                            'ORF', 'Domains'))
  ann <- annotat_file[annotat_file$Gene %in% rownames(stmp_sign_LFCreg),]
  ann <- ann[!duplicated(ann$Gene, by = Gene),]
  ann <- ann[order(ann$Gene),]
  ann$Description <- sub('[^ ]* ', '', ann$Description)
  ann$Description <- sub('\\[.*\\]', '', ann$Description)
  ann$Description <- sub('.*hypothetical protein.*', 'uncharacterized protein', ann$Description)
  ann$Description <- sub('.*[uU]ncharacterized protein.*', 'uncharacterized protein', ann$Description)
  ann$Description <- sub('.*predicted protein.*', 'uncharacterized protein', ann$Description)
  ann$Description <- sub('.*unnamed protein product.*', 'uncharacterized protein', ann$Description)
  
  ann$Description <- as.character(ann$Description)
  ann$Domains <- as.character(ann$Domains)
  ann$ID <- as.character(ann$ID)
  
  unchar_proteins <- c()
  ann$Category <- ''
  for (i in 1:nrow(ann)) {
    if (ann$Description[i] == 'uncharacterized protein' && 
        ann$Domains[i] != ''){
        ann$Category[i] <- 'Domain annotation'
    }
    else if (ann$Description[i] == 'uncharacterized protein' && 
             ann$Domains[i] == ''){
      ann$Category[i] <- 'Uncharacterized'
      ann$Description[i] <- ''
      unchar_proteins <- c(unchar_proteins, ann$ID[i])
    }
    else {ann$Category[i] <- 'Blastp annotation'}
  } 
  
  for (i in 1:nrow(ann)){
    if (ann$Description[i] == 'uncharacterized protein' && 
        ann$Domains[i] != ''){
        ann$Description[i] <- ann$Domains[i]
        ann$Description[i] <- paste0(ann$Description[i], ' *')
    }
  }
  
  pheatmap(stmp_sign_LFCreg,
           color = plasma(10),
           filename = file.path(dir_to_samples, 
                                paste0(name, '_heatmap_', adj_pv_thr, 
                                       '_LFC', lfc_wanted_down_limit, 
                                       '_', lfc_wanted_up_limit, '.png')),
           border_color = NA,
           clustering_distance_cols = "correlation",
           clustering_distance_rows = "correlation",
           fontsize_col = 15,
           fontsize_row = 15,
           show_rownames = gene_names,
           width = 19,
           height = 19,
           labels_row = ann$Description)
  if (dev.cur() != 1) {dev.off()}
  
  write.table(unchar_proteins, 
              file = file.path(dir, paste0('uncharacterized proteins_',
                                           adj_pv_thr, '_', 'lfc', 
                                           lfc_wanted_down_limit, '_', 
                                           lfc_wanted_up_limit,
                                           '.txt')),
              row.names = F, quote = F, col.names = F)
#  return(unchar_proteins)
}

while (dev.cur() != 1) { dev.off() }

heatmap_genes_descr(res, tr_counts, txi, annot_file, adj_pv_thr = 0.001, 
                    lfc_wanted_down_limit = 3, lfc_wanted_up_limit = 4, 
                    dir_to_samples = dir, name='description',
                    row_descr = 'Description')
heatmap_genes_descr(res, tr_counts, txi, annot_file, adj_pv_thr = 0.001, 
                    lfc_wanted = 3, dir_to_samples = dir, name='ORF',
                    row_descr = 'ORF')
heatmap_genes_descr(res, tr_counts, txi, annot_file, adj_pv_thr = 0.001, 
                    lfc_wanted = 3, dir_to_samples = dir, name='Domain',
                    row_descr = 'Domains')

### For functional annotated

heatmap_genes_descr_FA <- function(dds_results,
                                vst_output,
                                txi_table,
                                annotation_file,
                                adj_pv_thr=0.05,
                                lfc_wanted_down_limit=3,
                                lfc_wanted_up_limit=Inf,
                                dir_to_samples,
                                name=NULL,
                                gene_names=T, 
                                row_descr='Description'){
  tpm <- txi_table$abundance # take TPMs values
  stpm <- apply(tpm, 1, scale) # standardize data by each row
  stpm <- t(stpm) 
  
  resOrdered <- dds_results[order(dds_results$pvalue),]
  res_signif <- subset(resOrdered, padj < adj_pv_thr)
  res_signif_LFCreg <- subset(res_signif, abs(log2FoldChange) > lfc_wanted_down_limit
                              & abs(log2FoldChange) <= lfc_wanted_up_limit)
  stmp_sign_LFCreg <- stpm[row.names(stpm) %in% row.names(res_signif_LFCreg),] 
  stmp_sign_LFCreg <- as.matrix(stmp_sign_LFCreg)
  colnames(stmp_sign_LFCreg) <- paste(vst_output$condition, vst_output$replicate, sep = '_')
  
  annotat_file <- read.csv(annotation_file, sep = '\t', 
                           header = F, 
                           col.names = c('ID', 'Gene', 'GOs', 'Description'))
  annotat_file <- subset(annotat_file, Gene != '')
  ann <- annotat_file[annotat_file$Gene %in% rownames(stmp_sign_LFCreg),]
  ann <- ann[!duplicated(ann$Gene, by = Gene),]
  ann <- ann[order(ann$Gene),]
#  ann$Description <- sub('[^ ]* ', '', ann$Description)
#  ann$Description <- sub('\\[.*\\]', '', ann$Description)
#  ann$Description <- sub('.*hypothetical protein.*', 'uncharacterized protein', ann$Description)
#  ann$Description <- sub('.*[uU]ncharacterized protein.*', 'uncharacterized protein', ann$Description)
#  ann$Description <- sub('.*predicted protein.*', 'uncharacterized protein', ann$Description)
#  ann$Description <- sub('.*unnamed protein product.*', 'uncharacterized protein', ann$Description)
  
  ann$Description <- as.character(ann$Description)
#  ann$Domains <- as.character(ann$Domains)
  ann$ID <- as.character(ann$ID)
  
  pheatmap(stmp_sign_LFCreg,
           color = plasma(10),
           filename = file.path(dir_to_samples, 
                                paste0(name, '_heatmap_', adj_pv_thr, 
                                       '_LFC', lfc_wanted_down_limit, 
                                       '_', lfc_wanted_up_limit, '.png')),
           border_color = NA,
           clustering_distance_cols = "correlation",
           clustering_distance_rows = "correlation",
           fontsize_col = 15,
           fontsize_row = 15,
           show_rownames = gene_names,
           width = 19,
           height = 10,
           labels_row = ann$Description)
  if (dev.cur() != 1) {dev.off()}
}

heatmap_genes_descr_FA(res, tr_counts, txi, annot_file, adj_pv_thr = 0.05, 
                       lfc_wanted_down_limit = 3, #lfc_wanted_up_limit = 4, 
                       dir_to_samples = dir, name='FA_description',
                       row_descr = 'Description')

library(GOSemSim)
go_1 <- str_split(as.character(annotat_file$GOs[1]), ',')
go_2 <- str_split(annotat_file$GOs[2], ',')
mgoSim(go_1, go_2)

d <- godata('labeglo2/Transcriptomics/annotations/functional_annotation/go-basic.obo', ont="MF", computeIC=FALSE)

### Make heatmap for set of chosen contigs:

heatmap_gene_set <- function(txi_table, 
                          dds_results, 
                          vst_output,
                          gene_set,
                          adj_pv_thr=0.05,
                          lfc=2,
                          dir_to_samples,
                          name=NULL){
  tpm <- txi_table$abundance # take TPMs values
  stpm <- apply(tpm, 1, scale) # standardize data by each row
  stpm <- t(stpm) 
  
  resOrdered <- dds_results[order(dds_results$pvalue),]
  res_signif <- subset(resOrdered, padj < adj_pv_thr)
  res_signif <- subset(res_signif, abs(log2FoldChange) > lfc)
  
  stpm_sign <- stpm[row.names(stpm) %in% row.names(res_signif),] # take only genes with padj < 0.05 (of log2 fold change)
  stpm_sign <- as.matrix(stpm_sign)
  colnames(stpm_sign) <- paste(vst_output$condition, vst_output$replicate, sep = '_')

#  stpm_sign <- stpm_sign[match(gene_set, row.names(stpm_sign)),]
  stpm_sign <- stpm_sign[row.names(stpm_sign) %in% gene_set,]
  colnames(stpm_sign) <- rep(c('12°C', '1.5°C'), each = 4)
#  print(stpm_sign)

  pheatmap(stpm_sign,
#           color = rev(colorRampPalette(brewer.pal(8, "PiYG"))(25)),
#           color = cividis(10),
#           color = inferno(10),
#           color = magma(10),
           color = viridis(10),
           show_rownames = T,
           border_color = T,
           clustering_distance_cols = "euclidean",
           cluster_cols = F,
           cluster_rows = T,
           clustering_distance_rows = "euclidean",
           fontsize_row = 16,
           fontsize_col = 14,
           angle_col = 0,
           width = 7,
           height = 7,
           filename = file.path(dir_to_samples, paste0(name, '_heatmap_genes_', adj_pv_thr, 'lfc', lfc, '.png')))
  if (dev.cur() != 1) {dev.off()}
}

toheatmap <- names(chosen_contigs_to_heatmap_stat)
toheatmap <- unique(toheatmap)
#toheatmap <- c(toheatmap, 'Vg')
heatmap_gene_set(txi, res, tr_counts, toheatmap, dir_to_samples = dir, 
                 adj_pv_thr = 0.05, lfc = 0, name='gene_expression')
# RNA helicases:
toheatmap <- row.names(res[grepl('DDX|DHX|Cirbp|cirp|dbp3|ded1|eIF4A|csda|Dead|Dicer|Sub2p|UAP56|Mrh4|Mss116p|Mtt1p|ighmbp2|SETX|IBP160|MOV10|HELC1|Slh1p|SKIL2|Mtr4p', row.names(res), 
                                      ignore.case = T),])
toheatmap <- row.names(res[grepl('DDX', row.names(res), 
                                 ignore.case = T),])
# Proteins with cold-shock domain:
toheatmap <- row.names(res[grepl('UNR|CSDE1|YB-1|DBPB|YBX1|FRGY1|Yps|DjY1|FRGY2|contrin|MSY2|DBPC|DBPA|ZONAB|CSDA|Lin-28|eIF1A|eIF2|eIF5A|RRP44|Dis3|RRP40|RPR22|GRP2|AtCSP2', row.names(res), 
                                 ignore.case = T),])
# Desaturases:
toheatmap <- row.names(res[grepl('^(FADS.*|fat1|fat2|fat-.*|ERG3|DES1|SCD.*|desA|degs.*|desat1|ole1|fad8|fadB|pds1|efd1)$', row.names(res), 
                                 ignore.case = T),])

toheatmap <- row.names(res[grepl('AFP', row.names(res), 
                                 ignore.case = T),])
toheatmap <- row.names(res[grepl('DNAH|DNAL', row.names(res), 
                                 ignore.case = T),])
toheatmap <- row.names(res[grepl('^(Nrf1|NFE2L1|NRF2|DDI2|DDI1|VSM1|RELA|IKBKG|NFKB1|TIP60|Rpn4|Rpn.*|PSMD13|Rpt.*|NAS7|SKN1|YAP1|mtor|hsf1|pdr1|pdr3)$', row.names(res), 
                                 ignore.case = T),])
toheatmap <- row.names(res[grepl('K09040|K06693|K11831|K05638', row.names(res), 
                                 ignore.case = T),])
toheatmap <- row.names(res[grepl('^(UBB|UBC|UBA52|RPS27A)$', row.names(res), 
                                 ignore.case = T),])

### DE with miltiple conditions
library(multiDE)

count_data <- txi$counts
norm_count_data <- normalization(count_data, 'median')

condition <- c(rep('1.5C', 4), rep(c('12C'), 4), rep('6C', 3))

count_data_disp <- dispersion(M=norm_count_data, condition=condition, 
                              matched=F, method='DESeq2')
multi_deseq <- multiDE(count_data_disp)

tr_counts <- vst(count_data_disp$count)

ss_dist_heatmap <- function(vst_output, dir_to_samples, name=NULL){
  vst_m <- tr_counts_batch # to get matrix with transformed counts from vst() object
  vst_m <- vst_m[,colnames(vst_m) != '24C_1']
  sample_dist <- dist(t(vst_m), method = "euclidean") # t - transform data to samples-genes matrix (from gene-samples)
  sampleDistMatrix <- as.matrix(sample_dist)
  
  rownames(sampleDistMatrix) <- column_names[column_names != '24C_1']
  colnames(sampleDistMatrix) <- NULL
#  sampleDistMatrix <- sampleDistMatrix[row.names(sampleDistMatrix) != '24C_1',]
#  condition_name <- paste(levels(vst_output$condition)[2], 
#                          levels(vst_output$condition)[1], sep = '_vs_')
  
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sample_dist,
           clustering_distance_cols = sample_dist,
           color = plasma(10),
           filename = file.path(dir, paste0('1.5_12_6_24C', '_', '_heatmap_SampleDIST.png')),
           width = 7,
           height = 5)
  if (dev.cur() != 1) {dev.off()}
}

multiDE_pca <- prcomp(t(tr_counts))
df_out <- as.data.frame(multiDE_pca$x)
df_out$condition <- condition
ggplot(df_out, aes(x=PC1, y=PC2, color=condition))+
  geom_point()+
  theme_light()

heatmap_genes <- function(txi_table, 
                          dds_results, 
                          vst_output, 
                          adj_pv_thr=0.05,
                          dir_to_samples,
                          name=NULL){
  tpm <- txi$abundance # take TPMs values
  stpm <- apply(tpm, 1, scale) # standardize data by each row
  stpm <- t(stpm) 
  
#  resOrdered <- multi_deseq[order(multi_deseq$p.value),]
  multi_de_pvalue <- data.frame(p.value = multi_deseq$p.value)
  row.names(multi_de_pvalue) <- row.names(count_data)
  res_signif <- subset(multi_de_pvalue, p.value < 0.0001)
  
  stpm_sign <- stpm[row.names(stpm) %in% row.names(res_signif),] # take only genes with padj < 0.05 (of log2 fold change)
  stpm_sign <- as.matrix(stpm_sign)
  colnames(stpm_sign) <- condition
  #condition_name <- paste(levels(vst_output$condition)[2], 
  #                        levels(vst_output$condition)[1], sep = '_vs_')
  
  pheatmap(stpm_sign,
           color = plasma(10),
           show_rownames = F,
           border_color = NA,
           cluster_cols = T,
           clustering_distance_cols = "correlation",
           clustering_distance_rows = "correlation")
           #filename = file.path(dir_to_samples, 
          #                      paste0(name, '_', condition_name, 
           #                            '_heatmap_genes_', adj_pv_thr, '.png')))
  if (dev.cur() != 1) {dev.off()}
}

counts_batch <- txi$counts
column_names <- paste(rep(c('12C', '1.5C', '6C'), each = 4), rep(1:4, 3), sep = '_')
column_names <- column_names[column_names != '24C_4']
colnames(counts_batch) <- column_names

counts_batch_int <- round(counts_batch)

tr_counts_batch <- vst(counts_batch_int)

pca <- prcomp(t(tr_counts_batch))
#pca$x
df_out <- as.data.frame(pca$x)
df_out$Condition <- c(rep(c('12°C', '1.5°C'), each = 4), rep('6°C', 4))
#df_out <- df_out[row.names(df_out) != '24C_1',]
ggplot(df_out, aes(x=PC1, y=PC2, color=Condition))+
  geom_point()+
  theme_light()

ggsave(file.path(dir, 'pca_1.5_12_6.png'))

