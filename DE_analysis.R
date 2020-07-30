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
#dir <- 'labeglo2/Transcriptomics/quantification/AEA_GlaReads_EcyAssembly/GlaReads_EcyAssembly'
#dir <- 'labeglo2/Transcriptomics/quantification/GlaReads_EcyAssembly_all_samples/'

#files_names <- c('350', '351',  '352', '354', '632', '633',  '634',  '636',
#                 'Eve_6C_rep1', 'Eve_6C_rep3', 'Eve_6C_rep4',
#                 'Eve_24C_rep1','Eve_24C_rep2','Eve_24C_rep3', 'Eve_24C_rep4') # specify the names of folders with quant.sf data 

files_names <- c('311', '312', '313', '315', '594', '595', '596', '598') # specify the names of folders with quant.sf data, Gla
files_names <- c('129', '130', '131', '133', '560', '561', '562', '563') # Ecy
files_names <- c('350', '351',  '352', '354', '632', '633',  '634',  '636') # Eve
#files_names <- c('594_glaA', '595_glaA', '596_glaA', '598_glaA', 
#                 '560_glaA', '561_glaA', '562_glaA', '563_glaA') # Gla 1.5 and Ecy 1.5
#files_names <- c('594_ecyA', '595_ecyA', '596_ecyA', '598_ecyA', 
#                 '560_ecyA', '561_ecyA', '562_ecyA', '563_ecyA')
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
#txi <- txi_creation(dir, 'contigs_DAnnotation.txt', files)
#txi <- txi_creation(dir, 'Gla_kegg_pred_cds_annotation_wo_p_and_dupl.txt', files)
#txi <- txi_creation(dir, 'kegg_annot_Gla_assembly.txt', files)
#txi <- txi_creation(dir, 'pred_peptides_functional_annot_wo_p_and_Dupl.txt', files)
#txi$counts <- txi$counts[-1,]

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
#sample_description <- table_with_sample_descr(8, c('Gla', 'Ecy'))
#sample_description$condition <- factor(sample_description$condition, 
#                                       levels = c('Ecy', 'Gla'))
#sample_description$condition <- factor(sample_description$condition, 
#                                       levels = c('1.5C', '12C', '6C', '24C'))
#sample_description$condition <- ifelse(sample_description$condition == '1.5C', 
#                                       '1.5C', '6or12C')
#sample_description$replicate <- c(paste0('12C_', 1:4), paste0('1.5C_', 1:4), 
#                                  paste0('6C_', 1:3), paste0('24C_', 1:4))
#sample_description$condition <- factor(sample_description$condition, 
#                                       levels = c('6or12C','1.5C'))
#sample_description <- sample_description[!rownames(sample_description) == 'sample10',]
#row.names(sample_description) <- paste0('sample', 1:7)

### Run differential expression analysis based on the negative binominal distribution

deseq <- DESeqDataSetFromTximport(txi, colData = sample_description, 
                                      design = ~ condition)
deseq$condition

deseq <- DESeq(deseq)
res <- results(deseq, name = resultsNames(deseq)[2]) # results() extracts a results table with log2 fold changes, p values and adjusted p values  
#res[grepl('XP_015433778.1|RNA01364.1|RNA33404.1', row.names(res)),] # select set of genes

res[grepl('CYB5R3|DIA1|Cyb5r3|MCR1|YKL150W|YKL605|CBR1|CBR|CYB5R', row.names(res), ignore.case = T),] # NADH-cytochrome b5 reductase -> one of the component of aerobic desaturation (PUFA synthesis)
res[grepl('CYB5A|CYB5B|CYB5m|Omb5|CB5-A|N1949|YNL111C|Cyt-b5|CG2140', row.names(res), ignore.case = T),] # cytochrome b5 -> one of the component of aerobic desaturation (PUFA synthesis)
res[grepl('FADS|SCD|des6|Delta6FAD|FAD2|desat|OLE1|fat-6', row.names(res), ignore.case = T),] # acyl-coa desaturases -> one of the component of aerobic desaturation (PUFA synthesis)

res_sign <- subset(res, padj < 0.05)
res_ordered <- res_sign[order(res_sign$log2FoldChange),]
res_sign_lfc <- subset(res_sign, log2FoldChange > 2)
rownames(res_sign_lfc)

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

