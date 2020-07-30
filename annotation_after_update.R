library(stringr)
library(plyr)
library(ggplot2)
library(remotes)
library(ggsci)
library(scales)

# upload annotation file obtained by DIAMOND run (tabular format)
ecy_annot <- read.csv("labeglo2/Transcriptomics/annotations/EcySSv2_annotation.diamond.csv", 
                      sep = "\t", 
                      col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", 
                      "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore",
                      "staxids", "sscinames", "stitle", "salltitles"))
# number of annotated contigs
length(ecy_annot$qseqid) # 165409 -> not so much -> maybe to run Diamond with --very-sensitive mode?

# save taxids in new .csv to get high taxonomic levels 
tax_ids_Ecy <- data.frame(staxids = ecy_annot$staxids)
tax_ids_Ecy <- str_split_fixed(tax_ids_Ecy$staxids, ";", 2) # some of the ids have double id separated by ";" (the second one is strain)
tax_ids_Ecy <- data.frame(taxids = tax_ids_Ecy[,-2]) # get rid of the everything after ";"
tax_ids_Ecy$taxids <- as.numeric(as.vector(tax_ids_Ecy$taxids))

tax_ids_Ecy_wo_NA <- na.omit(tax_ids_Ecy) # remove NA 
length(tax_ids_Ecy_wo_NA$taxids) # number of contigs with known taxonomic ids -> 162535

write.table(tax_ids_Ecy_wo_NA, "labeglo2/Transcriptomics/annotations/taxids_Ecy.csv", row.names = F, col.names = F)
# then for taxonomy identification it needs to extract Lineage info from Taxonomy DB
# (I wrote python script to fetch from Taxonomy database the Lineage) 
# -> ~/labeglo2/Transcriptomics/annotations/taxonomy_extract.py

# upload lineages:
lineages_Ecy <- read.csv("labeglo2/Transcriptomics/annotations/taxonomy_issues/lineage_EcySSv2_scinames.csv", header = F)
lineages_Ecy <- str_split_fixed(lineages_Ecy$V1, ";", n = Inf) # separate every group into columns
lineages_Ecy <- as.data.frame(lineages_Ecy) 
count(lineages_Ecy$V1) # to know number of contigs in different groups of V1
lineages_Ecy_cleaned <- subset(lineages_Ecy, V1 != 'unclassified sequences')
lineages_Ecy_cleaned <- subset(lineages_Ecy_cleaned, V1 != 'other sequences') # now lineages_Ecy_cleaned has only Cellular organisms and Viruses

superkingdom_Ecy <- count(lineages_Ecy_cleaned$V1) # it is only 282 viruses vs 162242 of cellular organisms (0.17%)
lineages_Ecy_cleaned_wo_viruses <- subset(lineages_Ecy_cleaned, V1 == 'cellular organisms') # get rid of viruses
# The abundance of different kingdoms in assembly (and sorting the list of them):
kingdom_Ecy <- count(lineages_Ecy_cleaned_wo_viruses$V2) 
kingdom_Ecy <- kingdom_Ecy[order(-kingdom_Ecy$freq),]
# The abundance of different clades of 2nd level:
classes_Ecy <- count(lineages_Ecy_cleaned_wo_viruses$V3)
big_classes_Ecy <- subset(classes_Ecy, freq > 8000) # to take only huge groups
big_classes_Ecy <- big_classes_Ecy[order(-big_classes_Ecy$freq),] # to sort it

# to make dataframe with group 'Others' for 2nd level which are the sum of smaller groups
classes_Ecy_wothers <- data.frame(Classes = c(as.vector(big_classes_Ecy$x), "Others2"),
                                  Frequencies = c(big_classes_Ecy$freq, sum(subset(classes_Ecy, freq < 8000)[2])))
#
level3_Ecy <- count(lineages_Ecy_cleaned_wo_viruses$V4)
bid_level3_groups_Ecy <- subset(level3_Ecy, freq > 4000)
big_level3_groups_Ecy <- bid_level3_groups_Ecy[order(-bid_level3_groups_Ecy$freq),]
#
level3_Ecy_wothers <- data.frame(Classes = c(as.vector(big_level3_groups_Ecy$x), "Others3"),
                                 Frequencies = c(big_level3_groups_Ecy$freq, sum(subset(level3_Ecy, freq < 4000)[2])))
level3_Ecy_wothers <- level3_Ecy_wothers[order(-level3_Ecy_wothers$Frequencies),]
# to make dataframe for ggplot
taxon_freq_Ecy <- data.frame(Clades = c(rep('1st level', 3), rep('2nd level', 5), rep('3rd level', 6)),
                             Taxon = factor(c(as.vector(kingdom_Ecy$x), as.vector(classes_Ecy_wothers$Classes), as.vector(level3_Ecy_wothers$Classes)), 
                                            levels = c('Eukaryota', 'Bacteria', 'Archaea', 'Opisthokonta', 'Sar', 'Viridiplantae', 'Proteobacteria', 'Others2', 'Metazoa', 'Alveolata', 'Others3', 'Streptophyta', 'Alphaproteobacteria', 'Fungi')),
                             Frequency = c(kingdom_Ecy$freq, classes_Ecy_wothers$Frequencies, level3_Ecy_wothers$Frequencies))

ggplot(taxon_freq_Ecy, aes(x = Clades, y = Frequency)) +
  geom_col(aes(fill = Taxon)) +
  geom_text(aes(group = Frequency, label = Taxon), 
            position = position_stack(reverse = TRUE, vjust = .5))

# Metazoa from annotation of non-filtered contigs 

lineages_Metazoa_ecy <- subset(lineages_Ecy_cleaned_wo_viruses, V4 == 'Metazoa')
count(lineages_Metazoa_ecy$V5)
lineages_Metazoa_ecy_wo_Porifera <- subset(lineages_Metazoa_ecy, V5 == 'Eumetazoa')
count(lineages_Metazoa_ecy_wo_Porifera$V6)
lineages_Bilateria_ecy <- subset(lineages_Metazoa_ecy_wo_Porifera, V6 == 'Bilateria')
count(lineages_Bilateria_ecy$V7)
lineages_Protostomia_ecy <- subset(lineages_Bilateria_ecy, V7 == 'Protostomia')
count(lineages_Protostomia_ecy$V8)

# I obtained annotation of long ORFs (TransDecoder) of filtered contigs for Ecy assembly (2020-03-12)

ecy_orf_annot <- read.csv('labeglo2/Transcriptomics/annotations/Ecy_longORFs_woNAandDupl.diamond.csv', 
                          sep = '\t',
                          col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", 
                                        "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore",
                                        "staxids", "sscinames", "sskingdoms", "skingdoms", "sphylums", "stitle", "salltitles"))
length(ecy_orf_annot$qseqid) # 16099
#ecy_orf_annot <- subset(ecy_orf_annot, endsWith(as.character(ecy_orf_annot$qseqid), '.p1')) # because of the output of TransDecoder.LongOrfs, do not need anymore -> fixed it in command line 
#length(ecy_orf_annot$qseqid) # 16160

#which(ecy_orf_annot$sskingdoms == "0;Eukaryota")
count(ecy_orf_annot$sskingdoms)
subset(ecy_orf_annot, skingdoms == "0;Metazoa")
subset(ecy_orf_annot, sskingdoms == "0;Eukaryota") # for 14070th record we have mistake in record: '0;Eukaryota 0;Metazoa 0;Arthropoda'
ecy_orf_annot[which(ecy_orf_annot$sskingdoms == "0;Eukaryota"),]$skingdoms <- "Metazoa"
ecy_orf_annot[which(ecy_orf_annot$sskingdoms == "0;Eukaryota"),]$sphylums <- "Arthropoda"
ecy_orf_annot[which(ecy_orf_annot$sskingdoms == "0;Eukaryota"),]$sskingdoms <- "Eukaryota" # fix mistake

subset(ecy_orf_annot, sskingdoms == "0") # Expression vector pGDR11-9N -> that's why it is better to retrieve not only one targets from database
subset(ecy_orf_annot, sskingdoms == "Bacteria;Eukaryota")[,17] # For now I remove them, I don't know why it happened
ecy_orf_annot <- subset(ecy_orf_annot, sskingdoms != "Bacteria;Eukaryota")
ecy_orf_annot <- subset(ecy_orf_annot, sskingdoms != "0")
#subset(ecy_orf_annot, sskingdoms == "N/A")[,1:2] # some of them hypothetical proteins, but some is not; a lot of proteins from Portunus trituberculatus, maybe it makes sense to recover taxonomy which here is lost for no reasons !!! I fixed it with python script

taxon_freq_ORF_king <- count(ecy_orf_annot$skingdoms) # 14,342 of Metazoa
# Prepare data frame for plot
taxon_freq_ORF_king$x <- as.vector(taxon_freq_ORF_king$x)
taxon_freq_ORF_king$x[1] <- "Bacteria and Protozoa"
# unicel <- subset(ecy_orf_annot, skingdoms == '0')

taxon_freq_Ecy_orfs <- data.frame(Species = rep('Ecy', 4),
                             Taxon = factor(c(as.vector(taxon_freq_ORF_king$x)), 
                                            levels = c('Metazoa', 'Bacteria and Protozoa', 'Viridiplantae', 'Fungi')),
                             Frequency = c(taxon_freq_ORF_king$freq))

mypal <- pal_futurama()(4)
show_col(mypal)
mypal

ggplot(taxon_freq_Ecy_orfs, aes(x = Species, y = Frequency)) +
  geom_col(aes(fill = Taxon)) +
#  geom_text(aes(group = Frequency, label = Taxon), 
#            position = position_stack(reverse = TRUE, vjust = .5)) +
  coord_flip() +
  ylab("Number of long ORFs") +
  scale_fill_manual(values = c("#FF6F00FF", "#C71000FF", "#008EA0FF", "#8A4198FF")) +
  theme_bw()

# Find out the distribution of ORFs from Metazoa 
ecy_orf_metazoa <- subset(ecy_orf_annot, skingdoms == 'Metazoa')
count(ecy_orf_metazoa$sphylums)
metazoa_Ecy_orfs <- data.frame(Phylum = count(ecy_orf_metazoa$sphylums)$x,
                               Frequency = count(ecy_orf_metazoa$sphylums)$freq)

metazoa_Ecy_orfs <- metazoa_Ecy_orfs[order(metazoa_Ecy_orfs$Frequency),]
metazoa_Ecy_orfs$Phylum <- factor(metazoa_Ecy_orfs$Phylum, 
                                  levels = metazoa_Ecy_orfs$Phylum)

ggplot(metazoa_Ecy_orfs, aes(x = Frequency, y = Phylum)) +
  geom_col(fill = "#8A4198FF") +
  xlab('Number of long OFRs') +
  theme_light()

# Distribution of ORFs from Arthropoda
ecy_orf_arthropoda <- subset(ecy_orf_metazoa, sphylums == 'Arthropoda')
ecy_orf_arthropoda_freq <- count(ecy_orf_arthropoda$sscinames) # 13,348 records
more_10_arthropoda_ecy <- subset(ecy_orf_arthropoda_freq, freq > 10) # 12,852 records

more_10_arthropoda_ecy <- data.frame(Species = more_10_arthropoda_ecy$x,
                                     Frequency = more_10_arthropoda_ecy$freq)
more_10_arthropoda_ecy <- more_10_arthropoda_ecy[order(more_10_arthropoda_ecy$Frequency),]
more_10_arthropoda_ecy$Species <- factor(more_10_arthropoda_ecy$Species,
                                         levels = more_10_arthropoda_ecy$Species)

mypal2 <- pal_cosmic()(4)
show_col(mypal2)
mypal2

ggplot(more_10_arthropoda_ecy, aes(x = Frequency, y = Species)) +
  geom_col(fill = '#CF4E9CFF') +
  xlab('Number of long ORFs') +
  theme_light()

# let's take only Arthropoda for the transcripts quantification

length(ecy_orf_arthropoda$qseqid) # 13,348 records

write.table(ecy_orf_arthropoda, file = 'labeglo2/Transcriptomics/annotations/arthropoda_orf_ecy.csv', sep = '\t', quote = F, col.names = F)

# And if take Arthropoda, Mollusca and Chordata:

ecy_orf_Arthr_Chordata_Mollusca <- subset(ecy_orf_metazoa, 
                                                 sphylums == 'Arthropoda' |
                                                   sphylums == 'Chordata' |
                                                   sphylums == 'Mollusca')
count(ecy_orf_Arthr_Chordata_Mollusca$sphylums)
write.table(ecy_orf_Arthr_Chordata_Mollusca, file = 'labeglo2/Transcriptomics/annotations/arthr_chord_mollus_orf_ecy.csv', sep = '\t', quote = F, col.names = F)

# I obtained annotation of filtered contigs for Ecy assembly (2020-03-13)
# And I recover taxonomy for N/A fields
#'labeglo2/Transcriptomics/annotations/Ecy_filtered_contigs_woNA.diamond.csv'

annotation_data <- function(path_to_csv) 
{read.csv(path_to_csv, 
          sep = '\t',
          col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", 
                        "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore",
                        "staxids", "sscinames", "sskingdoms", "skingdoms", "sphylums", "stitle", "salltitles"))
}

species_filt_annot <- annotation_data('labeglo2/Transcriptomics/annotations/Ecy_filtered_contigs_woNA.diamond.csv') # Ecy
species_filt_annot <- annotation_data('labeglo2/Transcriptomics/annotations/Eve_filtered_contigs_woNAandDupl.diamond.csv') # Eve
species_filt_annot <- annotation_data('labeglo2/Transcriptomics/annotations/Gla_filtered_contigs_woNAandDupl.diamond.csv') # Gla
species_filt_annot <- annotation_data('labeglo2/Transcriptomics/annotations/EveNEB_muscle_filtered_woNAandDupl.diamond.csv') # EveNEB, muscle

count(species_filt_annot$sskingdoms)
species_filt_annot <- subset(species_filt_annot, sskingdoms != "0;Eukaryota")
species_filt_annot <- subset(species_filt_annot, sskingdoms != "Bacteria;Archaea")
species_filt_annot <- subset(species_filt_annot, sskingdoms != "Bacteria;Eukaryota")
species_filt_annot <- subset(species_filt_annot, sskingdoms != "0")
count(species_filt_annot$sskingdoms)

count(species_filt_annot$skingdoms)
species_filt_annot <- subset(species_filt_annot, skingdoms != "0;Metazoa")
species_filt_annot <- subset(species_filt_annot, skingdoms != "Fungi;Metazoa")
count(species_filt_annot$skingdoms)

taxon_freq_FiltContigs_king <- count(species_filt_annot$skingdoms)

taxon_freq_FiltContigs_king$x <- as.vector(taxon_freq_FiltContigs_king$x)
#taxon_freq_FiltContigs_king$x[1] <- "Bacteria and Protozoa"
unicel <- subset(species_filt_annot, skingdoms == '0')
unicel_count <- count(unicel$sskingdoms)
unicel_count$x <- sub('Eukaryota', 'Unicellular eukaryota', unicel_count$x)

taxon_freq_species_FiltContigs <- data.frame(Species = rep('Eve, muscle', 7),
                                             Taxon = factor(c(as.vector(taxon_freq_FiltContigs_king$x)[-1], as.vector(unicel_count$x)), 
                                                            levels = c('Metazoa', 'Unicellular eukaryota', 'Bacteria', 
                                                                       'Viridiplantae', 'Fungi', 'Archaea', 'Viruses')),
                                             Frequency = c(taxon_freq_FiltContigs_king$freq[-1], unicel_count$freq))

#taxon_freq_species_FiltContigs <- data.frame(Species = rep('Ecy', 4),
#                                  Taxon = factor(c(as.vector(taxon_freq_FiltContigs_king$x)), 
#                                                 levels = c('Metazoa', 'Bacteria and Protozoa', 'Viridiplantae', 'Fungi')),
#                                  Frequency = c(taxon_freq_FiltContigs_king$freq))

cbPalette <- c("#E69F00", "#56B4E9", "#F0E442", "#009E73", "#CC79A7", "#0072B2", "#D55E00")
ggplot(taxon_freq_species_FiltContigs, aes(x = Species, y = Frequency)) +
  geom_col(aes(fill = Taxon)) +
  #  geom_text(aes(group = Frequency, label = Taxon), 
  #            position = position_stack(reverse = TRUE, vjust = .5)) +
  coord_flip() +
  ylab("Number of contigs") +
  scale_fill_manual(values = cbPalette) +
  theme_bw()

## Different colours
ggplot(taxon_freq_species_FiltContigs, aes(x = Species, y = Frequency)) +
  geom_col(aes(fill = Taxon)) +
  #  geom_text(aes(group = Frequency, label = Taxon), 
  #            position = position_stack(reverse = TRUE, vjust = .5)) +
  coord_flip() +
  ylab("Number of contigs") +
  scale_fill_manual(values = c("#FF6F00FF", "#C71000FF", "#008EA0FF", "#8A4198FF")) +
  theme_bw()

# Find out the distribution of filtered contigs from Metazoa 
species_filtContigs_metazoa <- subset(species_filt_annot, skingdoms == 'Metazoa')
count(species_filtContigs_metazoa$sphylums)
species_filtContigs_metazoa[which(species_filtContigs_metazoa$sphylums != 'Arthropoda;Chordata;Hemichordata'),]
species_filtContigs_metazoa <- subset(species_filtContigs_metazoa, sphylums != 'Arthropoda;Chordata;Hemichordata')
species_filtContigs_metazoa <- subset(species_filtContigs_metazoa, sphylums != 'Platyhelminthes;Mollusca')

count(species_filtContigs_metazoa$sphylums)

metazoa_species_filtContigs <- data.frame(Phylum = count(species_filtContigs_metazoa$sphylums)$x,
                               Frequency = count(species_filtContigs_metazoa$sphylums)$freq)

metazoa_species_filtContigs <- metazoa_species_filtContigs[order(metazoa_species_filtContigs$Frequency),]
metazoa_species_filtContigs$Phylum <- factor(metazoa_species_filtContigs$Phylum, 
                                  levels = metazoa_species_filtContigs$Phylum)

ggplot(metazoa_species_filtContigs, aes(x = Frequency, y = Phylum)) +
  geom_col(fill = "#8A4198FF") +
  xlab('Number of contigs') +
  theme_light()

##  Save csv with contigs assigned to Metazoa:
write.table(species_filtContigs_metazoa, file = 'labeglo2/Transcriptomics/annotations/metazoa_contigs_EVE_MUSCLE.csv', sep = '\t', quote = F, col.names = F)

# Distribution of contigs assigned by Arthropoda
species_FiltContig_arthropoda <- subset(species_filtContigs_metazoa, sphylums == 'Arthropoda')
species_FiltContig_arthropoda_freq <- count(species_FiltContig_arthropoda$sscinames) # 31,202 records for Ecy; 36,542 for Eve; 37,365 for Gla; 21,949 for Eve Muscle
more_10_arthropoda_FiltContig_species <- subset(species_FiltContig_arthropoda_freq, freq > 20) # 30,493 records for Ecy (>10); 34,886 for Eve (> 20); 36,100 for Gla (>20);20,993 for Eve Muscle

species_FiltContig_arthropoda_freq <- data.frame(Species = species_FiltContig_arthropoda_freq$x,
                                     Frequency = species_FiltContig_arthropoda_freq$freq)
species_FiltContig_arthropoda_freq <- species_FiltContig_arthropoda_freq[order(species_FiltContig_arthropoda_freq$Frequency),]
species_FiltContig_arthropoda_freq$Species <- factor(species_FiltContig_arthropoda_freq$Species,
                                         levels = species_FiltContig_arthropoda_freq$Species)
species_FiltContig_arthropoda_freq$Species1 <- factor(
  ifelse(species_FiltContig_arthropoda_freq$Frequency >= 50,
         as.character(species_FiltContig_arthropoda_freq$Species),
         'Other'),
  levels=c('Other', levels(species_FiltContig_arthropoda_freq$Species)))

ggplot(species_FiltContig_arthropoda_freq) +
  geom_col(aes(x = Frequency, y = Species1),
    fill = '#CF4E9CFF') +
  xlab('Number of contigs') +
  ylab('Species') +
  theme_light()

species_FiltContig_mollusca <- subset(species_filtContigs_metazoa, sphylums == 'Mollusca')
count(species_FiltContig_mollusca$sscinames)

species_FiltContig_chordata <- subset(species_filtContigs_metazoa, sphylums == 'Chordata')
View(count(species_FiltContig_chordata$sscinames))

## let's take only Arthropoda, Chordata and Mollusca for the transcripts quantification
species_FiltContig_Arthr_Chotdata_Mollusca <- subset(species_filtContigs_metazoa, 
                                                 sphylums == 'Arthropoda' |
                                                   sphylums == 'Chordata' |
                                                   sphylums == 'Mollusca')
count(species_FiltContig_Arthr_Chotdata_Mollusca$sphylums)
write.table(species_FiltContig_Arthr_Chotdata_Mollusca, file = 'labeglo2/Transcriptomics/annotations/arthr_chord_mollus_contigs_Gla.csv', sep = '\t', quote = F, col.names = F)

## Plot distribution for all three species
library(formattable)

#taxon_freq_species_FiltContigs_Gla <- data.frame(Species = rep('Gla', 4),
#                                             Taxon = factor(c(as.vector(taxon_freq_FiltContigs_king$x)), 
#                                                            levels = c('Metazoa', 'Bacteria and Protozoa', 'Viridiplantae', 'Fungi')),
#                                             Frequency = c(taxon_freq_FiltContigs_king$freq))
taxon_freq_species_FiltContigs_Gla <- data.frame(Species = rep('Gla', 7),
                                             Taxon = factor(c(as.vector(taxon_freq_FiltContigs_king$x)[-1], as.vector(unicel_count$x)), 
                                                            levels = c('Metazoa', 'Unicellular eukaryota', 'Bacteria', 
                                                                       'Viridiplantae', 'Fungi', 'Archaea', 'Viruses')),
                                             Frequency = c(taxon_freq_FiltContigs_king$freq[-1], unicel_count$freq))
taxon_freq_species_FiltContigs_Gla$Perc <- round(sapply(taxon_freq_species_FiltContigs_Gla$Frequency, 
                                                  function(x) x/sum(taxon_freq_species_FiltContigs_Gla$Frequency) * 100), 1)
#save(taxon_freq_species_FiltContigs_Gla, file = 'labeglo2/Transcriptomics/annotations/taxonomy_issues/Gla_taxons_counts.RData')
load(file = 'labeglo2/Transcriptomics/annotations/taxonomy_issues/Gla_taxons_counts.RData')

#taxon_freq_species_FiltContigs_Ecy <- data.frame(Species = rep('Ecy', 4),
#                                                 Taxon = factor(c(as.vector(taxon_freq_FiltContigs_king$x)), 
#                                                                levels = c('Metazoa', 'Bacteria and Protozoa', 'Viridiplantae', 'Fungi')),
#                                                 Frequency = c(taxon_freq_FiltContigs_king$freq))
taxon_freq_species_FiltContigs_Ecy <- data.frame(Species = rep('Ecy', 7),
                                             Taxon = factor(c(as.vector(taxon_freq_FiltContigs_king$x)[-1], as.vector(unicel_count$x)), 
                                                            levels = c('Metazoa', 'Unicellular eukaryota', 'Bacteria', 
                                                                       'Viridiplantae', 'Fungi', 'Archaea', 'Viruses')),
                                             Frequency = c(taxon_freq_FiltContigs_king$freq[-1], unicel_count$freq))
taxon_freq_species_FiltContigs_Ecy$Perc <- round(sapply(taxon_freq_species_FiltContigs_Ecy$Frequency, 
                                                        function(x) x/sum(taxon_freq_species_FiltContigs_Ecy$Frequency) * 100), 1)
#save(taxon_freq_species_FiltContigs_Ecy, file = 'labeglo2/Transcriptomics/annotations/taxonomy_issues/Ecy_taxons_counts.RData')
load(file = 'labeglo2/Transcriptomics/annotations/taxonomy_issues/Ecy_taxons_counts.RData')

#taxon_freq_species_FiltContigs_Eve <- data.frame(Species = rep('Eve', 4),
#                                                 Taxon = factor(c(as.vector(taxon_freq_FiltContigs_king$x)), 
#                                                                levels = c('Metazoa', 'Bacteria and Protozoa', 'Viridiplantae', 'Fungi')),
#                                                 Frequency = c(taxon_freq_FiltContigs_king$freq))
taxon_freq_species_FiltContigs_Eve <- data.frame(Species = rep('Eve', 7),
                                             Taxon = factor(c(as.vector(taxon_freq_FiltContigs_king$x)[-1], as.vector(unicel_count$x)), 
                                                            levels = c('Metazoa', 'Unicellular eukaryota', 'Bacteria', 
                                                                       'Viridiplantae', 'Fungi', 'Archaea', 'Viruses')),
                                             Frequency = c(taxon_freq_FiltContigs_king$freq[-1], unicel_count$freq))
taxon_freq_species_FiltContigs_Eve$Perc <- round(sapply(taxon_freq_species_FiltContigs_Eve$Frequency, 
                                                  function(x) x/sum(taxon_freq_species_FiltContigs_Eve$Frequency) * 100), 1)
#save(taxon_freq_species_FiltContigs_Eve, file = 'labeglo2/Transcriptomics/annotations/taxonomy_issues/Eve_taxons_counts.RData')
load(file = 'labeglo2/Transcriptomics/annotations/taxonomy_issues/Eve_taxons_counts.RData')

#taxon_freq_species_FiltContigs_Eve_muscle <- data.frame(Species = rep('Eve, muscle', 4),
#                                                        Taxon = factor(c(as.vector(taxon_freq_FiltContigs_king$x)), 
#                                                                       levels = c('Metazoa', 'Bacteria and Protozoa', 'Viridiplantae', 'Fungi')),
#                                                        Frequency = c(taxon_freq_FiltContigs_king$freq))
taxon_freq_species_FiltContigs_Eve_muscle <- data.frame(Species = rep('Eve, muscle', 7),
                                             Taxon = factor(c(as.vector(taxon_freq_FiltContigs_king$x)[-1], as.vector(unicel_count$x)), 
                                                            levels = c('Metazoa', 'Unicellular eukaryota', 'Bacteria', 
                                                                       'Viridiplantae', 'Fungi', 'Archaea', 'Viruses')),
                                             Frequency = c(taxon_freq_FiltContigs_king$freq[-1], unicel_count$freq))
taxon_freq_species_FiltContigs_Eve_muscle$Perc <- round(sapply(taxon_freq_species_FiltContigs_Eve_muscle$Frequency, 
                                                  function(x) x/sum(taxon_freq_species_FiltContigs_Eve_muscle$Frequency) * 100), 1)
#save(taxon_freq_species_FiltContigs_Eve_muscle, file = 'labeglo2/Transcriptomics/annotations/taxonomy_issues/EveMuscle_taxons_counts.RData')
load(file = 'labeglo2/Transcriptomics/annotations/taxonomy_issues/EveMuscle_taxons_counts.RData')

taxon_freq_species_FiltContigs_Gla <- rbind(taxon_freq_species_FiltContigs_Gla, 
                                            data.frame(Species = 'Gla',
                                                       Taxon = 'Bacteria and Archaea',
                                              Frequency = taxon_freq_species_FiltContigs_Gla$Frequency[4] + 
                                              taxon_freq_species_FiltContigs_Gla$Frequency[5],
                                            Perc = taxon_freq_species_FiltContigs_Gla$Perc[4] + 
                                              taxon_freq_species_FiltContigs_Gla$Perc[5]))
taxon_freq_species_FiltContigs_Gla <- taxon_freq_species_FiltContigs_Gla[-c(4,5),]
taxon_freq_species_FiltContigs_Gla$Taxon <- factor(taxon_freq_species_FiltContigs_Gla$Taxon,
                                                      levels = c('Metazoa', 'Unicellular eukaryota', 'Bacteria and Archaea', 
                                                                 'Viridiplantae', 'Fungi', 'Viruses'))
taxon_freq_species_FiltContigs_Ecy <- rbind(taxon_freq_species_FiltContigs_Ecy, 
                                            data.frame(Species = 'Ecy',
                                                       Taxon = 'Bacteria and Archaea',
                                                       Frequency = taxon_freq_species_FiltContigs_Ecy$Frequency[4] + 
                                                         taxon_freq_species_FiltContigs_Ecy$Frequency[5],
                                                       Perc = taxon_freq_species_FiltContigs_Ecy$Perc[4] + 
                                                         taxon_freq_species_FiltContigs_Ecy$Perc[5]))
taxon_freq_species_FiltContigs_Ecy <- taxon_freq_species_FiltContigs_Ecy[-c(4,5),]
taxon_freq_species_FiltContigs_Ecy$Taxon <- factor(taxon_freq_species_FiltContigs_Ecy$Taxon,
                                                   levels = c('Metazoa', 'Unicellular eukaryota', 'Bacteria and Archaea', 
                                                              'Viridiplantae', 'Fungi', 'Viruses'))
taxon_freq_species_FiltContigs_Eve <- rbind(taxon_freq_species_FiltContigs_Eve, 
                                            data.frame(Species = 'Eve',
                                                       Taxon = 'Bacteria and Archaea',
                                                       Frequency = taxon_freq_species_FiltContigs_Eve$Frequency[4] + 
                                                         taxon_freq_species_FiltContigs_Eve$Frequency[5],
                                                       Perc = taxon_freq_species_FiltContigs_Eve$Perc[4] + 
                                                         taxon_freq_species_FiltContigs_Eve$Perc[5]))
taxon_freq_species_FiltContigs_Eve <- taxon_freq_species_FiltContigs_Eve[-c(4,5),]
taxon_freq_species_FiltContigs_Eve$Taxon <- factor(taxon_freq_species_FiltContigs_Eve$Taxon,
                                                   levels = c('Metazoa', 'Unicellular eukaryota', 'Bacteria and Archaea', 
                                                              'Viridiplantae', 'Fungi', 'Viruses'))
taxon_freq_species_FiltContigs_Eve_muscle <- rbind(taxon_freq_species_FiltContigs_Eve_muscle, 
                                            data.frame(Species = 'Eve, muscle',
                                                       Taxon = 'Bacteria and Archaea',
                                                       Frequency = taxon_freq_species_FiltContigs_Eve_muscle$Frequency[4] + 
                                                         taxon_freq_species_FiltContigs_Eve_muscle$Frequency[5],
                                                       Perc = taxon_freq_species_FiltContigs_Eve_muscle$Perc[4] + 
                                                         taxon_freq_species_FiltContigs_Eve_muscle$Perc[5]))
taxon_freq_species_FiltContigs_Eve_muscle <- taxon_freq_species_FiltContigs_Eve_muscle[-c(4,5),]
taxon_freq_species_FiltContigs_Eve_muscle$Taxon <- factor(taxon_freq_species_FiltContigs_Eve_muscle$Taxon,
                                                   levels = c('Metazoa', 'Unicellular eukaryota', 'Bacteria and Archaea', 
                                                              'Viridiplantae', 'Fungi', 'Viruses'))

all_species_freq <- rbind(taxon_freq_species_FiltContigs_Gla, taxon_freq_species_FiltContigs_Ecy)
all_species_freq <- rbind(all_species_freq, taxon_freq_species_FiltContigs_Eve)  
all_species_freq <- rbind(all_species_freq, taxon_freq_species_FiltContigs_Eve_muscle)  

all_species_freq$Species <- factor(all_species_freq$Species, 
                                   levels = c('Eve, muscle', 'Eve', 'Ecy', 'Gla'))
all_species_freq <- subset(all_species_freq, Taxon != 'Viruses')
cbPalette <- c("#E69F00", "#56B4E9", "#F0E442", "#009E73", "#CC79A7", "#0072B2")
ggplot(all_species_freq, aes(Species, Frequency)) +
  geom_col(aes(fill = Taxon)) +
  geom_text(aes(group = Taxon,
                label = comma(ifelse(Frequency > 20000, Frequency, ''))),
            position = position_stack(reverse = F, vjust = .5), colour = "black", size = 3, na.rm = T) +
  coord_flip() +
#  scale_x_discrete(position = position_stack(reverse = TRUE)) +
  xlab("") +
  ylab("Number of contigs") +
  scale_fill_manual("Taxon", values = cbPalette) +
  theme_light() +
  theme(axis.text.y = element_text(hjust=0.5))

ggsave('labeglo2/Transcriptomics/annotations/all_species_freq_enlarged.png', scale=1)

ggplot(all_species_freq, aes(Species, Perc)) +
  geom_col(aes(fill = Taxon)) +
  geom_text(aes(group = Taxon, label = ifelse(Perc > 12, Perc, '')),
            position = position_stack(reverse = F, vjust = .5), colour = "black", size = 3) +
  coord_flip() +
  xlab("") +
  ylab("Contigs (%)") +
  scale_fill_manual("Taxon", values = cbPalette) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_text(hjust=0.5))
ggsave('labeglo2/Transcriptomics/annotations/all_species_perc_enlarged.png', scale=1)

