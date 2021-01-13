library(stringr)
library(plyr)
library(ggplot2)
library(remotes)
library(ggsci)
library(scales)
library(formattable)

# The annotation file using here is came after annotation by Diamond.   
# Beforehand, it is essential to get rid of duplicates from the annotation file.
# The recover of the taxonomy for N/A fields are done with fix_NA_in_Annotation.py

### 1. Upload the data
annotation_data <- function(path_to_csv) # read the annotation file
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
species_filt_annot <- annotation_data('labeglo2/Transcriptomics/labeglo1_HS/Ecy_filtered_contigs_HS_rnaspades_woNA.diamond.csv') # Ecy, HS, RNAspades
species_filt_annot <- annotation_data('labeglo2/Transcriptomics/labeglo1_HS/Eve_filtered_contigs_HS_rnaspades_woNA.diamond.csv') # Eve, HS, RNAspades
species_filt_annot <- annotation_data('labeglo2/Transcriptomics/labeglo1_HS/Gla_filtered_contigs_HS_rnaspades_woNA.diamond.csv') # Gla, HS, RNAspades
species_filt_annot <- annotation_data('/home/polina/labeglo2/Transcriptomics/annotations/ltce_RNAspades/Ecy_filtered_contigs_ltce_rnaspades_woNA.diamond.csv')
species_filt_annot <- annotation_data('/home/polina/labeglo2/Transcriptomics/annotations/ltce_RNAspades/Gla_filtered_contigs_ltce_rnaspades_woNA.diamond.csv')
species_filt_annot <- annotation_data('/home/polina/labeglo2/Transcriptomics/annotations/ltce_RNAspades/Eve_filtered_contigs_ltce_rnaspades_woNA.diamond.csv')

### 2. Get rid of the strange annotations:
count(species_filt_annot$sskingdoms)
species_filt_annot <- subset(species_filt_annot, sskingdoms != "0;Eukaryota")
species_filt_annot <- subset(species_filt_annot, sskingdoms != "Bacteria;Archaea")
species_filt_annot <- subset(species_filt_annot, sskingdoms != "Bacteria;Eukaryota")
species_filt_annot <- subset(species_filt_annot, sskingdoms != "0")
species_filt_annot <- subset(species_filt_annot, sskingdoms != "0;Bacteria")

count(species_filt_annot$sskingdoms)

count(species_filt_annot$skingdoms)
species_filt_annot <- subset(species_filt_annot, skingdoms != "0;Metazoa")
species_filt_annot <- subset(species_filt_annot, skingdoms != "Fungi;Metazoa")
species_filt_annot <- subset(species_filt_annot, skingdoms != "0;Viridiplantae;Metazoa")
count(species_filt_annot$skingdoms)

### 3. Prepare info for the statistic plot about the ratio of different taxons 
taxon_freq_FiltContigs_king <- count(species_filt_annot$skingdoms)

taxon_freq_FiltContigs_king$x <- as.vector(taxon_freq_FiltContigs_king$x)
#taxon_freq_FiltContigs_king$x[1] <- "Bacteria and Protozoa"
unicel <- subset(species_filt_annot, skingdoms == '0')
unicel_count <- count(unicel$sskingdoms)
unicel_count$x <- sub('Eukaryota', 'Unicellular eukaryota', unicel_count$x)

unicel_count$freq[4] <- unicel_count$freq[4] + 116 + 22 + 6 # for Ecy HS RNAspades, we have Bamfordvirae, Orthornavirae, and Pararnavirae, which are viruses, but did not appear in the Viruses
unicel_count$freq[4] <- unicel_count$freq[4] + 78 + 13 + 4 # for Eve HS RNAspades, we have Bamfordvirae, Orthornavirae, and Pararnavirae, which are viruses, but did not appear in the Viruses
unicel_count$freq[4] <- unicel_count$freq[4] + 55 + 12 + 8  # for Eve HS RNAspades, we have Bamfordvirae, Orthornavirae, and Pararnavirae, which are viruses, but did not appear in the Viruses
unicel_count$freq[4] <- unicel_count$freq[4] + 127 + 19 + 10 + 9 + 2 # for Ecy ltce RNAspades, we have Bamfordvirae, Orthornavirae, Pararnavirae, Sangervirae, and Shotokuvirae, which are viruses
unicel_count$freq[4] <- unicel_count$freq[4] + 177 + 18 + 6 # for Gla ltce RNAspades, we have Bamfordvirae, Orthornavirae, Pararnavirae, which are viruses
unicel_count$freq[4] <- unicel_count$freq[4] + 244 + 26 + 7 + 4 # for Eve ltce RNAspades, we have Bamfordvirae, Orthornavirae, Pararnavirae, and Shotokuvirae which are viruses

taxon_freq_FiltContigs_king <- subset(taxon_freq_FiltContigs_king, 
                                      x != 'Bamfordvirae' & x != 'Orthornavirae' & x != 'Pararnavirae' & x != 'Shotokuvirae' & x != 'Sangervirae')

taxon_freq_species_FiltContigs <- 
  data.frame(Species = rep('Gla', 7),
  Taxon = factor(c(as.vector(taxon_freq_FiltContigs_king$x)[-1], 
                   as.vector(unicel_count$x)), 
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

### 4. Use this plot if you don't split unicellular organisms from bacteria and archeae
ggplot(taxon_freq_species_FiltContigs, aes(x = Species, y = Frequency)) +
  geom_col(aes(fill = Taxon)) +
  #  geom_text(aes(group = Frequency, label = Taxon), 
  #            position = position_stack(reverse = TRUE, vjust = .5)) +
  coord_flip() +
  ylab("Number of contigs") +
  scale_fill_manual(values = c("#FF6F00FF", "#C71000FF", "#008EA0FF", "#8A4198FF")) +
  theme_bw()

### 5. To plot the distribution of filtered contigs from Metazoa 
species_filtContigs_metazoa <- subset(species_filt_annot, skingdoms == 'Metazoa')
count(species_filtContigs_metazoa$sphylums)
#species_filtContigs_metazoa[which(species_filtContigs_metazoa$sphylums != 'Arthropoda;Chordata;Hemichordata'),]
species_filtContigs_metazoa <- subset(species_filtContigs_metazoa, sphylums != 'Arthropoda;Chordata;Hemichordata')
species_filtContigs_metazoa <- subset(species_filtContigs_metazoa, sphylums != 'Platyhelminthes;Mollusca')
#species_filtContigs_metazoa <- subset(species_filtContigs_metazoa, sphylums != 'Brachiopoda;Hemichordata')
species_filtContigs_metazoa <- subset(species_filtContigs_metazoa, sphylums != 'Arthropoda;Chordata')

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

ggsave('labeglo2/Transcriptomics/annotations/ltce_RNAspades/metazoa_Eve.png',
       scale = 1.2)

### 6. Save csv with contigs assigned to Metazoa:
write.table(species_filtContigs_metazoa, 
  file = 'labeglo2/Transcriptomics/annotations/ltce_RNAspades/metazoa_contigs_Eve_ltce_rnaspades.csv', 
  sep = '\t', quote = F, col.names = F)

### 7. Plot the distribution of contigs assigned to Arthropoda
species_FiltContig_arthropoda <- subset(species_filtContigs_metazoa, 
                                        sphylums == 'Arthropoda')
species_FiltContig_arthropoda_freq <- count(species_FiltContig_arthropoda$sscinames) # 31,202 records for Ecy; 36,542 for Eve; 37,365 for Gla; 21,949 for Eve Muscle
more_10_arthropoda_FiltContig_species <- subset(species_FiltContig_arthropoda_freq, 
                                                freq > 20) # 30,493 records for Ecy (>10); 34,886 for Eve (> 20); 36,100 for Gla (>20);20,993 for Eve Muscle

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

ggsave('labeglo2/Transcriptomics/annotations/ltce_RNAspades/arthropoda_eve.png',
       scale = 1.3)
species_FiltContig_mollusca <- subset(species_filtContigs_metazoa, sphylums == 'Mollusca')
count(species_FiltContig_mollusca$sscinames)

species_FiltContig_chordata <- subset(species_filtContigs_metazoa, sphylums == 'Chordata')
View(count(species_FiltContig_chordata$sscinames))

### 8. Let's take only Arthropoda, Chordata and Mollusca for the transcripts quantification
species_FiltContig_Arthr_Chotdata_Mollusca <- subset(species_filtContigs_metazoa, 
                                                 sphylums == 'Arthropoda' |
                                                   sphylums == 'Chordata' |
                                                   sphylums == 'Mollusca')
count(species_FiltContig_Arthr_Chotdata_Mollusca$sphylums)
write.table(species_FiltContig_Arthr_Chotdata_Mollusca, file = 'labeglo2/Transcriptomics/annotations/arthr_chord_mollus_contigs_Gla.csv', sep = '\t', quote = F, col.names = F)

### 9. Plot distribution for all three species

#taxon_freq_species_FiltContigs_Gla <- data.frame(Species = rep('Gla', 4),
#                                             Taxon = factor(c(as.vector(taxon_freq_FiltContigs_king$x)), 
#                                                            levels = c('Metazoa', 'Bacteria and Protozoa', 'Viridiplantae', 'Fungi')),
#                                             Frequency = c(taxon_freq_FiltContigs_king$freq))
taxon_freq_species_FiltContigs_Gla <- data.frame(Species = rep('Gla', 7),
                                             Taxon = factor(c(as.vector(taxon_freq_FiltContigs_king$x)[-1], as.vector(unicel_count$x)), 
                                                            levels = c('Metazoa', 'Unicellular eukaryota', 'Bacteria', 
                                                                       'Viridiplantae', 'Fungi', 'Archaea', 'Viruses')),
                                             Frequency = c(taxon_freq_FiltContigs_king$freq[-1], unicel_count$freq))
taxon_freq_species_FiltContigs_Gla$Perc <- 
  round(sapply(taxon_freq_species_FiltContigs_Gla$Frequency, 
  function(x) x/sum(taxon_freq_species_FiltContigs_Gla$Frequency) * 100), 1)
#save(taxon_freq_species_FiltContigs_Gla, 
#     file = 'labeglo2/Transcriptomics/annotations/ltce_RNAspades/Gla_taxons_counts.RData')
load(file = 'labeglo2/Transcriptomics/annotations/taxonomy_issues/Eve_taxons_counts.RData')

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
#save(taxon_freq_species_FiltContigs_Ecy, file = 'labeglo2/Transcriptomics/annotations/ltce_RNAspades/Ecy_taxons_counts.RData')
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
#save(taxon_freq_species_FiltContigs_Eve, file = 'labeglo2/Transcriptomics/annotations/ltce_RNAspades/Eve_taxons_counts.RData')
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

### figure S2B (with number of contigs)
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

ggsave('labeglo2/Transcriptomics/annotations/ltce_RNAspades/all_species_freq.png', scale=1)

### figure S2B (with % of contigs)
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
ggsave('labeglo2/Transcriptomics/annotations/ltce_RNAspades/all_species_perc.png', scale=1)

