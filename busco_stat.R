library(ggplot2)

ecy_busco <- c(32.3, 51.0, 7.6, 9.1) # complete and single-copy, complete and duplicated, fragmented, missing BUSCOs
ecy_filter <- c(61.0, 11.0, 11.2, 16.8)
ecy_taxfilter <- c(61.1, 9.6, 11.3, 18.0)
gl_busco <- c(30.1, 53.9, 6.4, 9.6)
gl_filter <- c(52.5, 30.4, 6.1, 11.0)
gl_taxfilter <- c(52.5, 28.5, 7.5, 11.5)
eve_busco <- c(33.8, 49.9, 6.8, 9.5)
eve_filter <- c(63.8, 18.8, 7.0, 10.4)
eve_taxfilter <- c(63.9, 16.8, 7.7, 11.6)

buscos <- data.frame(
  Species = rep(c("Ecy", "Ecy", "Ecy", "Gla", "Gla", "Gla", "Eve", "Eve", "Eve"),
                each = 4),
  Filtered = rep(c("Before filtering", "Filtered", "TaxFiltered", 
                   "Before filtering", "Filtered", "TaxFiltered",
                   "Before filtering", "Filtered", "TaxFiltered"), each = 4),
  BUSCO_terms = rep(c("Complete and single-copy", "Complete and duplicated", "Fragmented", "Missing"), 9), 
  Values = c(32.3, 51.0, 7.6, 9.1, 61.0, 11.0, 11.2, 16.8, 61.1, 9.6, 11.3, 18.0, 
             30.1, 53.9, 6.4, 9.6, 52.5, 30.4, 6.1, 11.0, 52.5, 28.5, 7.5, 11.5,
             33.8, 49.9, 6.8, 9.5, 63.8, 18.8, 7.0, 10.4, 63.9, 16.8, 7.7, 11.6))

buscos$BUSCO_terms <- factor(buscos$BUSCO_terms,
                             levels=c("Complete and single-copy", "Complete and duplicated", "Fragmented", "Missing"))
buscos$Filtered <- factor(buscos$Filtered, 
                          levels=c("Before filtering", "Filtered", "TaxFiltered"))

ggplot(buscos, aes(Species, Values)) +
  geom_col(aes(fill = BUSCO_terms), position = position_stack(reverse = TRUE)) +
  geom_text(aes(group = BUSCO_terms, label = Values),
            position = position_stack(reverse = TRUE, vjust = .5), colour = "white", size = 3) +
  coord_flip() +
  xlab("") +
  ylab("% BUSCO") +
  scale_fill_manual("BUSCO terms (%)", values = c("#665191", "#a05195", "#ff7c43", "#ffa600")) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank()) +
  facet_wrap(~ Filtered)

ggsave('labeglo2/Transcriptomics/BUSCO_BeforeFilt_Filt_TaxFilt_3species.png', 
       scale = 1.3, height = 2.5, width = 8)

#####################################################
buscos <- data.frame(
  Assembler = rep(c("Trinity", "RNAspades"),
                each = 12),
  Filtered = rep(rep(c("Before filtering", "Filtered", "Wo contigs < 200bp"), 
                     each = 4), 2),
  BUSCO_terms = rep(c("Complete and single-copy", "Complete and duplicated", 
                      "Fragmented", "Missing"), 6), 
  Values = c(36.2, 44.1, 10.1, 9.6,
             64.7, 13.9, 10.3, 11.1,
             64.7, 13.9, 10.3, 11.1,
             62.7, 17.6, 9.1, 10.6,
             69.2, 11.1, 9.0, 10.7,
             69.2, 11.1, 9.0, 10.7))

buscos$BUSCO_terms <- factor(buscos$BUSCO_terms,
                             levels=c("Complete and single-copy", 
                                      "Complete and duplicated", 
                                      "Fragmented", "Missing"))
buscos$Filtered <- factor(buscos$Filtered, 
                          levels=c("Before filtering", "Filtered", 
                                   "Wo contigs < 200bp"))

ggplot(buscos, aes(Assembler, Values)) +
  geom_col(aes(fill = BUSCO_terms), position = position_stack(reverse = TRUE)) +
  geom_text(aes(group = BUSCO_terms, label = Values),
            position = position_stack(reverse = TRUE, vjust = .5), 
            colour = "white", size = 3) +
  coord_flip() +
  xlab("") +
  ylab("% BUSCO") +
  scale_fill_manual("BUSCO terms (%)", 
                    values = c("#665191", "#a05195", "#ff7c43", "#ffa600")) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank()) +
  facet_wrap(~ Filtered)

ggsave('labeglo2/Transcriptomics/BUSCO_BeforeFilt_Filt_Filt_ECY_2assemblers.png', 
       scale = 1.3, height = 1.5, width = 8.5)


