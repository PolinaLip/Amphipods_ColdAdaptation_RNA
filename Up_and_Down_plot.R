library(ggplot2)

overrepr <- data.frame(Species = rep(c("G. lacustris","E. cyaneus","E. verrucosus"), 
                                     each = 5),
                       Fold_change = rep(c("Up, >4x", "Up, >2x", "Not changed",
                                           "Down, >2x", "Down, >4x"), 3),
                       Percentage = c(1.19, 9.55, 100-9.55-1.19-17.93-12.5, 12.5, 17.93,
                                      0.2, 0.52, 100-0.52-0.2-0.77-0.089, 0.77, 0.089,
                                      0.66, 0.94, 100-0.66-0.94-0.83-0.48, 0.83, 0.48))
overrepr$Fold_change <- factor(overrepr$Fold_change, 
                              levels = c("Up, >4x", "Up, >2x", "Not changed",
                                          "Down, >2x", "Down, >4x"))


ggplot(overrepr, aes(Species, Percentage)) +
  geom_col(aes(fill = Fold_change, alpha = Fold_change),
           position = position_stack(reverse = TRUE)) +
  coord_flip() +
  xlab("") +
  scale_y_continuous("Genes, %", breaks=seq(0, 100, 20)) +
  scale_fill_manual("Fold change:",
    values = c("#ff3333", "#ffb3b3", "#808080", "#b3d9ff", "#0081cc")) +
#  values = c("#E64241", "#FF9999", "#354E6C", "#99CCFF", "#EAAC31")) +
  scale_alpha_manual("Fold change:", values=c(1, 1, 0.75, 1, 1)) +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
#        panel.grid.minor = element_blank(),
#        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(face = "italic", size = 14))

ggsave('labeglo2/Transcriptomics/UpandDown_genes_all_species5.png', 
        scale = 0.8, height = 3, width = 9)

### To rebuild the plot for the better reading:
overrepr <- subset(overrepr, Fold_change != 'Not changed')
overrepr$perc_signed <- ifelse(grepl('Down', overrepr$Fold_change),
                               -overrepr$Percentage, overrepr$Percentage)

overrepr$large <- ifelse(grepl('>4x', overrepr$Fold_change), 
                         '>4x', '>2x')

ggplot(overrepr, aes(Species, perc_signed)) +
  geom_bar(aes(fill = Species, alpha = large, colour = Species), stat = "identity") +
  coord_flip() +
  scale_alpha_manual("Fold change:", values = c(.5, 1)) +
  scale_color_manual(values = c("#56B4E9", "#009E73", "#E69F00")) +
  scale_fill_manual(values = c("#56B4E9", "#009E73", "#E69F00")) +
  scale_y_continuous('Genes, %', labels = function(x) abs(x)) +
  geom_hline(yintercept = 0) +
  theme_light() +
  guides(fill = F, colour = F) +
  theme(axis.text.y = element_text(face = "italic", size = 14),
        axis.title.y = element_blank())

overrepr$combined_labels <- factor(sprintf('%s, %s', overrepr$Species, 
                                           overrepr$large),
                                   levels = c("G. lacustris, >2x",
                                              "E. cyaneus, >2x",
                                              "E. verrucosus, >2x",
                                              "",
                                              "G. lacustris, >4x",
                                              "E. cyaneus, >4x",
                                              "E. verrucosus, >4x"))  

ggplot(overrepr, aes(Species, perc_signed)) +
  geom_bar(aes(fill = combined_labels, colour = combined_labels, 
               alpha = combined_labels), 
           stat = "identity", size = 0.3) +
  coord_flip() +
  scale_color_manual('Fold change:', values = c("#E69F00", "#56B4E9", "#009E73",
                                                "white",
                                                "#E69F00", "#56B4E9", "#009E73"),
                     labels=c('', '>2x', '', '\n', '', '>4x', ''), drop = F) +
  scale_fill_manual('Fold change:', values = c("#E69F00", "#56B4E9", "#009E73",
                                               "white",
                                               "#E69F00", "#56B4E9", "#009E73"),
                     labels=c('', '>2x', '', '\n', '', '>4x', ''), drop = F) +
  scale_alpha_manual('Fold change:', values = c(.3,.3,.3,1,1,1,1),
                     labels=c('', '>2x', '', '\n', '', '>4x', ''), drop = F) +
  scale_y_continuous('Genes, %', labels = function(x) abs(x)) +
  geom_hline(yintercept = 0, size = 0.2) +
  theme_light() +
  theme(axis.text.y = element_text(face = "italic", size = 14),
        axis.title.y = element_blank())

ggsave('labeglo2/Transcriptomics/UpandDown_genes_all_species_newversion.pdf',
       scale = 0.8, height = 3, width = 7)

                                  