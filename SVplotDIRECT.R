setwd('~/Desktop')

library(ggplot2); library(tidyverse)

abd.mets <- read.table('SVsmetabolites_var.txt', header=T, sep='\t')
gc()

pdf('SVs_variable_metabolites.pdf', height = 25, width = 12)
abd.mets %>%
  filter(fdr<.1) %>%
  mutate('delta' = as.numeric(delta)) %>%
  ggplot(aes(met.cluster, bacteria, fill = delta)) +
  geom_tile() +
  scale_fill_gradient2(low='blue', high='red') +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust = .5), panel.border = element_rect(fill=NA),
        panel.grid.major = element_line(linetype='dashed', color='gray75')) +
  labs(x = 'Metabolite Clusters', y = 'Bacterial regions (abbreviated)')
dev.off()
