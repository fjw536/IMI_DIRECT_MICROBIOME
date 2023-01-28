###########################################################################################
################################## Biochem data analyses ################################## 
###########################################################################################

setwd('/home/scratch/margar') # directory
# libraries
library(ggplot2); library(reshape2); library(dplyr); library(tidyverse
                                                             )
## Baseline ----

## Follow-up ----

## Progression ----
all.data <- read.table('data/curatedmetadata/WP2.1_BlFw-641_MetadataImputed_OnlyBIOCHEM.txt', header = T)
all.data$sample.renamed[1:641] == all.data$sample.renamed[642:1282]
progression <- data.frame('baseline' = all.data$Gly.Cat[1:641], 'followup' = all.data$Gly.Cat[642:1282])
progression <- progression %>% 
  group_by(baseline, followup) %>% 
  tally()
progression$followup <- factor(progression$followup, levels=c('T2D', 'IGR', 'NGR'))
alluvial::alluvial(progression[,1:2], freq = progression$n, 
                   col = c('red', 'darkolivegreen4', 'darkorchid'), cex = 1.5, cex.axis = 1.5, blocks = T)

progression <- data.frame('baseline' = all.data$Gly.Cat.2[1:641], 'followup' = all.data$Gly.Cat.2[642:1282])
progression <- progression %>% 
  group_by(baseline, followup) %>% 
  tally()
progression$baseline <- factor(progression$baseline, levels=c('IGT', 'IFG'))
progression$followup <- factor(progression$followup, levels=c('T2D', 'IGT', 'IFG', 'NGR'))
alluvial::alluvial(progression[,1:2], freq = progression$n, 
                   col = c('red', 'darkred', 'darkolivegreen4', 'darkorchid'), cex = 1.5, cex.axis = 1.5, blocks = T)

# delta of log values

delta.df <- readRDS('DeltaofLogValues.rds')
delta.df$CenterID <- as.factor(delta.df$CenterID)
delta.df$Gender <- as.factor(delta.df$Gender)

delta.df.m <- melt(delta.df)
delta.df.m$Gly.Cat <- factor(delta.df.m$Gly.Cat, levels = c('NGR', 'IFG', 'IGT', 'T2D'))
ggplot(delta.df.m, aes(Gly.Cat, value), group = interaction(Gly.Cat)) + 
  geom_violin(alpha = .1, aes(fill=Gly.Cat, color=Gly.Cat)) +
  geom_quasirandom(alpha = .6, dodge.width = .9, aes(color = Gly.Cat)) +
  geom_boxplot(width = .1, alpha = .7, outlier.colour = NA, position = position_dodge(.9), fill='white') +
  facet_wrap(~variable, scales='free') +
  stat_compare_means(comparisons = list(c('NGR', 'IFG'),
                                        c('NGR', 'IGT'),
                                        c('NGR', 'T2D'), 
                                        c('IFG', 'T2D'),
                                        c('IGT', 'T2D'),
                                        c('IFG', 'IGT')), label = 'p.signif') +
  scale_color_manual(values = c('darkolivegreen4', 'red', 'darkred', 'purple')) +
  scale_fill_manual(values = c('darkolivegreen4', 'red', 'darkred', 'purple')) + 
  theme(panel.grid.major = element_line(colour = "gray85",linetype = "dashed"), 
        axis.text = element_text(size = 11, face = "bold", colour = "black"), 
        panel.background = element_rect(fill = NA), legend.title = element_text(face='bold'),
        legend.text = element_text(face='bold'), strip.background = element_rect(fill=NA),
        strip.text = element_text(face='bold', size =12), panel.border = element_rect(fill=NA), legend.position = 'none')

delta.df <- delta.df[complete.cases(delta.df),]
delta.df %>%
  select(-c('Age', 'CenterID')) %>% 
  group_by(Gly.Cat.2) %>% 
  dplyr::summarise(across(is.numeric, list(mean=mean, sd=sd, median=median))) %>% 
  pivot_longer(ends_with(c("_mean", "_sd", "_median")), names_to = c("name", ".value"), names_sep = "_") %>% 
  mutate(Gly.Cat = factor(Gly.Cat.2, levels=c('NGR', 'IFG', 'IGT', 'T2D'))) %>% 
  mutate(name = factor(name, levels=as.character(unique(name)))) %>% 
  ggplot(aes(reorder(name, desc(name)), mean, fill=Gly.Cat, color=Gly.Cat)) + geom_col(position = position_dodge2(reverse = TRUE), alpha = .6) +
  scale_color_manual(values = c('darkolivegreen4', 'red', 'darkred', 'purple')) +
  scale_fill_manual(values = c('darkolivegreen4', 'red', 'darkred', 'purple')) +
  coord_flip() + ylab(label = 'Mean delta value') + xlab(NULL) +
  theme(panel.grid.major = element_line(colour = "gray85",linetype = "dashed"), 
        axis.text = element_text(size = 11, face = "bold", colour = "black"), 
        axis.title = element_text(size=14, face='bold', color='black'),
        panel.background = element_rect(fill = NA), legend.title = element_text(face='bold'),
        legend.text = element_text(face='bold'), strip.background = element_rect(fill=NA),
        strip.text = element_text(face='bold', size =12), panel.border = element_rect(fill=NA))

delta.df$GlyCat3g <- delta.df$Gly.Cat
delta.df$GlyCat3g <- gsub('IFG', 'IGR', delta.df$GlyCat3g)
delta.df$GlyCat3g <- gsub('IGT', 'IGR', delta.df$GlyCat3g)

delta.df %>%
  select(-c('Age', "CenterID", 'Gly.Cat.2')) %>% 
  group_by(GlyCat3g) %>% 
  dplyr::summarise(across(is.numeric, list(mean=mean, sd=sd, median=median))) %>% 
  pivot_longer(ends_with(c("_mean", "_sd", "_median")), names_to = c("name", ".value"), names_sep = "_") %>% 
  mutate(GlyCat3g = factor(GlyCat3g, levels=c('NGR', 'IGR', 'T2D'))) %>% 
  mutate(name = factor(name, levels=as.character(unique(name)))) %>% 
  ggplot(aes(reorder(name, desc(name)), mean, fill=GlyCat3g, color=GlyCat3g)) + geom_col(position = position_dodge2(reverse = TRUE), alpha = .6) +
  scale_color_manual(values = c('darkolivegreen4', 'red', 'purple')) +
  scale_fill_manual(values = c('darkolivegreen4', 'red', 'purple')) +
  coord_flip() + ylab(label = 'Mean delta value') + xlab(NULL) +
  theme(panel.grid.major = element_line(colour = "gray85",linetype = "dashed"), 
        axis.text = element_text(size = 11, face = "bold", colour = "black"), 
        axis.title = element_text(size=14, face='bold', color='black'),
        panel.background = element_rect(fill = NA), legend.title = element_text(face='bold'),
        legend.text = element_text(face='bold'), strip.background = element_rect(fill=NA),
        strip.text = element_text(face='bold', size =12), panel.border = element_rect(fill=NA))

