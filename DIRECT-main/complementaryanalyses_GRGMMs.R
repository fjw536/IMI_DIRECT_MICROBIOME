## Commplementary analyses

# setwd
setwd('/home/scratch/margar')

# libraries
library(ggplot2); library(ggpmisc); library(ggbeeswarm); library(ggpubr); library(tidyverse)

## Is gene richness associated to reduced MGS? ----
mtdt.bsl <- read.table('data/curatedmetadata/WP2.1_BL_MetadataWithInfl.txt', header=T, sep='\t')
qmp.bsl <- read.table('data/WP2.1_775_IGR_MGS-QMP_igrHA_Ent_M.csv', header=T, sep='\t', row.names=1)
qmp.bsl <- qmp.bsl[,7:727]

qmp.bsl <- data.frame(t(qmp.bsl))
mgs.rich <- colSums(qmp.bsl!=0)
names(mgs.rich) == mtdt.bsl$samplerenamed
  
mtdt.bsl$MGSRich <- mgs.rich
mtdt.bsl$DSGARLTertileGroup <- factor(mtdt.bsl$DSGARLTertileGroup, levels=c('LGC', 'HGC', 'EGC'))

p1 <- ggplot(mtdt.bsl, aes(MGSRich, DSGeneAbundRichness)) + 
  geom_point(aes(color=DSGARLTertileGroup)) + stat_smooth(method = 'lm', color='black') +
  stat_cor(method='pearson') + 
  scale_color_manual(values=c('darkblue', 'gray90', 'yellow')) + 
  theme(panel.grid.major = element_line(colour = "gray85",  linetype = "dashed"), 
        axis.title = element_text(size = 15, face = "bold"), axis.text = element_text(size = 12, face = "bold"),
        panel.background = element_rect(fill = NA), legend.key = element_rect(fill = NA), 
        legend.background = element_rect(fill = NA), legend.position = c(0.9, 0.15), 
        panel.border=element_rect(fill=NA)) +
  labs(x = "MGS Richness", y = "Gene richness",  colour = "Gene richness")

ggExtra::ggMarginal(
  p = p1,
  type = 'density',
  margins = 'both',
  size = 5,
  colour = 'black',
  groupFill = TRUE
)

p2 <- ggplot(mtdt.bsl, aes(DSGARLTertileGroup, MGSRich), group = interaction(DSGARLTertileGroup)) + 
  geom_violin(alpha = .1, aes(fill=DSGARLTertileGroup, color=DSGARLTertileGroup)) +
  geom_quasirandom(alpha = .6, dodge.width = .9, aes(color = DSGARLTertileGroup)) +
  geom_boxplot(width = .1, alpha = .7, outlier.colour = NA, position = position_dodge(.9), fill='white')  +
  stat_compare_means() + scale_color_manual(values=c('darkblue', 'gray90', 'yellow')) +
  scale_fill_manual(values=c('darkblue', 'gray90', 'yellow')) +
  stat_compare_means(comparisons = list(c('LGC', 'HGC'),
                                        c('LGC', 'EGC'),
                                        c('HGC', 'EGC'))) +
  theme(panel.grid.major = element_line(colour = "gray85",  linetype = "dashed"), 
        axis.title = element_text(size = 15, face = "bold"), axis.text.y = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 15, face = "bold"),
        panel.background = element_rect(fill = NA), legend.key = element_rect(fill = NA), 
        legend.background = element_rect(fill = NA), legend.position = 'none', 
        panel.border=element_rect(fill=NA)) +
  labs(y = "MGS Richness", x = NULL,  colour = "Gene richness", fill='Gene richness')

pdf('/home/margar/download/export_supplFIGRICH.pdf', height = 6, width = 15)
cowplot::plot_grid(ggExtra::ggMarginal(
  p = p1,
  type = 'density',
  margins = 'both',
  size = 5,
  colour = 'black',
  groupFill = TRUE
), p2)
dev.off()

mtdt.bsl %>% 
  select(c('DSGARLTertileGroup', 'MGSRich')) %>% 
  group_by(DSGARLTertileGroup) %>% 
  summarise(across(where(is.numeric), list(mean=mean, sd=sd, median=median)))

## and what about enterotypes? ----
p2 <- ggplot(mtdt.bsl, aes(enterotype, MGSRich), group = interaction(enterotype)) + 
  geom_violin(alpha = .1, aes(fill=enterotype, color=enterotype)) +
  geom_quasirandom(alpha = .6, dodge.width = .9, aes(color = enterotype)) +
  geom_boxplot(width = .1, alpha = .7, outlier.colour = NA, position = position_dodge(.9), fill='white')  +
  stat_compare_means() + scale_color_manual(values=c('darkolivegreen4', 'red', 'cadetblue3', 'darkgoldenrod2')) +
  scale_fill_manual(values=c('darkolivegreen4', 'red', 'cadetblue3', 'darkgoldenrod2')) +
  stat_compare_means(comparisons = list(c('Bact1', 'Bact2'),
                                        c('Bact1', 'Prev'),
                                        c('Bact1', 'Rum'),
                                        c('Bact2', 'Prev'),
                                        c('Bact2', 'Rum'),
                                        c('Prev', 'Rum'))) +
  theme(panel.grid.major = element_line(colour = "gray85",  linetype = "dashed"), 
        axis.title = element_text(size = 15, face = "bold"), axis.text.y = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 15, face = "bold"),
        panel.background = element_rect(fill = NA), legend.key = element_rect(fill = NA), 
        legend.background = element_rect(fill = NA), legend.position = 'right', 
        panel.border=element_rect(fill=NA)) +
  labs(y = "MGS Richness", x = NULL,  colour = "Enterotype", fill='Enterotype')

mtdt.bsl %>% 
  select(c('enterotype', 'MGSRich')) %>% 
  group_by(enterotype) %>% 
  summarise(across(where(is.numeric), list(mean=mean, sd=sd, median=median)))

## Differential GMMs ----
gmm <- read.table('data/WP2.1_775_IGR_GMM-QMP_igrHA_Ent_M.csv', header=T, row.names=1)
gmm.annotation <- read.table('data/GMMs.v1.07.names', header=F, sep='\t')

selected.gmms <- c('MF0030', 'MF0031', 'MF0032',
                   'MF0088', 'MF0089', 'MF0093',
                   'MF0094', 'MF0095', 'MF0099', 
                   'MF0100', 'MF0101', 'MF0081')

gmm <- gmm[,names(gmm)%in%selected.gmms]
row.names(gmm) == mtdt.bsl$samplerenamed
gmm <- data.frame(gmm, 'GR' = mtdt.bsl$DSGARLTertileGroup, 'Entero'=mtdt.bsl$enterotype)
gmm.annotation <- gmm.annotation[gmm.annotation$V1%in%names(gmm),]
gmm.annotation <- gmm.annotation[match(names(gmm)[1:10], gmm.annotation$V1),]
names(gmm)[1:10] <- gmm.annotation$V2

gmm$GR <- factor(gmm$GR, levels=c('LGC', 'HGC', 'EGC'))
gmm.m <- reshape2::melt(gmm)
ggplot(gmm.m, aes(GR, log10(value), group = interaction(GR))) + 
         geom_violin(alpha = .1, aes(fill=GR, color=GR)) +
         geom_quasirandom(alpha = .6, dodge.width = .9, aes(color = GR)) +
         geom_boxplot(width = .1, alpha = .7, outlier.colour = NA, position = position_dodge(.9), fill='white')  +
  facet_wrap(~variable, scales='free') +
         stat_compare_means() + scale_color_manual(values=c('darkblue', 'gray90', 'yellow')) +
         scale_fill_manual(values=c('darkblue', 'gray90', 'yellow')) +
         stat_compare_means(comparisons = list(c('LGC', 'HGC'),
                                               c('LGC', 'EGC'),
                                               c('HGC', 'EGC'))) +
         theme(panel.grid.major = element_line(colour = "gray85",  linetype = "dashed"), 
               axis.title = element_text(size = 15, face = "bold"), axis.text.y = element_text(size = 12, face = "bold"),
               axis.text.x = element_text(size = 15, face = "bold"),
               panel.background = element_rect(fill = NA), legend.key = element_rect(fill = NA), 
               legend.background = element_rect(fill = NA), legend.position = 'right', 
               panel.border=element_rect(fill=NA)) +
         labs(y = "log10(GMM counts)", x = NULL,  colour = "Gene richness", fill='Gene richness')
       
ggplot(gmm.m, aes(Entero, log10(value), group = interaction(Entero))) + 
  geom_violin(alpha = .1, aes(fill=Entero, color=Entero)) +
  geom_quasirandom(alpha = .6, dodge.width = .9, aes(color = Entero)) +
  geom_boxplot(width = .1, alpha = .7, outlier.colour = NA, position = position_dodge(.9), fill='white')  +
  facet_wrap(~variable, scales='free') +
  stat_compare_means() + scale_color_manual(values=c('darkolivegreen4', 'red', 'cadetblue3', 'darkgoldenrod2')) +
  scale_fill_manual(values=c('darkolivegreen4', 'red', 'cadetblue3', 'darkgoldenrod2')) +
  stat_compare_means(comparisons = list(c('Bact1', 'Bact2'),
                                        c('Bact1', 'Prev'),
                                        c('Bact1', 'Rum'),
                                        c('Bact2', 'Prev'),
                                        c('Bact2', 'Rum'),
                                        c('Prev', 'Rum'))) +
  theme(panel.grid.major = element_line(colour = "gray85",  linetype = "dashed"), 
        axis.title = element_text(size = 15, face = "bold"), axis.text.y = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 15, face = "bold"),
        panel.background = element_rect(fill = NA), legend.key = element_rect(fill = NA), 
        legend.background = element_rect(fill = NA), legend.position = 'right', 
        panel.border=element_rect(fill=NA)) +
  labs(y = "log10(GMM counts)", x = NULL,  colour = "Enterotype", fill='Enterotype')

mtdt.fw <- read.table('data/curatedmetadata/WP2.1_BlFw-641_MetadataImputed_OnlyBIOCHEM_2.txt', header=T, sep = '\t')
mtdt.fw <- mtdt.fw[mtdt.fw$timepoint=='Follow-up',]
gmm <- gmm[row.names(gmm)%in%mtdt.fw$sample.renamed,]
gmm$progression <- factor(mtdt.fw$Gly.Cat, levels=c('NGR', 'IGR', 'T2D'))
gmm.m <- reshape2::melt(gmm)
ggplot(gmm.m, aes(progression, log10(value), group = interaction(progression))) + 
  geom_violin(alpha = .1, aes(fill=progression, color=progression)) +
  geom_quasirandom(alpha = .6, dodge.width = .9, aes(color = progression)) +
  geom_boxplot(width = .1, alpha = .7, outlier.colour = NA, position = position_dodge(.9), fill='white')  +
  facet_wrap(~variable, scales='free') +
  stat_compare_means() + scale_color_manual(values=c('darkolivegreen3', 'red', 'darkorchid')) +
  scale_fill_manual(values=c('darkolivegreen3', 'red', 'darkorchid')) +
  stat_compare_means(comparisons = list(c('NGR', 'IGR'),
                                        c('NGR', 'T2D'),
                                        c('IGR', 'T2D'))) +
  theme(panel.grid.major = element_line(colour = "gray85",  linetype = "dashed"), 
        axis.title = element_text(size = 15, face = "bold"), axis.text.y = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 15, face = "bold"),
        panel.background = element_rect(fill = NA), legend.key = element_rect(fill = NA), 
        legend.background = element_rect(fill = NA), legend.position = 'right', 
        panel.border=element_rect(fill=NA)) +
  labs(y = "log10(GMM counts)", x = NULL,  colour = "Gene richness", fill='Gene richness')

## Tpred score - alpha diversity association to lifestyle
reb <- read.table('/home/Teams/teamVIP/Joram/Results/JMP_Results_NIHR_Predict_OriginalFileNamed DIRECT_WP2_join_TargetedMetabolites_Filtered_04022018.csv', header = T, sep='\t')
head(reb)
reb <- reb[reb$DIRECT.ID%in%mtdt.bsl$StudyID,]
reb <- reb[match(mtdt.bsl$StudyID, reb$DIRECT.ID),]

mtdt.bsl$Tpred <- reb$Tpred..mccv.
mtdt.bsl$Smoking.Status.BL <- plyr::revalue(as.factor(mtdt.bsl$Smoking.Status.BL), c('ex-smoker'='non-smoker', 'never'='non-smoker', 'current-smoker'='smoker'))

reg.model.alphadiv <- lm(DSMGSAbundShannon~HDI*Value.vm.hpf.mean.BL*Smoking.Status.BL*Tpred + Age + Gender + CenterID, data = mtdt.bsl)
coef(reg.model.alphadiv)
summary(reg.model.alphadiv)


reg.model.alphadiv2 <- lm(DSMGSAbundShannon~HDI+Value.vm.hpf.mean.BL+Smoking.Status.BL+Tpred + Age + Gender + CenterID, data = mtdt.bsl)
summary(reg.model.alphadiv2)

reg.model.alphadiv3 <- lm(DSMGSAbundShannon~HDI*Value.vm.hpf.mean.BL*Smoking.Status.BL*Tpred, data = mtdt.bsl)
coef(reg.model.alphadiv3)
summary(reg.model.alphadiv3)

# stats differences depending on tertiles
tpred.aov <- aov(Tpred~DSGARLTertileGroup, data=mtdt.bsl)
summary(tpred.aov)
TukeyHSD(tpred.aov)
ggplot(mtdt.bsl, aes(DSGARLTertileGroup, Tpred)) + geom_boxplot()

HDI.aov <- aov(HDI~DSGARLTertileGroup, data=mtdt.bsl)
summary(HDI.aov)
TukeyHSD(HDI.aov)
ggplot(mtdt.bsl, aes(DSGARLTertileGroup, HDI)) + geom_boxplot()

physical.aov <- aov(Value.vm.hpf.mean.BL~DSGARLTertileGroup, data=mtdt.bsl)
summary(physical.aov)
TukeyHSD(physical.aov)
ggplot(mtdt.bsl, aes(DSGARLTertileGroup, Value.vm.hpf.mean.BL)) + geom_boxplot()

smoker.table <- mtdt.bsl %>% 
  select(c('Smoking.Status.BL', 'DSGARLTertileGroup')) %>% 
  group_by(DSGARLTertileGroup, Smoking.Status.BL) %>% 
  tally() %>% 
  pivot_wider(values_from=n, names_from=Smoking.Status.BL)

row.names(smoker.table) <- smoker.table$DSGARLTertileGroup
chisq.test(smoker.table[,2:3])
smoker.table %>% 
  pivot_longer(cols=c('smoker', 'non-smoker')) %>% 
  mutate(name=factor(name, levels=c('smoker', 'non-smoker'))) %>% 
  ggplot(aes(DSGARLTertileGroup, value, fill=name)) + geom_col(position = 'dodge')

##