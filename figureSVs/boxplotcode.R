## figure SVs
setwd('~/Desktop/01_Papers/04_DIRECT/figureSVs')

library(tidyverse); library(ggplot2); library(ggpubr); library(ggbeeswarm)

dels <- read.table('lmprDELS.txt', header=T, sep='\t')
dels <- dels[dels$fdr<.1,]
dels <- dels[complete.cases(dels),]

dels$direction <- ifelse(dels$delta>0, 'Positive', 'Negative')
dels$direction <- factor(dels$direction, levels=c('Positive', 'Negative'))
dels$pheno <-  ifelse(dels$delta>0, 'NGR', 'T2D')

pdf('progressionSV.pdf', width = 10, height = 4)
ggplot(dels, aes(delta, mgs, shape = direction, fill =pheno)) + 
  geom_errorbar(aes(xmin=delta-delta.low, xmax=delta+delta.up), width=.1) +
  geom_point(size=7) +
  geom_vline(xintercept = 0, linetype='dashed') +
  scale_shape_manual(values=c(24,25)) +
  scale_fill_manual(values=c('darkolivegreen3', 'darkorchid')) +
  guides(fill=F, shape=F) +
  theme(panel.background = element_blank(), panel.border = element_rect(fill=NA),
        axis.text.y = element_text(size=15, face='italic'),
        axis.text.x = element_text(size = 15),
        axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        panel.grid.major = element_line(linetype = 'dashed', color = 'gray90')) +
  labs(x= "Cliff's Delta Effect size", y = 'Bacterial regions')
dev.off()


eub.rect <- read.table('SVsplotsBXs.txt', header=T) 
eub.rect2 <- read.table('SVsplotsBXs2.txt', header=T) 

p2 <- eub.rect %>% 
  filter(!is.na(E.rectale.1706.1708)) %>% 
  mutate('E.rectale.1706.1708' = factor(E.rectale.1706.1708,levels = c(1,0))) %>% 
  mutate('E.rectale.1706.1708' = plyr::revalue(E.rectale.1706.1708, c('0'='With deletion', '1'='Without deletions'))) %>% 
  ggplot(aes(E.rectale.1706.1708, indolepropionate, color=E.rectale.1706.1708)) +
  geom_quasirandom() +
  geom_boxplot(alpha=.3) +
  geom_hline(yintercept = 0, linetype='dashed', color='gray80') +
  stat_compare_means() + theme(axis.line = element_line(size = 0.5, 
                                                        linetype = "solid"), axis.title = element_text(size = 15), 
                               panel.background = element_rect(fill = NA), legend.position = "none", legend.direction = "horizontal", 
                               plot.caption = element_text(face = "italic"), axis.text = element_text(size=15)) +
  labs(x = NULL, y = "Fasting plasma indolepropionate\n(normalized peak intensity)")

p1 <- eub.rect2 %>% 
  filter(!is.na(E.rectale.1701)) %>% 
  mutate('E.rectale.1701' = factor(E.rectale.1701,levels = c(1,0))) %>% 
  mutate('E.rectale.1701' = plyr::revalue(E.rectale.1701, c('0'='With deletion', '1'='Without deletions'))) %>% 
  ggplot(aes(E.rectale.1701, mean.gluc, color=E.rectale.1701)) +
  geom_quasirandom() +
  geom_boxplot(alpha=.3) +
  stat_compare_means() + 
  theme(axis.line = element_line(size = 0.5, linetype = "solid"), axis.title = element_text(size = 15), 
        panel.background = element_rect(fill = NA), legend.position = "none", legend.direction = "horizontal", 
        plot.caption = element_text(face = "italic"), axis.text = element_text(size=15)) +
  labs(x = NULL, y = "Mean plasma glucose\nunder OGTT (mmol/L)")

pdf('SVboxplots.pdf', width = 10, height = 6)
ggarrange(p1, p2, ncol=2)
dev.off()
