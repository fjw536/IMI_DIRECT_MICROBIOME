##################################################################################################
#################### baseline metabolomics vs metaG stratifications ##############################
##################################################################################################

## set directory
setwd('/home/scratch/margar')

## libraries
library(ggplot2); library(tidyverse)

## data upload
mets <- readRDS('data/metabolomics/MetabolitesResiduals.rds')
mtdt.bsl <- read.table('data/curatedmetadata/WP2.1_BL_MetadataWithInfl.txt', header =T, sep = '\t')
annotation <- read.table('data/metabolomics/Chemical_Annotation_Metabolone_DTU_sh.txt', row.names = 1, header=T)
## association to gene richness tertiles
mets.rich <- data.frame(mets, 'GR' = mtdt.bsl$DSGARLTertileGroup)

mets.rich %>% 
  group_by(GR) %>% 
  summarise(across(where(is.numeric), list(mean=mean)))  %>% 
  pivot_longer(ends_with(c('mean', 'sd')),
               names_to = c('name', '.value'),
               names_sep='_') %>% 
  pivot_wider()
  
anova.results <- sapply(mets.rich[,1:973], function(a){
  aov.model <- lm(a~GR, data=mets.rich)
  p.val <- anova(aov.model)$`Pr(>F)`[1]
  return(p.val)
})
anova.fdr <- p.adjust(anova.results, method ='fdr')

pairwise.results <- lapply(mets.rich[,1:973], function(a){
  aov.model <- aov(a~GR, data=mets.rich)
  pairwise.p <- TukeyHSD(aov.model)$GR[,4] # order HGCvsEGC, LGC-EGC, LGC-HGC
  return(pairwise.p)
})

pairwise.results.df <- do.call(rbind.data.frame, pairwise.results)
names(pairwise.results.df) <- c('HGCvsEGC', 'LGCvsEGC', 'LGCvsHGC')
row.names(pairwise.results.df) <- names(pairwise.results)
fc.df <- mets.rich %>% 
  group_by(GR) %>% 
  summarise(across(where(is.numeric), list(mean=mean))) %>% 
  pivot_longer(ends_with(c('mean'))) %>% 
  pivot_wider(names_from = GR, values_from=value) 
fc.df <- as.data.frame(fc.df)
fc.df$name <- gsub('_mean', "", fc.df$name); fc.df$name <- gsub('X', '', fc.df$name)
fc.df$name == row.names(annotation)

mets.rich$GR <- factor(mets.rich$GR, levels=c('LGC', 'HGC', 'EGC'))
delta.lgcegc <- lapply(mets.rich[mets.rich$GR!='HGC',1:973], function(h){
  deltas <- effsize::cliff.delta(h, droplevels(mets.rich[mets.rich$GR!='HGC','GR']))
  delta.res <- data.frame('delta.lgcegc' = deltas$estimate, 
                          'down.lgcegc' = deltas$conf.int[[1]],
                          'up.lgcegc' = deltas$conf.int[[2]])
  return(delta.res)
})
delta.lgcegc <- do.call(rbind.data.frame, delta.lgcegc)

delta.hgcegc <- lapply(mets.rich[mets.rich$GR!='EGC',1:973], function(h){
  deltas <- effsize::cliff.delta(h, droplevels(mets.rich[mets.rich$GR!='EGC','GR']))
  delta.res <- data.frame('delta.hgcegc' = deltas$estimate, 
                          'down.hgcegc' = deltas$conf.int[[1]],
                          'up.hgcegc' = deltas$conf.int[[2]])
  return(delta.res)
})
delta.hgcegc <- do.call(rbind.data.frame, delta.hgcegc)

delta.lgchgc <- lapply(mets.rich[mets.rich$GR!='EGC',1:973], function(h){
  deltas <- effsize::cliff.delta(h, droplevels(mets.rich[mets.rich$GR!='EGC','GR']))
  delta.res <- data.frame('delta.lghegc' = deltas$estimate, 
                          'down.lgchgc' = deltas$conf.int[[1]],
                          'up.lgchgc' = deltas$conf.int[[2]])
  return(delta.res)
})
delta.lgchgc <- do.call(rbind.data.frame, delta.lgchgc)


export.sig <- data.frame(annotation[,4:7], annotation[,14:16], fc.df[,-1], delta.lgcegc, delta.lgchgc, delta.hgcegc, 'anova.fdr' = anova.fdr, pairwise.results.df)
write.table(export.sig, 'data/metabolomics/residualssignificanceGR.txt', sep = '\t')

keep.mets <- names(anova.fdr[anova.fdr<.1]) 
mets.rich.plot <- mets.rich[,names(mets.rich)%in%keep.mets]
annotation.red <- annotation[make.names(row.names(annotation))%in%keep.mets,]
names(mets.rich.plot) <- c(make.names(annotation.red$SHORT_NAME))
mets.rich.plot$GR <- mtdt.bsl$DSGARLTertileGroup

mets.rich.plot.m <- reshape2::melt(mets.rich.plot)
mets.rich.plot.m$family <- rep(annotation.red$SUPER_PATHWAY, each = 775)
mets.rich.plot.m$pathway <- rep(annotation.red$SUB_PATHWAY, each = 775)
mets.rich.plot.m$name <- rep(annotation.red$SHORT_NAME, each = 775)

mets.rich.plot.m$name <- gsub('2-hydroxybutyrate/2-hydroxyisobutyrate',  '2-hydroxybutyrate\n/2-hydroxyisobutyrate', mets.rich.plot.m$name)
mets.rich.plot.m$name <- gsub('palmitoyl-sphingosine-phosphoethanolamine',  'palmitoyl-sphingosine-\nphosphoethanolamine', mets.rich.plot.m$name)
mets.rich.plot.m$name <- gsub('1-\\(1-enyl-palmitoyl)-2-',  '1-\\(1-enyl-palmitoyl\\)-2-\n', mets.rich.plot.m$name)
mets.rich.plot.m$name <- gsub('1-oleoyl-2-docosahexaenoyl-',  '1-oleoyl-2-docosahexaenoyl-\n', mets.rich.plot.m$name)
mets.rich.plot.m$name <- gsub('1-stearoyl-2-docosahexaenoyl-',  '1-stearoyl-2-docosahexaenoyl-\n', mets.rich.plot.m$name)
mets.rich.plot.m$name <- gsub('glycosyl-N-palmitoyl-',  'glycosyl-N-palmitoyl-\n', mets.rich.plot.m$name)
mets.rich.plot.m$name <- gsub('lactosyl-N-palmitoyl-',  'lactosyl-N-palmitoyl-\n', mets.rich.plot.m$name)
mets.rich.plot.m$name <- gsub('linoleoyl-N-arachidnoyl-',  'linoleoyl-N-arachidnoyl-\n', mets.rich.plot.m$name)
mets.rich.plot.m$name <- gsub('sphingomyelin',  'sphingomyelin\n', mets.rich.plot.m$name)
mets.rich.plot.m$name <- gsub('glycerol',  'glycerol\n', mets.rich.plot.m$name)
mets.rich.plot.m$name <- gsub('glycosyl ceramide',  'glycosyl ceramide\n', mets.rich.plot.m$name)

mets.rich.plot.m$GR <- factor(mets.rich.plot.m$GR, levels = c('LGC', 'HGC', 'EGC'))

mets.rich.plot.m %>% 
  filter(family!='Lipid'&family!='Amino Acid') %>% 
  ggplot(aes(name, value, fill=GR)) + 
  geom_boxplot(position='dodge', outlier.size=.4) + facet_wrap(~family, scales = 'free', ncol=4) +
  scale_x_discrete(labels = scales::label_wrap(10)) +
  scale_fill_manual(values=c('gold2', 'gray50', 'blue')) +
  theme(panel.grid.major = element_line(colour = "gray85", linetype = "dashed"), 
        axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90), 
        legend.text = element_text(face = "bold"), legend.title = element_text(face = "bold"), panel.background = element_rect(fill = NA), 
        legend.key = element_rect(fill = NA), legend.background = element_rect(fill = NA), panel.border = element_rect(fill=NA),
        strip.background = element_rect(fill=NA), strip.text.x = element_text(face='bold')) +
  labs(x = NULL, y = "Residuals", fill = "Gene richness tertiles")

bxs.lipids <- mets.rich.plot.m %>% 
  filter(family=='Lipid') %>% 
  ggplot(aes(name, value, fill=GR)) + 
  geom_boxplot(position='dodge', outlier.size = .5) + facet_wrap(~pathway, scales = 'free_x', nrow = 1) + 
  scale_fill_manual(values=c('gold2', 'gray50', 'blue')) +
  theme(panel.grid.major = element_line(colour = "gray85", linetype = "dashed"), 
        axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90), 
        legend.text = element_text(face = "bold"), legend.title = element_text(face = "bold"), panel.background = element_rect(fill = NA), 
        legend.key = element_rect(fill = NA), legend.background = element_rect(fill = NA), panel.border = element_rect(fill=NA),
        strip.background = element_rect(fill=NA), strip.text.x = element_text(face='bold'), legend.position = 'none') +
  labs(x = NULL, y = "Residuals", fill = "Gene richness tertiles", title = 'Lipids') + ylim(-6, 5)
gp.lipids <- ggplotGrob(bxs.lipids)
# gtable::gtable_show_layout(gp.lipids)
facet.columns <- gp.lipids$layout$l[grepl("panel", gp.lipids$layout$name)]
x.var <- sapply(ggplot_build(bxs.lipids)$layout$panel_scales_x,
                function(l) length(l$range$range))
gp.lipids$widths[facet.columns] <- gp.lipids$widths[facet.columns] * x.var
grid::grid.draw(gp.lipids)

bxs.amino <- mets.rich.plot.m %>% 
  filter(family=='Amino Acid') %>% 
  ggplot(aes(name, value, fill=GR)) + 
  geom_boxplot(position='dodge', outlier.size = .5) + facet_wrap(~pathway, nrow = 1, scales = 'free_x') +
  scale_fill_manual(values=c('gold2', 'gray50', 'blue')) +
  theme(panel.grid.major = element_line(colour = "gray85", linetype = "dashed"), 
        axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90), 
        legend.text = element_text(face = "bold"), legend.title = element_text(face = "bold"), panel.background = element_rect(fill = NA), 
        legend.key = element_rect(fill = NA), legend.background = element_rect(fill = NA), panel.border = element_rect(fill=NA),
        strip.background = element_rect(fill=NA), strip.text.x = element_text(face='bold', margin=margin(b=5)), legend.position = 'none') +
  labs(x = NULL, y = "Residuals", fill = "Gene richness tertiles", title = 'Amino Acids')
gp.amino<- ggplotGrob(bxs.amino)
# gtable::gtable_show_layout(gp.amino)
facet.columns <- gp.amino$layout$l[grepl("panel", gp.amino$layout$name)]
x.var <- sapply(ggplot_build(bxs.amino)$layout$panel_scales_x,
                function(l) length(l$range$range))
gp.amino$widths[facet.columns] <- gp.amino$widths[facet.columns] * x.var
grid::grid.draw(gp.amino)

## enterotype differences -----
mets.ent <- data.frame(mets, 'entero' = mtdt.bsl$enterotype)

mets.ent %>% 
  group_by(entero) %>% 
  summarise(across(where(is.numeric), list(mean=mean)))  %>% 
  pivot_longer(ends_with(c('mean', 'sd')),
               names_to = c('name', '.value'),
               names_sep='_') %>% 
  pivot_wider()

anova.results <- sapply(mets.ent[,1:973], function(a){
  aov.model <- lm(a~entero, data=mets.ent)
  p.val <- anova(aov.model)$`Pr(>F)`[1]
  return(p.val)
})
anova.fdr <- p.adjust(anova.results, method ='fdr')

pairwise.results <- lapply(mets.ent[,1:973], function(a){
  aov.model <- aov(a~entero, data=mets.ent)
  pairwise.p <- TukeyHSD(aov.model)$entero[,4] # order B2-B1, P-B1, R-B1, P-B2, R-B2, R-P
  return(pairwise.p)
})

pairwise.results.df <- do.call(rbind.data.frame, pairwise.results)
names(pairwise.results.df) <- c('Bact2vsBact1', 'PrevotellavsBact1', 'RumvsBact1', 'PrevvcBact2', 'RumvsBact2', 'RumvsPrev')
row.names(pairwise.results.df) <- names(pairwise.results)

mets.ent$entero2 <- plyr::revalue(mets.ent$entero, c('Bact1' = 'NoBact2', 
                                                     'Prev' = 'NoBact2',
                                                     'Rum' = 'NoBact2'))


delta.bact2 <- lapply(mets.ent[,1:973], function(h){
  deltas <- effsize::cliff.delta(h, mets.ent$entero2)
  delta.res <- data.frame('delta.bact2' = deltas$estimate, 
                          'down.bact2' = deltas$conf.int[[1]],
                          'up.bact2' = deltas$conf.int[[2]])
  return(delta.res)
})
delta.bact2 <- do.call(rbind.data.frame, delta.bact2)
pairwise.results2<- sapply(mets.ent[,1:973], function(a){
  t.test(a~entero2, data=mets.ent)$p.value
})
pairwise.results2 <- data.frame('pval' = pairwise.results2, 'fdr' = p.adjust(pairwise.results2, method='fdr'))

fc.df <- mets.ent %>% 
  group_by(entero) %>% 
  summarise(across(where(is.numeric), list(mean=mean))) %>% 
  pivot_longer(ends_with(c('mean'))) %>% 
  pivot_wider(names_from = entero, values_from=value) 
fc.df2 <-   mets.ent %>% 
  group_by(entero2) %>% 
  summarise(across(where(is.numeric), list(mean=mean))) %>% 
  pivot_longer(ends_with(c('mean'))) %>% 
  pivot_wider(names_from = entero2, values_from=value)

fc.df <- as.data.frame(fc.df)
fc.df2 <- as.data.frame(fc.df2)
fc.df$name <- gsub('_mean', "", fc.df$name); fc.df$name <- gsub('X', '', fc.df$name)
fc.df2$name <- gsub('_mean', "", fc.df2$name); fc.df2$name <- gsub('X', '', fc.df2$name)
fc.df$name == row.names(annotation)
export.sig <- data.frame(annotation[,4:7], annotation[,14:16], fc.df[,-1], fc.df2[,-1], delta.bact2, 'anova.fdr' = anova.fdr, pairwise.results.df, pairwise.results2)
write.table(export.sig, 'data/metabolomics/residualssignificanceentero.txt', sep = '\t')
