setwd('~/Desktop')

models <- read.table('modelssummaryfigure.txt', header=T, sep='\t')

library(tidyverse)
models$Model <- factor(models$Model, levels=c('Host+Microbiome+Metabolome', 'Host+Microbiome', 'Host'))
models <- models %>% 
  filter(Variable=='Fasting plasma glucose (48m) (mmol/L)'|Variable=='Delta Fasting P. Glucose (mmol/L)'|Variable=='Delta Fasting P. Insulin (pmol/L)')

models$Variable <- factor(models$Variable, levels=c('Delta Fasting P. Insulin (pmol/L)', 'Delta Fasting P. Glucose (mmol/L)', 'Fasting plasma glucose (48m) (mmol/L)'))
pdf('ExplVars_Test.pdf', height = 10, width = 15)
ggplot() +
  geom_col(data=models[models$type=='Train',], aes(Variable, var.explained*100, fill = Model)) +#, position = 'dodge') +
  geom_col(data=models[models$type=='Test',], aes(Variable, var.explained*100, fill = Model), position = 'dodge') +
  scale_fill_manual(values=c('darkblue', 'red', 'gray55')) +
  # scale_color_manual(values = c('purple', 'green')) + 
  facet_wrap(~Approach, scales ='free') +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill=NA),
        panel.grid.major = element_line(linetype = 'dashed', color = 'gray85'),
        strip.background = element_rect(fill=NA),
        strip.text = element_text(size=12.5)) + 
  ylim(0,40) +
  coord_flip() + labs(x= 'Progression of the phenotypic variables over 4 years', y = 'Variance explained (%)')
dev.off()


feats <- read.table('featureimportances_finalmodels.txt', header=T, sep='\t')
tax <- read.table('annotation_modelsfeats.txt', header=T, sep ='\t')
feat.names <- c()
for (i in seq(1:nrow(feats))){
  a <- tax[which(tax$New_Name == feats$Feature[i]),]$Description 
  if (is.null(a)==T){
    feat.names[i] <- feats$Feature[i]
  } else {
    feat.names[i] <- paste0(feats$Feature[i], ':', a)
  }
}
feats$newname <- feat.names

feats %>% 
  mutate("Type" = factor(Type, levels = c('Host', 'Metabolome', 'Microbiome'))) %>% 
  arrange(desc(Type), desc(Feature)) %>% 
  mutate("newname" = factor(newname, levels = unique(newname))) %>% 
  ggplot(aes(Feature.Importance, newname, fill = Type, color = Type)) +
    geom_point() +
    facet_wrap(~Variable)  



feats2 <- split(feats, feats$Variable)
feats2 <- lapply(feats2, function(x){
  x <- x[order(-x$Feature.Importance),]
  x[1:25,]
})

feats2 <- do.call(rbind.data.frame, feats2)
pdf('FIGURE5_feat.imp.pdf', width = 14, height = 8)
feats2 %>% 
  mutate("Type" = factor(Type, levels = c('Host', 'Metabolome', 'Microbiome'))) %>% 
  arrange(desc(Type), desc(Feature)) %>% 
  mutate("newname" = factor(newname, levels = unique(newname))) %>% 
  ggplot(aes(Feature.Importance, newname, fill = Type, color = Type, shape = Type)) +
    geom_point(size=3) +
    geom_segment(aes(x = 0, xend = Feature.Importance, y = newname, yend = newname)) +
    facet_wrap(~Variable, scales = 'free_x') +
    scale_color_manual(values = c('gray50', 'darkblue', 'red')) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5),
          panel.grid.major = element_line(linetype = 'dashed')) +
    labs(x = 'RF feature importance (reduction in tree impurity)', 
         y = 'Feature')
dev.off()
