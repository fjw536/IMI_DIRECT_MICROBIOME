################################################################################################################
################################## Antibiotic resistance genes prevalence ######################################
################################################################################################################

## directory
setwd('/home/scratch/margar')

## libraries
library(ggplot2); library(ggpubr); library(ggbeeswarm); library(vegan); library(factoextra); library(FactoMineR); library(tidyverse)
library(effsize); library(VennDiagram); library(phyloseq)
## data
# amr.genes <- read.table('/home/scratch/vogt/AMR/blast_IGC-CARD.txt', sep='\t')
# load('/home/Data/Repository/Microbiome/WP/Baseline/03.processed/2018_09_11/output.1592samples.genecounts.RData')
# total.genes <- as.data.frame(output.1592samples.genecounts)
# rm(output.1592samples.genecounts)
# 
# total.genes2 <- total.genes[,names(total.genes)%in%amr.genes$V1]
# rm(total.genes)
# saveRDS(total.genes2, 'data/AMRGenes_prediab.rds')
total.genes2 <- readRDS('data/AMRGenes_prediab.rds')
mtdt.bsl <- read.table('data/curatedmetadata/WP2.1_BL_MetadataWithInfl.txt', header=T, sep='\t')
amr.genes <- read.table('data/amr_annotation.txt.csv', sep='\t')
amr.genes <- amr.genes[amr.genes$V1%in%names(total.genes2),]
amr.genes %>% 
  group_by(V3) %>% 
  tally() %>% 
  arrange(-n) %>% 
  ggplot(aes(V3, n)) + geom_col() + 
    theme(axis.text.x = element_text(angle=90))

## subset
row.names(total.genes2) <- gsub('X', '', row.names(total.genes2))
total.genes2 <- total.genes2[row.names(total.genes2)%in%mtdt.bsl$SampleID,]
total.genes2 <- total.genes2[match(mtdt.bsl$SampleID, row.names(total.genes2)),]

#remove duplicates in gene annotations
amr.genes <- amr.genes %>% 
  distinct(V1, .keep_all=T) %>% 

# cluster to common annotations
amr.genes <- amr.genes[amr.genes$V1%in%names(total.genes2),]
amr.genes$V1 == names(total.genes2)
total.genes2 <- as.data.frame(t(total.genes2))
total.genes2 <- total.genes2[match(amr.genes$V1, row.names(total.genes2)),]
total.genes2 <- total.genes2 %>% 
  mutate('gene.name' = amr.genes$V3) %>% 
  group_by(gene.name) %>% 
  summarise(across(where(is.numeric), sum))

total.genes2 <- as.data.frame(total.genes2)
row.names(total.genes2) <- total.genes2$gene.name
total.genes2$gene.name <- NULL
total.genes2 <- as.data.frame(t(total.genes2))

row.names(total.genes2) == mtdt.bsl$SampleID
## AMR richness vs stratification ----
amr.richness <- data.frame('amr.rich'=rowSums(total.genes2!=0), 
                           'prediab' = mtdt.bsl$Gly.Cat2, 
                           'GR' = factor(mtdt.bsl$DSGARLTertileGroup, levels=c('LGC', 'HGC', 'EGC')), 
                           'richness'=mtdt.bsl$DSGeneAbundRichness, 
                           'entero' = mtdt.bsl$enterotype)

amr.richness %>% 
  group_by(prediab) %>% 
  summarise(across(where(is.numeric), list(mean=mean, sd=sd, median=median)))

amr.richness %>% 
  group_by(GR) %>% 
  summarise(across(where(is.numeric), list(mean=mean, sd=sd, median=median)))

amr.richness %>% 
  group_by(entero) %>% 
  summarise(across(where(is.numeric), list(mean=mean, sd=sd, median=median)))
  
plot.pre <- ggplot(amr.richness, aes(prediab, amr.rich), group = interaction(prediab)) + 
  geom_violin(alpha = .1, aes(fill=prediab, color=prediab)) +
  geom_quasirandom(alpha = .6, dodge.width = .9, aes(color = prediab)) +
  geom_boxplot(width = .1, alpha = .7, outlier.colour = NA, position = position_dodge(.9), fill='white')  +
  stat_compare_means() + scale_color_manual(values=c('red', 'darkred')) +
  scale_fill_manual(values=c('red', 'darkred')) +
  theme(panel.grid.major = element_line(colour = "gray85",  linetype = "dashed"), 
        axis.title = element_text(size = 15, face = "bold"), axis.text.y = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 15, face = "bold"),
        panel.background = element_rect(fill = NA), legend.key = element_rect(fill = NA), 
        legend.background = element_rect(fill = NA), legend.position = 'none', 
        panel.border=element_rect(fill=NA)) +
  labs(y = "AMR richness", x = NULL,  colour = "Prediabetic", fill='Prediabetic')

plot.gr1 <- ggplot(amr.richness, aes(GR, amr.rich), group = interaction(GR)) + 
  geom_violin(alpha = .1, aes(fill=GR, color=GR)) +
  geom_quasirandom(alpha = .6, dodge.width = .9, aes(color = GR)) +
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
  labs(y = "AMR Richness", x = NULL,  colour = "Gene richness", fill='Gene richness') 

plot.gr2 <- ggplot(amr.richness, aes(richness, amr.rich)) + geom_point(aes(color=GR)) +
 stat_smooth(method = 'lm', color='black') +
  stat_cor(method='pearson') + 
  scale_color_manual(values=c('darkblue', 'gray90', 'yellow')) + 
  theme(panel.grid.major = element_line(colour = "gray85",  linetype = "dashed"), 
        axis.title = element_text(size = 15, face = "bold"), axis.text = element_text(size = 12, face = "bold"),
        panel.background = element_rect(fill = NA), legend.key = element_rect(fill = NA), 
        legend.background = element_rect(fill = NA), legend.position = c(0.9, 0.15), 
        panel.border=element_rect(fill=NA)) +
  labs(x = "Gene Richness", y = "AMR richness",  colour = "Gene richness")     

plot.ent <- ggplot(amr.richness, aes(entero, amr.rich), group = interaction(entero)) + 
  geom_violin(alpha = .1, aes(fill=entero, color=entero)) +
  geom_quasirandom(alpha = .6, dodge.width = .9, aes(color = entero)) +
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
        legend.background = element_rect(fill = NA), legend.position = 'none', 
        panel.border=element_rect(fill=NA)) +
  labs(y = "AMR Richness", x = NULL,  colour = "Enterotype", fill='Enterotype')  

cowplot::plot_grid(plot.pre, plot.ent,
                   plot.gr1, plot.gr2,
                  labels = c('A', 'B', 'C', ""), label_size = 14)
## models
data.models <- data.frame(total.genes2, 'prediab' = mtdt.bsl$Gly.Cat2, 'GR' = mtdt.bsl$DSGARLTertileGroup, 'entero' = mtdt.bsl$enterotype, 
                     'Gender' = mtdt.bsl$Gender, 'Age' = mtdt.bsl$Age, 'CenterID' = mtdt.bsl$CenterID)
data.models$entero <- plyr::revalue(factor(data.models$entero), c('Bact1' = 'NoBact2', 'Prev'= 'NoBact2', 'Rum'='NoBact2'))

# prediabetics
result.pred <- list()
for (x in seq(from = 1, to = 279)) {
  model <-
    data.frame('amr' = data.models[, x], data.models[, 280:285])
  mod.a <-
    lm(amr ~ prediab + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['prediabIGT', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['prediabIGT', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$prediab)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(data.models)[x], p.val, delta, delta.low, delta.up)
  result.pred[[x]] <- result
}
pred.lm <- do.call(rbind.data.frame, result.pred)
names(pred.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
pred.lm$fdr <- p.adjust(pred.lm$pval, method = 'fdr')
pred.lm <- pred.lm[!is.na(pred.lm$fdr),]
pred.lm[pred.lm$fdr<.1,]

# enterotypes
result.ent <- list()
for (x in seq(from = 1, to = 279)) {
  model <-
    data.frame('amr' = data.models[, x], data.models[, 280:285])
  mod.a <-
    lm(amr ~ entero + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['enteroBact2', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['enteroBact2', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$entero)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(data.models)[x], p.val, delta, delta.low, delta.up)
  result.ent[[x]] <- result
}
entero.lm <- do.call(rbind.data.frame, result.ent)
names(entero.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
entero.lm$fdr <- p.adjust(entero.lm$pval, method = 'fdr')
entero.lm <- entero.lm[!is.na(entero.lm$fdr),]
dim(entero.lm[entero.lm$fdr<.1,])

# generichness 
data.models.gr <- data.models[data.models$GR!='HGC',]
result.gr <- list()
for (x in seq(from = 1, to = 279)) {
  model <-
    data.frame('amr' = data.models.gr[, x], data.models.gr[, 280:285])
  mod.a <-
    lm(amr ~ GR + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['GRLGC', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['GRLGC', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$GR)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(data.models.gr)[x], p.val, delta, delta.low, delta.up)
  result.gr[[x]] <- result
}
GR.lm <- do.call(rbind.data.frame, result.gr)
names(GR.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
GR.lm$fdr <- p.adjust(GR.lm$pval, method = 'fdr')
GR.lm <- GR.lm[!is.na(GR.lm$fdr),]
GR.lm[GR.lm$fdr<.1,]

# cross data between enterotypes and GR
delta.pos.bact2 <- entero.lm[entero.lm$fdr<.1 & entero.lm$delta>0,]$amr
delta.pos.lgc <- GR.lm[GR.lm$fdr<.1 & GR.lm$delta>0,]$amr
plt <- venn.diagram(x=list(delta.pos.bact2, delta.pos.lgc),
                    category.names = c('Bact2', 'Low gene richness'),
                    filename = NULL,
                    cex = 1,
                    cat.cex = 1,
                    lwd = 2,
)
grid::grid.draw(plt)

delta.neg.bact2 <- entero.lm[entero.lm$fdr<.1 & entero.lm$delta<0,]$amr
delta.neg.lgc <- GR.lm[GR.lm$fdr<.1 & GR.lm$delta<0,]$amr
plt <- venn.diagram(x=list(delta.neg.bact2, delta.neg.lgc),
                    category.names = c('Not-Bact2', 'High gene richness'),
                    filename = NULL,
                    cex = 1,
                    cat.cex = 1,
                    lwd = 2,
)
grid::grid.draw(plt)

plt <- venn.diagram(x=list(delta.pos.bact2, delta.pos.lgc, delta.neg.bact2, delta.neg.lgc),
                    category.names = c('Bact2 - delta positive', 'Low gene richness - delta positive', 
                                       'Bact2 - delta negative', 'Low gene richness - delta negative'),
                    filename = NULL,
                    cex = 1,
                    cat.cex = 1,
                    lwd = 2,
)
grid::grid.draw(plt)

list.of.assoc <- list('Bact2.pos' = delta.pos.bact2, 
                     'LGC.pos' = delta.pos.lgc, 
                     'Bact2.neg' = delta.neg.bact2, 
                     'LGC.neg' = delta.neg.lgc)
upset(fromList(list.of.assoc), order.by = 'freq', text.scale=2)

## subset for progression ----
mtdt.all <- read.table('data/curatedmetadata/WP2.1_BlFw-641_MetadataImputed_OnlyBIOCHEM_2.txt', header=T, sep='\t')
mtdt.all[mtdt.all$timepoint=='Baseline',]$StudyID == mtdt.all[mtdt.all$timepoint!='Baseline',]$StudyID
total.genes3 <- total.genes2[row.names(total.genes2)%in%mtdt.all$SampleID_baseline,]
mtdt.all$progression <- rep(mtdt.all[mtdt.all$timepoint=='Follow-up',]$Gly.Cat)
mtdt.all <- mtdt.all[mtdt.all$timepoint=='Baseline',]  

amr.richness2 <- data.frame('amr.rich'=rowSums(total.genes3!=0), 'progression' = factor(mtdt.all$progression, levels=c('NGR', 'IGR', 'T2D')))
ggplot(amr.richness2, aes(progression, amr.rich), group = interaction(progression)) + 
  geom_violin(alpha = .1, aes(fill=progression, color=progression)) +
  geom_quasirandom(alpha = .6, dodge.width = .9, aes(color = progression)) +
  geom_boxplot(width = .1, alpha = .7, outlier.colour = NA, position = position_dodge(.9), fill='white')  +
  stat_compare_means() + scale_color_manual(values=c('darkolivegreen4', 'red', 'darkorchid')) +
  scale_fill_manual(values=c('darkolivegreen4', 'red', 'darkorchid')) +
  stat_compare_means(comparisons = list(c('NGR', 'IGR'),
                                        c('NGR', 'T2D'),
                                        c('IGR', 'T2D'))) +
  theme(panel.grid.major = element_line(colour = "gray85",  linetype = "dashed"), 
        axis.title = element_text(size = 15, face = "bold"), axis.text.y = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 15, face = "bold"),
        panel.background = element_rect(fill = NA), legend.key = element_rect(fill = NA), 
        legend.background = element_rect(fill = NA), legend.position = 'none', 
        panel.border=element_rect(fill=NA)) +
  labs(y = "AMR Richness", x = NULL,  colour = "Progression", fill='Progression')

amr.richness2 %>% 
  group_by(progression) %>% 
  summarise(across(where(is.numeric), list(mean=mean, sd=sd, median=median)))

data.models2 <- data.frame(total.genes3, 'progression' = mtdt.all$progression, 
                          'Gender' = mtdt.all$Gender, 'Age' = mtdt.all$Age, 'CenterID' = mtdt.all$CenterID)
data.models2$progression <- plyr::revalue(factor(data.models2$progression), c('NGR' = 'noT2D', 'IGR'='noT2D'))

# model
result.prg <- list()
for (x in seq(from = 1, to = 279)) {
  model <-
    data.frame('amr' = data.models2[, x], data.models2[, 280:283])
  mod.a <-
    lm(amr ~ progression + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['progressionT2D', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['progressionT2D', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$progression)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(data.models.gr)[x], p.val, delta, delta.low, delta.up)
  result.prg[[x]] <- result
}
progression.lm <- do.call(rbind.data.frame, result.prg)
names(progression.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
progression.lm$fdr <- p.adjust(progression.lm$pval, method = 'fdr')
progression.lm <- progression.lm[!is.na(progression.lm$fdr),]
t2d.assoc <- progression.lm[progression.lm$fdr<.1,]
t2d.assoc$amr <- gsub('X', '', t2d.assoc$amr)

t2d.assoc$gene <- sapply(t2d.assoc$amr, function(a){
  amr.genes[grep(a, amr.genes$V1),]$V3
})

ggplot(t2d.assoc, aes(amr, as.numeric(delta))) + geom_col() + coord_flip() +
  geom_vline(xintercept = 0) +
  theme(axis.text.x = element_text(size=15), axis.text.y=element_text(size=10)) +
  labs(x=NULL, y = "Cliff's Effect size")

data.models2 <- data.frame(total.genes3, 'progression' = mtdt.all$progression, 
                           'Gender' = mtdt.all$Gender, 'Age' = mtdt.all$Age, 'CenterID' = mtdt.all$CenterID)
data.models2$progression <- plyr::revalue(factor(data.models2$progression), c('T2D' = 'noNGR', 'IGR'='noNGR'))
result.ngr <- list()
for (x in seq(from = 1, to = 279)) {
  model <-
    data.frame('amr' = data.models2[, x], data.models2[, 280:283])
  mod.a <-
    lm(amr ~ progression + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['progressionNGR', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['progressionNGR', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$progression)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(data.models.gr)[x], p.val, delta, delta.low, delta.up)
  result.ngr[[x]] <- result
}
ngr.lm <- do.call(rbind.data.frame, result.ngr)
names(ngr.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
ngr.lm$fdr <- p.adjust(ngr.lm$pval, method = 'fdr')
ngr.lm <- ngr.lm[!is.na(ngr.lm$fdr),]
ngr.assoc <- ngr.lm[ngr.lm$fdr<.1,]

## beta diversity ----
# bray-curtis
total.genes.dist <- vegdist(total.genes2, method='bray')
total.genes.pcoa <- prcomp(total.genes.dist, center=T, scale=T)
fviz_screeplot(total.genes.pcoa)
fviz_pca_ind(total.genes.pcoa, geom = 'point', col.ind = mtdt.bsl$Gly.Cat2)
pairwise.adonis2(total.genes2, as.factor(mtdt.bsl$Gly.Cat2))
fviz_pca_ind(total.genes.pcoa, geom = 'point', col.ind = mtdt.bsl$DSGARLTertileGroup)
pairwise.adonis(total.genes2, as.factor(mtdt.bsl$DSGARLTertileGroup))
fviz_pca_ind(total.genes.pcoa, geom = 'point', col.ind = mtdt.bsl$enterotype, addEllipses = F)
pairwise.adonis2(total.genes2, as.factor(mtdt.bsl$enterotype))

total.genes.dist.red <- vegdist(total.genes3, method='bray')
total.genes.pcoa.red <- prcomp(total.genes.dist.red, center=T, scale=T)
fviz_screeplot(total.genes.pcoa.red)
fviz_pca_ind(total.genes.pcoa.red, geom = 'point', col.ind = mtdt.all$progression)
adonis(total.genes.dist.red~mtdt.all$progression)
pairwise.adonis(total.genes3, as.factor(mtdt.all$progression))

## phenotypic models ----
progression.model <- cbind(data.models, mtdt.bsl)
progression.model$Smoking.Status.BL <- plyr::revalue(factor(progression.model$Smoking.Status.BL), c('ex-smoker'='non-smoker', 'never'='non-smoker'))

# amr vs smoking
result.smoking <- list()
for (x in seq(from = 1, to = 279)) {
  print(names(progression.model[x]))
  model <-
    data.frame('amr' = progression.model[, x], progression.model[, 280:357])
  mod.a <-
    lm(amr ~ Smoking.Status.BL + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Smoking.Status.BLnon-smoker', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Smoking.Status.BLnon-smoker', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$Smoking.Status.BL)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.smoking[[x]] <- result
}
smoking.lm <- do.call(rbind.data.frame, result.smoking)
names(smoking.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
smoking.lm$fdr <- p.adjust(smoking.lm$pval, method = 'fdr')


# amr vs Height
result.height <- list()
for (x in seq(from = 1, to = 279)) {
  print(names(progression.model[x]))
  # result.model$amr[x] <- names(amr.model.bsl[x])
  model <-
    data.frame('amr' = progression.model[, x], progression.model[, 280:357])
  mod.a <-
    lm(amr ~ Height.cm + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Height.cm', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Height.cm', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$Height.cm)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.height[[x]] <- result
}
height.lm <- do.call(rbind.data.frame, result.height)
names(height.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
height.lm$fdr <- p.adjust(height.lm$pval, method = 'fdr')

# amr vs Weight
result.weight <- list()
for (x in seq(from = 1, to = 279)) {
  print(names(progression.model[x]))
  model <-
    data.frame('amr' = progression.model[, x], progression.model[, 280:357])
  mod.a <-
    lm(amr ~ Weight.kg + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Weight.kg', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Weight.kg', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$Weight.kg)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.weight[[x]] <- result
}
weight.lm <- do.call(rbind.data.frame, result.weight)
names(weight.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
weight.lm$fdr <- p.adjust(weight.lm$pval, method = 'fdr')

# amr vs Waist
result.waist <- list()
for (x in seq(from = 1, to = 279)) {
  print(names(progression.model[x]))
  model <-
    data.frame('amr' = progression.model[, x], progression.model[, 280:357])
  mod.a <-
    lm(amr ~ Waist.cm + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Waist.cm', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Waist.cm', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$Waist.cm)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.waist[[x]] <- result
}
waist.lm <- do.call(rbind.data.frame, result.waist)
names(waist.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
waist.lm$fdr <- p.adjust(waist.lm$pval, method = 'fdr')

# amr vs hip
result.hip <- list()
for (x in seq(from = 1, to = 279)) {
  print(names(progression.model[x]))
  model <-
    data.frame('amr' = progression.model[, x], progression.model[, 280:357])
  mod.a <-
    lm(amr ~ Hip.cm + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Hip.cm', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Hip.cm', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$Hip.cm)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.hip[[x]] <- result
}
hip.lm <- do.call(rbind.data.frame, result.hip)
names(hip.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
hip.lm$fdr <- p.adjust(hip.lm$pval, method = 'fdr')

# amr vs Waist:Hip
result.waist.hip <- list()
for (x in seq(from = 1, to = 279)) {
  print(names(progression.model[x]))
  model <-
    data.frame('amr' = progression.model[, x], progression.model[, 280:357])
  mod.a <-
    lm(amr ~ Waist.Hip + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Waist.Hip', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Waist.Hip', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$Waist.Hip)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.waist.hip[[x]] <- result
}
waist.hip.lm <- do.call(rbind.data.frame, result.waist.hip)
names(waist.hip.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
waist.hip.lm$fdr <- p.adjust(waist.hip.lm$pval, method = 'fdr')

# amr vs BMI
result.bmi <- list()
for (x in seq(from = 1, to = 279)) {
  print(names(progression.model[x]))
  model <-
    data.frame('amr' = progression.model[, x], progression.model[, 280:357])
  mod.a <-
    lm(amr ~ BMI + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['BMI', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['BMI', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$BMI)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.bmi[[x]] <- result
}
bmi.lm <- do.call(rbind.data.frame, result.bmi)
names(bmi.lm) <- c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
bmi.lm$fdr <- p.adjust(bmi.lm$pval, method = 'fdr')

# amr vs Fast.Glucose
result.Fast.Glucose <- list()
for (x in seq(from = 1, to = 279)) {
  print(names(progression.model[x]))
  model <-
    data.frame('amr' = progression.model[, x], progression.model[, 280:357])
  mod.a <-
    lm(amr ~ Fasting.Glucose + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Fasting.Glucose', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Fasting.Glucose', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$Fasting.Glucose)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.Fast.Glucose[[x]] <- result
}
Fast.Glucose.lm <- do.call(rbind.data.frame, result.Fast.Glucose)
names(Fast.Glucose.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
Fast.Glucose.lm$fdr <- p.adjust(Fast.Glucose.lm$pval, method = 'fdr')
Fast.Glucose.lm[Fast.Glucose.lm$fdr < .1, ]

# amr vs Basal.Glucose
result.Basal.Glucose <- list()
for (x in seq(from = 1, to = 279)) {
  print(names(progression.model[x]))
  # result.model$amr[x] <- names(amr.model.bsl[x])
  model <-
    data.frame('amr' = progression.model[, x], progression.model[, 280:357])
  mod.a <-
    lm(amr ~ Basal.Glucose + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Basal.Glucose', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Basal.Glucose', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$Basal.Glucose)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.Basal.Glucose[[x]] <- result
}
Basal.Glucose.lm <- do.call(rbind.data.frame, result.Basal.Glucose)
names(Basal.Glucose.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
Basal.Glucose.lm$fdr <-
  p.adjust(Basal.Glucose.lm$pval, method = 'fdr')
Basal.Glucose.lm[Basal.Glucose.lm$fdr < .1, ]

# amr vs Mean.Glucose
result.Mean.Glucose <- list()
for (x in seq(from = 1, to = 279)) {
  print(names(progression.model[x]))
  # result.model$amr[x] <- names(amr.model.bsl[x])
  model <-
    data.frame('amr' = progression.model[, x], progression.model[, 280:357])
  mod.a <-
    lm(amr ~ Mean.Glucose + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Mean.Glucose', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Mean.Glucose', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$Mean.Glucose)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.Mean.Glucose[[x]] <- result
}
Mean.Glucose.lm <- do.call(rbind.data.frame, result.Mean.Glucose)
names(Mean.Glucose.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
Mean.Glucose.lm$fdr <- p.adjust(Mean.Glucose.lm$pval, method = 'fdr')
Mean.Glucose.lm[Mean.Glucose.lm$fdr < .1, ]

# amr vs Glucose.2h
result.Glucose.2h <- list()
for (x in seq(from = 1, to = 279)) {
  print(names(progression.model[x]))
  # result.model$amr[x] <- names(amr.model.bsl[x])
  model <-
    data.frame('amr' = progression.model[, x], progression.model[, 280:357])
  mod.a <-
    lm(amr ~ glucose.2h + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['glucose.2h', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['glucose.2h', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$glucose.2h)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.Glucose.2h[[x]] <- result
}
Glucose.2h.lm <- do.call(rbind.data.frame, result.Glucose.2h)
names(Glucose.2h.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
Glucose.2h.lm$fdr <- p.adjust(Glucose.2h.lm$pval, method = 'fdr')

# amr vs Fast.Insulin
result.Fast.Insulin <- list()
for (x in seq(from = 1, to = 279)) {
  print(names(progression.model[x]))
  model <-
    data.frame('amr' = progression.model[, x], progression.model[, 280:357])
  mod.a <-
    lm(amr ~ Fasting.Insulin + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Fasting.Insulin', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Fasting.Insulin', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$Fasting.Insulin)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.Fast.Insulin[[x]] <- result
}
Fast.Insulin.lm <- do.call(rbind.data.frame, result.Fast.Insulin)
names(Fast.Insulin.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
Fast.Insulin.lm$fdr <- p.adjust(Fast.Insulin.lm$pval, method = 'fdr')
Fast.Insulin.lm[Fast.Insulin.lm$fdr < .1, ]

# amr vs Basal.Insulin
result.Basal.Insulin <- list()
for (x in seq(from = 1, to = 279)) {
  print(names(progression.model[x]))
  # result.model$amr[x] <- names(amr.model.bsl[x])
  model <-
    data.frame('amr' = progression.model[, x], progression.model[, 280:357])
  mod.a <-
    lm(amr ~ Basal.Insulin + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Basal.Insulin', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Basal.Insulin', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$Basal.Insulin)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.Basal.Insulin[[x]] <- result
}
Basal.Insulin.lm <- do.call(rbind.data.frame, result.Basal.Insulin)
names(Basal.Insulin.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
Basal.Insulin.lm$fdr <-
  p.adjust(Basal.Insulin.lm$pval, method = 'fdr')
Basal.Insulin.lm[Basal.Insulin.lm$fdr < .1, ]

# amr vs Mean.Insulin
result.Mean.Insulin <- list()
for (x in seq(from = 1, to = 279)) {
  print(names(progression.model[x]))
  # result.model$amr[x] <- names(amr.model.bsl[x])
  model <-
    data.frame('amr' = progression.model[, x], progression.model[, 280:357])
  mod.a <-
    lm(amr ~ Mean.Insulin + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Mean.Insulin', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Mean.Insulin', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$Mean.Insulin)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  
  result.Mean.Insulin[[x]] <- result
}
Mean.Insulin.lm <- do.call(rbind.data.frame, result.Mean.Insulin)
names(Mean.Insulin.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
Mean.Insulin.lm$fdr <- p.adjust(Mean.Insulin.lm$pval, method = 'fdr')
Mean.Insulin.lm[Mean.Insulin.lm$fdr < .1, ]

# amr vs Basal.InsSecr
result.Basal.InsSecr <- list()
for (x in seq(from = 1, to = 279)) {
  print(names(progression.model[x]))
  # result.model$amr[x] <- names(amr.model.bsl[x])
  model <-
    data.frame('amr' = progression.model[, x], progression.model[, 280:357])
  mod.a <-
    lm(amr ~ Basal.InsulinSecretion + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Basal.InsulinSecretion', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Basal.InsulinSecretion', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$Basal.InsulinSecretion)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.Basal.InsSecr[[x]] <- result
}
Basal.InsSecr.lm <- do.call(rbind.data.frame, result.Basal.InsSecr)
names(Basal.InsSecr.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
Basal.InsSecr.lm$fdr <-
  p.adjust(Basal.InsSecr.lm$pval, method = 'fdr')
Basal.InsSecr.lm[Basal.InsSecr.lm$fdr < .1, ]

# amr vs Total.InsSecr
result.Total.InsSecr <- list()
for (x in seq(from = 1, to = 279)) {
  print(names(progression.model[x]))
  # result.model$amr[x] <- names(amr.model.bsl[x])
  model <-
    data.frame('amr' = progression.model[, x], progression.model[, 280:357])
  mod.a <-
    lm(amr ~ Total.InsulinSecretion + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Total.InsulinSecretion', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Total.InsulinSecretion', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$Total.InsulinSecretion)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.Total.InsSecr[[x]] <- result
}
Total.InsSecr.lm <- do.call(rbind.data.frame, result.Total.InsSecr)
names(Total.InsSecr.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
Total.InsSecr.lm$fdr <-
  p.adjust(Total.InsSecr.lm$pval, method = 'fdr')
Total.InsSecr.lm[Total.InsSecr.lm$fdr < .1, ]

# amr vs Gluc.Sensitivity
result.Gluc.Sensitivity <- list()
for (x in seq(from = 1, to = 279)) {
  print(names(progression.model[x]))
  # result.model$amr[x] <- names(amr.model.bsl[x])
  model <-
    data.frame('amr' = progression.model[, x], progression.model[, 280:357])
  mod.a <-
    lm(
      amr ~ Glucose.Sensitivity + Gender + CenterID + Age,
      data = model,
      na.action = na.omit
    )
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Glucose.Sensitivity', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Glucose.Sensitivity', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$Glucose.Sensitivity)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.Gluc.Sensitivity[[x]] <- result
}
Gluc.Sensitivity.lm <-
  do.call(rbind.data.frame, result.Gluc.Sensitivity)
names(Gluc.Sensitivity.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
Gluc.Sensitivity.lm$fdr <-
  p.adjust(Gluc.Sensitivity.lm$pval, method = 'fdr')
Gluc.Sensitivity.lm[Gluc.Sensitivity.lm$fdr < .1, ]

# amr vs OGIS.2h
result.OGIS.2h <- list()
for (x in seq(from = 1, to = 279)) {
  print(names(progression.model[x]))
  # result.model$amr[x] <- names(amr.model.bsl[x])
  model <-
    data.frame('amr' = progression.model[, x], progression.model[, 280:357])
  mod.a <-
    lm(amr ~ OGIS.2h + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['OGIS.2h', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['OGIS.2h', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$OGIS.2h)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  
  result.OGIS.2h[[x]] <- result
}
OGIS.2h.lm <- do.call(rbind.data.frame, result.OGIS.2h)
names(OGIS.2h.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
OGIS.2h.lm$fdr <- p.adjust(OGIS.2h.lm$pval, method = 'fdr')
OGIS.2h.lm[OGIS.2h.lm$fdr < .1, ]

# amr vs Stumvoll
result.Stumvoll <- list()
for (x in seq(from = 1, to = 279)) {
  print(names(progression.model[x]))
  # result.model$amr[x] <- names(amr.model.bsl[x])
  model <-
    data.frame('amr' = progression.model[, x], progression.model[, 280:357])
  mod.a <-
    lm(amr ~ Stumvoll + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Stumvoll', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Stumvoll', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$Stumvoll)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  
  result.Stumvoll[[x]] <- result
}
Stumvoll.lm <- do.call(rbind.data.frame, result.Stumvoll)
names(Stumvoll.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
Stumvoll.lm$fdr <- p.adjust(Stumvoll.lm$pval, method = 'fdr')

# amr vs Matsuda
result.Matsuda <- list()
for (x in seq(from = 1, to = 279)) {
  print(names(progression.model[x]))
  # result.model$amr[x] <- names(amr.model.bsl[x])
  model <-
    data.frame('amr' = progression.model[, x], progression.model[, 280:357])
  mod.a <-
    lm(amr ~ Matsuda + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Matsuda', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Matsuda', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$Matsuda)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.Matsuda[[x]] <- result
}
Matsuda.lm <- do.call(rbind.data.frame, result.Matsuda)
names(Matsuda.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
Matsuda.lm$fdr <- p.adjust(Matsuda.lm$pval, method = 'fdr')
Matsuda.lm[Matsuda.lm$fdr < .1, ]

# amr vs Clinsb
result.Clinsb <- list()
for (x in seq(from = 1, to = 279)) {
  print(names(progression.model[x]))
  # result.model$amr[x] <- names(amr.model.bsl[x])
  model <-
    data.frame('amr' = progression.model[, x], progression.model[, 280:357])
  mod.a <-
    lm(amr ~ Clinsb + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Clinsb', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Clinsb', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$Clinsb)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.Clinsb[[x]] <- result
}
Clinsb.lm <- do.call(rbind.data.frame, result.Clinsb)
names(Clinsb.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
Clinsb.lm$fdr <- p.adjust(Clinsb.lm$pval, method = 'fdr')

# amr vs Clins
result.Clins <- list()
for (x in seq(from = 1, to = 279)) {
  print(names(progression.model[x]))
  model <-
    data.frame('amr' = progression.model[, x], progression.model[, 280:357])
  mod.a <-
    lm(amr ~ Clins + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Clins', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Clins', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$Clins)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.Clins[[x]] <- result
}
Clins.lm <- do.call(rbind.data.frame, result.Clins)
names(Clins.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
Clins.lm$fdr <- p.adjust(Clins.lm$pval, method = 'fdr')

# amr vs activeGLP1
result.activeGLP1 <- list()
for (x in seq(from = 1, to = 279)) {
  print(names(progression.model[x]))
  model <-
    data.frame('amr' = progression.model[, x], progression.model[, 280:357])
  mod.a <-
    lm(amr ~ Active.GLP1.Conc.0m + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Active.GLP1.Conc.0m', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Active.GLP1.Conc.0m', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$Active.GLP1.Conc.0m)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.activeGLP1[[x]] <- result
}
activeGLP1.lm <- do.call(rbind.data.frame, result.activeGLP1)
names(activeGLP1.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
activeGLP1.lm$fdr <- p.adjust(activeGLP1.lm$pval, method = 'fdr')

# amr vs totalGLP1.0
result.totalGLP1.0 <- list()
for (x in seq(from = 1, to = 279)) {
  print(names(progression.model[x]))
  model <-
    data.frame('amr' = progression.model[, x], progression.model[, 280:357])
  mod.a <-
    lm(amr ~ Total.GLP1.Conc.0m + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Total.GLP1.Conc.0m', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Total.GLP1.Conc.0m', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$Total.GLP1.Conc.0m)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.totalGLP1.0[[x]] <- result
}
totalGLP1.0.lm <- do.call(rbind.data.frame, result.totalGLP1.0)
names(totalGLP1.0.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
totalGLP1.0.lm$fdr <- p.adjust(totalGLP1.0.lm$pval, method = 'fdr')

# amr vs totalGLP1.0
result.totalGLP1.6 <- list()
for (x in seq(from = 1, to = 279)) {
  print(names(progression.model[x]))
  model <-
    data.frame('amr' = progression.model[, x], progression.model[, 280:357])
  mod.a <-
    lm(amr ~ Total.GLP1.Conc.60m + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Total.GLP1.Conc.60m', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Total.GLP1.Conc.60m', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$Total.GLP1.Conc.60m)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.totalGLP1.60[[x]] <- result
}
totalGLP1.60.lm <- do.call(rbind.data.frame, result.totalGLP1.60)
names(totalGLP1.60.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
totalGLP1.60.lm$fdr <- p.adjust(totalGLP1.60.lm$pval, method = 'fdr')

# amr vs hsCRP
result.hscrp <- list()
for (x in seq(from = 1, to = 279)) {
  print(names(progression.model[x]))
  model <-
    data.frame('amr' = progression.model[, x], progression.model[, 280:357])
  mod.a <-
    lm(amr ~ hsCRP + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['hsCRP', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['hsCRP', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$hsCRP)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.hscrp[[x]] <- result
}
hscrp.lm <- do.call(rbind.data.frame, result.hscrp)
names(hscrp.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
hscrp.lm$fdr <- p.adjust(hscrp.lm$pval, method = 'fdr')

# amr vs Fasting.HDL
result.hdl <- list()
for (x in seq(from = 1, to = 279)) {
  print(names(progression.model[x]))
  model <-
    data.frame('amr' = progression.model[, x], progression.model[, 280:357])
  mod.a <-
    lm(amr ~ Fasting.HDL + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Fasting.HDL', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Fasting.HDL', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$Fasting.HDL)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.hdl[[x]] <- result
}
hdl.lm <- do.call(rbind.data.frame, result.hdl)
names(hdl.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
hdl.lm$fdr <- p.adjust(hdl.lm$pval, method = 'fdr')

# amr vs Fasting.LDL
result.ldl <- list()
for (x in seq(from = 1, to = 279)) {
  print(names(progression.model[x]))
  model <-
    data.frame('amr' = progression.model[, x], progression.model[, 280:357])
  mod.a <-
    lm(amr ~ Fasting.LDL + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Fasting.LDL', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Fasting.LDL', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$Fasting.LDL)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.ldl[[x]] <- result
}
ldl.lm <- do.call(rbind.data.frame, result.ldl)
names(ldl.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
ldl.lm$fdr <- p.adjust(ldl.lm$pval, method = 'fdr')

# amr vs Fasting.TG
result.tg <- list()
for (x in seq(from = 1, to = 279)) {
  print(names(progression.model[x]))
  model <-
    data.frame('amr' = progression.model[, x], progression.model[, 280:357])
  mod.a <-
    lm(amr ~ Fasting.TG + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Fasting.TG', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Fasting.TG', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$Fasting.TG)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.tg[[x]] <- result
}
tg.lm <- do.call(rbind.data.frame, result.tg)
names(tg.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
tg.lm$fdr <- p.adjust(tg.lm$pval, method = 'fdr')

# amr vs Fasting.AST
result.ast <- list()
for (x in seq(from = 1, to = 279)) {
  print(names(progression.model[x]))
  model <-
    data.frame('amr' = progression.model[, x], progression.model[, 280:357])
  mod.a <-
    lm(amr ~ Fasting.AST + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Fasting.AST', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Fasting.AST', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$Fasting.AST)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.ast[[x]] <- result
}
ast.lm <- do.call(rbind.data.frame, result.ast)
names(ast.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
ast.lm$fdr <- p.adjust(ast.lm$pval, method = 'fdr')

# amr vs Fasting.ALT
result.alt <- list()
for (x in seq(from = 1, to = 279)) {
  print(names(progression.model[x]))
  model <-
    data.frame('amr' = progression.model[, x], progression.model[, 280:357])
  mod.a <-
    lm(amr ~ Fasting.ALT + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Fasting.ALT', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Fasting.ALT', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$Fasting.ALT)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.alt[[x]] <- result
}
alt.lm <- do.call(rbind.data.frame, result.alt)
names(alt.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
asl.lm$fdr <- p.adjust(alt.lm$pval, method = 'fdr')

# amr vs Fasting.Chol
result.chol <- list()
for (x in seq(from = 1, to = 279)) {
  print(names(progression.model[x]))
  model <-
    data.frame('amr' = progression.model[, x], progression.model[, 280:357])
  mod.a <-
    lm(amr ~ Fasting.Chol + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Fasting.Chol', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Fasting.Chol', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$Fasting.Chol)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.chol[[x]] <- result
}
chol.lm <- do.call(rbind.data.frame, result.chol)
names(chol.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
chol.lm$fdr <- p.adjust(chol.lm$pval, method = 'fdr')

# amr vs Liver.Iron
result.liver.iron <- list()
for (x in seq(from = 1, to = 279)) {
  print(names(progression.model[x]))
  model <-
    data.frame('amr' = progression.model[, x], progression.model[, 280:357])
  mod.a <-
    lm(amr ~ Liver.Iron + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Liver.Iron', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Liver.Iron', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$Liver.Iron)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.liver.iron[[x]] <- result
}
liver.iron.lm <- do.call(rbind.data.frame, result.liver.iron)
names(liver.iron.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
liver.iron.lm$fdr <- p.adjust(liver.iron.lm$pval, method = 'fdr')

# amr vs Panc.Iron
result.panc.iron <- list()
for (x in seq(from = 1, to = 279)) {
  print(names(progression.model[x]))
  model <-
    data.frame('amr' = progression.model[, x], progression.model[, 280:357])
  mod.a <-
    lm(amr ~ Panc.Iron + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Panc.Iron', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Panc.Iron', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$Panc.Iron)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.panc.iron[[x]] <- result
}
panc.iron.lm <- do.call(rbind.data.frame, result.panc.iron)
names(panc.iron.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
panc.iron.lm$fdr <- p.adjust(panc.iron.lm$pval, method = 'fdr')

# amr vs IAAT
result.IAAT <- list()
for (x in seq(from = 1, to = 279)) {
  print(names(progression.model[x]))
  model <-
    data.frame('amr' = progression.model[, x], progression.model[, 280:357])
  mod.a <-
    lm(amr ~ IAAT + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['IAAT', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['IAAT', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$IAAT)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.IAAT[[x]] <- result
}
IAAT.lm <- do.call(rbind.data.frame, result.IAAT)
names(IAAT.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
IAAT.lm$fdr <- p.adjust(IAAT.lm$pval, method = 'fdr')

# amr vs ASAT
result.ASAT <- list()
for (x in seq(from = 1, to = 279)) {
  print(names(progression.model[x]))
  model <-
    data.frame('amr' = progression.model[, x], progression.model[, 280:357])
  mod.a <-
    lm(amr ~ ASAT + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['ASAT', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['ASAT', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$ASAT)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.ASAT[[x]] <- result
}
ASAT.lm <- do.call(rbind.data.frame, result.ASAT)
names(ASAT.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
ASAT.lm$fdr <- p.adjust(ASAT.lm$pval, method = 'fdr')

# amr vs TAAT
result.TAAT <- list()
for (x in seq(from = 1, to = 279)) {
  print(names(progression.model[x]))
  model <-
    data.frame('amr' = progression.model[, x], progression.model[, 280:357])
  mod.a <-
    lm(amr ~ TAAT + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['TAAT', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['TAAT', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$TAAT)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.TAAT[[x]] <- result
}
TAAT.lm <- do.call(rbind.data.frame, result.TAAT)
names(TAAT.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
TAAT.lm$fdr <- p.adjust(TAAT.lm$pval, method = 'fdr')

# amr vs LiverFat..
result.liverfat <- list()
for (x in seq(from = 1, to = 279)) {
  print(names(progression.model[x]))
  model <-
    data.frame('amr' = progression.model[, x], progression.model[, 280:357])
  mod.a <-
    lm(amr ~ LiverFat.. + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['LiverFat..', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['LiverFat..', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$LiverFat..)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.liverfat[[x]] <- result
}
liverfat.lm <- do.call(rbind.data.frame, result.liverfat)
names(liverfat.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
liverfat.lm$fdr <- p.adjust(liverfat.lm$pval, method = 'fdr')

# amr vs PancreasFat..
result.pancreasfat <- list()
for (x in seq(from = 1, to = 279)) {
  print(names(progression.model[x]))
  model <-
    data.frame('amr' = progression.model[, x], progression.model[, 280:357])
  mod.a <-
    lm(amr ~ PancreasFat.. + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['PancreasFat..', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['PancreasFat..', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$PancreasFat..)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.pancreasfat[[x]] <- result
}
pancreasfat.lm <- do.call(rbind.data.frame, result.pancreasfat)
names(pancreasfat.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
pancreasfat.lm$fdr <- p.adjust(pancreasfat.lm$pval, method = 'fdr')

# amr vs BP.Sys
result.BP.Sys <- list()
for (x in seq(from = 1, to = 279)) {
  print(names(progression.model[x]))
  # result.model$amr[x] <- names(amr.model.bsl[x])
  model <-
    data.frame('amr' = progression.model[, x], progression.model[, 280:357])
  mod.a <-
    lm(amr ~ BP.S.Mean + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['BP.S.Mean', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['BP.S.Mean', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$BP.S.Mean)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  
  result.BP.Sys[[x]] <- result
}
BP.Sys.lm <- do.call(rbind.data.frame, result.BP.Sys)
names(BP.Sys.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
BP.Sys.lm$fdr <- p.adjust(BP.Sys.lm$pval, method = 'fdr')
BP.Sys.lm[BP.Sys.lm$fdr < .1, ]

# amr vs BP.Dias
result.BP.Dias <- list()
for (x in seq(from = 1, to = 279)) {
  print(names(progression.model[x]))
  # result.model$amr[x] <- names(amr.model.bsl[x])
  model <-
    data.frame('amr' = progression.model[, x], progression.model[, 280:357])
  mod.a <-
    lm(amr ~ BP.D.Mean + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['BP.D.Mean', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['BP.D.Mean', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$BP.D.Mean)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.BP.Dias[[x]] <- result
}
BP.Dias.lm <- do.call(rbind.data.frame, result.BP.Dias)
names(BP.Dias.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
BP.Dias.lm$fdr <- p.adjust(BP.Dias.lm$pval, method = 'fdr')
BP.Dias.lm[BP.Dias.lm$fdr < .1, ]

# amr vs Value.vm.hpf.mean.BL
result.physical <- list()
for (x in seq(from = 1, to = 279)) {
  print(names(progression.model[x]))
  # result.model$amr[x] <- names(amr.model.bsl[x])
  model <-
    data.frame('amr' = progression.model[, x], progression.model[, 280:357])
  mod.a <-
    lm(amr ~ Value.vm.hpf.mean.BL + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Value.vm.hpf.mean.BL', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Value.vm.hpf.mean.BL', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$Value.vm.hpf.mean.BL)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.physical[[x]] <- result
}
physical.lm <- do.call(rbind.data.frame, result.physical)
names(physical.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
physical.lm$fdr <- p.adjust(physical.lm$pval, method = 'fdr')

# amr vs HDI
result.hdi<- list()
for (x in seq(from = 1, to = 279)) {
  print(names(progression.model[x]))
  # result.model$amr[x] <- names(amr.model.bsl[x])
  model <-
    data.frame('amr' = progression.model[, x], progression.model[, 280:357])
  mod.a <-
    lm(amr ~ HDI + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['HDI', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['HDI', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$HDI)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.hdi[[x]] <- result
}
hdi.lm <- do.call(rbind.data.frame, result.hdi)
names(hdi.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
hdi.lm$fdr <- p.adjust(hdi.lm$pval, method = 'fdr')

# extract significant ones - amr ----
lm.outputs<-grep(".lm",names(.GlobalEnv),value=TRUE)
results.dfs <- do.call("list",mget(lm.outputs))
results.dfs.sig <- lapply(results.dfs, function(x){
  x <- x[x$fdr<.1,]
  return(data.frame('amr' = x$amr, 'delta' = x$delta, 'fdr' = x$fdr))
})
amr.associations <-  do.call(rbind.data.frame, results.dfs.sig)
amr.associations <- amr.associations[complete.cases(amr.associations),]
amr.associations$association <- gsub('.lm.[0-9][0-9]', '', row.names(amr.associations))
amr.associations$association <- gsub('.lm.[0-9]', '', amr.associations$association)
amr.associations$association <- gsub('.lm', '', amr.associations$association)

amr.associations$delta <- as.numeric(amr.associations$delta)
direction <- c()
for (i in 1:nrow(amr.associations)){
  if(amr.associations[i,]$delta > 0) {
  direction[[i]] <- c("Positive")
} else {
  if(amr.associations[i,]$delta  == 0) {
    direction[[i]] <- c("Zero")
  } else {
    direction[[i]] <- c("Negative")
  }
}}
amr.associations$direction <- unlist(direction)
amr.associations$direction <- factor(amr.associations$direction, levels=c('Positive', 'Negative'))
amr.associations <-amr.associations[amr.associations$association!='GR',]
amr.associations <-amr.associations[amr.associations$association!='entero',]
amr.associations <-amr.associations[amr.associations$association!='pred',]
amr.associations <-amr.associations[amr.associations$association!='progression',]

# amr.associations$name <- paste(amr.associations$amr, amr.associations$spps, sep=':')
# 
# pdf('amrplot.pdf')
ggplot(amr.associations, aes(association, amr, fill=delta, shape=direction)) +
  geom_point(size=3.9)  + scale_shape_manual(values = c(24, 25)) +
  scale_fill_gradient2(low='blue', high = 'red') +
  theme(panel.border = element_rect(fill=NA), panel.background = element_blank(),
        panel.grid.major = element_line(color='gray85', linetype='dashed'),
        axis.text.x = element_text(angle=90, vjust=0.5, face='bold', color='black', size=12),
        axis.text.y = element_text(face='bold.italic'),
        axis.title = element_text(size = 13, face = "bold"), axis.text = element_text(face = "bold"),
        legend.text = element_text(face = "bold"), legend.title = element_text(face = "bold"),
        legend.background = element_blank(), legend.box.background = element_blank(), legend.key = element_blank()) +
  labs(x = NULL, y = NULL, fill = "Cliff's delta", shape = "Direction", title = 'amr associations')

## enterotype model controlling by AMR richness ----
data.models$amr.rich <- amr.richness$amr.rich
result.ent2 <- list()
for (x in seq(from = 1, to = 279)) {
  model <-
    data.frame('amr' = data.models[, x], data.models[, 280:286])
  mod.a <-
    lm(amr ~ entero + amr.rich + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['enteroBact2', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['enteroBact2', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$amr, model$entero)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(data.models)[x], p.val, delta, delta.low, delta.up)
  result.ent2[[x]] <- result
}
entero2.lm <- do.call(rbind.data.frame, result.ent2)
names(entero2.lm) <-
  c('amr', 'pval', 'delta', 'delta.low', 'delta.up')
entero2.lm$fdr <- p.adjust(entero2.lm$pval, method = 'fdr')
entero2.lm <- entero2.lm[!is.na(entero2.lm$fdr),]
dim(entero2.lm[entero2.lm$fdr<.1,])
