## phages
setwd('/home/scratch/margar')

##
library(tidyverse); library(ggpubr)

##
phages <- read.table('data/phageome.txt', header=T, row.names=1)
phages <- as.data.frame(t(phages))
phages[is.na(phages)] <- 0
phages.rar <- GUniFrac::Rarefy(phages)#, min(colSums(phages, na.rm=T)))
phages.rar <- as.data.frame(phages.rar$otu.tab.rff)
#
t2 <- phages.rar
Phageome_before_filtering <- t2
Phageome_logic <- Phageome_before_filtering  == 0
colsum <- colSums(Phageome_logic)
phage_sig <- names(colsum[colsum<700])
Phageome_after_filtering <- Phageome_before_filtering[, phage_sig]
phages.rar <- Phageome_after_filtering
rm(list=c('t2', 'Phageome_before_filtering', 'Phageome_logic', 'colsum', 'phage_sig', 'Phageome_after_filtering'))
#
phages.rar.clr <- as.data.frame(compositions::clr(phages.rar))

# rename
mtdt <- read.table('data/curatedmetadata/WP2.1_BL_MetadataWithInfl2.txt', header=T)
row.names(phages.rar.clr) <- gsub('X', '', row.names(phages.rar.clr))
mtdt$SampleID == row.names(phages.rar.clr)
phages.rar.clr <- phages.rar.clr[match(mtdt$SampleID, row.names(phages.rar.clr)),]
row.names(phages.rar.clr) <- row.names(mtdt)
saveRDS(phages.rar.clr, 'Phages_CLR.rds')
phages.rar.clr <- readRDS('/home/scratch/margar/data/Phages_CLR.rds')

## control for age, sex and center
mtdt.bsl <- read.table('data/curatedmetadata/WP2.1_BL_MetadataWithInfl2.txt', header=T)
phages.rar.clr <- cbind(phages.rar.clr, 'Center' = mtdt.bsl$CenterID, 'Age' = mtdt.bsl$Age, 'Gender' = mtdt.bsl$Gender)
residuals.df <- list()
for (phage in seq(1:347)){
  model <- lm(phages.rar.clr[,phage] ~ Center + Age + Gender, data = phages.rar.clr)
  residuals.phage <- residuals(model)
  residuals.df[[phage]]<- residuals.phage
}

residuals.df <- do.call(cbind.data.frame, residuals.df)
names(residuals.df) <- names(phages.rar.clr)[1:347]
phages.rar.clr <- residuals.df

## phages and baseline phenotype ----
# mtdt.bsl <- read.table('data/curatedmetadata/WP2.1_BL_MetadataWithInfl2.txt', header=T, row.names=1)
row.names(mtdt.bsl) == row.names(phages.rar.clr)

## things to load before ----
source('scripts/lm_and_correlations_functions.R')
blank_theme <- function() {
  theme(
    # add border 1)
    panel.border = element_rect(colour = "black", fill = NA),
    # color background 2)
    panel.background = element_blank(),
    # modify grid 3)
    panel.grid.major = element_line(colour = "gray85", linetype = 'dashed'),
    # modify text, axis and colour 4) and 5)
    axis.text = element_text(colour = "black", family = "Times New Roman"),
  )
}

# some plot previously ----
# PCA
phages.byc <- vegdist(phages.rar, method='bray')
phages.byc.pcoa <- ape::pcoa(phages.byc)
phages.byc.pcoa.scores <- as.data.frame(phages.byc.pcoa$vectors)
row.names(phages.byc.pcoa.scores) == row.names(mtdt.bsl)
ggplot(phages.byc.pcoa.scores, aes(Axis.1, Axis.2, color=as.factor(mtdt.bsl$CenterID))) +
  geom_point() +
  stat_ellipse()
adonis(phages.rar~as.factor(mtdt.bsl$CenterID))

phages.ait <- coda.base::dist(phages.rar, method='aitchison')
phages.ait.pcoa <- ape::pcoa(phages.ait)
phages.ait.pcoa.scores <- as.data.frame(phages.ait.pcoa$vectors)
row.names(phages.ait.pcoa.scores) == row.names(mtdt.bsl)
ggplot(phages.ait.pcoa.scores, aes(Axis.1, Axis.2, color=as.factor(mtdt.bsl$CenterID))) +
  geom_point()


# richness
row.names(phages.rar) <- gsub('X', '', row.names(phages.rar))
row.names(phages.rar) == mtdt.bsl$SampleID
phages.rar <- phages.rar[match(mtdt.bsl$SampleID, row.names(phages.rar)),]
row.names(phages.rar) <- row.names(mtdt.bsl)
phages.rich <- rowSums(phages.rar!=0)

qmp.bsl <- read.table('data/WP2.1_775_IGR_MGS-QMP_igrHA_Ent_M.csv', header=T, sep='\t', row.names=1)
qmp.bsl <- qmp.bsl[,7:727]
row.names(qmp.bsl) == mtdt.bsl$samplerenamed
qmp.bsl <- data.frame(t(qmp.bsl))
mgs.rich <- colSums(qmp.bsl!=0)

phages.rich <- data.frame(phages.rich, 'GR' = mtdt.bsl$DSGeneAbundRichness, 'tertile' = mtdt.bsl$DSGARLTertileGroup, 'enterotype' = mtdt.bsl$enterotype, 'mgr.rich' = mgs.rich)
phages.rich$tertile <- plyr::revalue(phages.rich$tertile, c('EGC' = 'HGR', 'HGC' = 'MGR', 'LGC' = 'LGR'))
phages.rich$tertile  <- factor(phages.rich$tertile, levels=c('HGR', 'MGR', 'LGR'))
cor.gr <- ggplot(phages.rich, aes(GR, phages.rich)) + 
  geom_point(aes(color = tertile)) +
  geom_smooth(method='lm') +
  stat_cor() +
  scale_color_manual(values=c('blue', 'gray75', 'gold2')) +
  blank_theme() +
  labs(x='Bacterial gene counts', y = 'Number of different phages', color = 'Gene richness\ntertile')

cor.mgs <- ggplot(phages.rich, aes(mgs.rich, phages.rich)) + 
  geom_point() +
  geom_smooth(method='lm') +
  stat_cor() +
  # scale_color_manual(values=c('blue', 'gray75', 'gold2')) +
  blank_theme() +
  labs(x='Numer of different MGS', y = 'Number of different phages', color = 'Gene richness\ntertile')

cowplot::plot_grid(cor.gr, cor.mgs, rel_widths = c(1,0.75))


phages.rich$tertile  <- factor(phages.rich$tertile, levels=c('LGR', 'MGR', 'HGR'))
gr.bxs <- ggplot(phages.rich, aes(tertile, phages.rich, group=interaction(tertile))) +
  geom_violin(alpha = .1, aes(fill=tertile, color=tertile)) +
  geom_quasirandom(alpha = .6, dodge.width = .9, aes(color = tertile)) +
  geom_boxplot(width = .1, alpha = .7, outlier.colour = NA, position = position_dodge(.9), fill='white') +
  stat_compare_means() +
  stat_compare_means(comparisons = list(c('LGR', 'MGR'), c('LGR', 'HGR'), c('MGR', 'HGR'))) + 
  scale_fill_manual(values = c('blue', 'gray75', 'gold2')) +
  scale_color_manual(values = c('blue', 'gray75', 'gold2')) +
  theme_bw() +
  labs(x = NULL, y = 'Number of different phages', colour = 'Gene richness\ntertile', fill = 'Gene richness\ntertile')


ent.bxs <- ggplot(phages.rich, aes(enterotype, phages.rich, group=interaction(enterotype))) +
  geom_violin(alpha = .1, aes(fill=enterotype, color=enterotype)) +
  geom_quasirandom(alpha = .6, dodge.width = .9, aes(color = enterotype)) +
  geom_boxplot(width = .1, alpha = .7, outlier.colour = NA, position = position_dodge(.9), fill='white') +
  stat_compare_means() +
  stat_compare_means(comparisons = list(c('Bact2', 'Bact1'), c('Bact2', 'Rum'), c('Bact2', 'Prev'))) + 
  scale_fill_manual(values = c('darkolivegreen4', 'red', 'cadetblue3', 'darkgoldenrod2', 'blue3')) +
  scale_color_manual(values = c('darkolivegreen4', 'red', 'cadetblue3', 'darkgoldenrod2', 'blue3')) +
  theme_bw() +
  labs(x = NULL, y = 'Number of different phages', colour = "Enterotype\nstratifications", fill = "Enterotype\nstratifications")

cowplot::plot_grid(gr.bxs, ent.bxs)
# phages.bsl <- lm.associations(mtdt.bsl, phages.rar.clr, names(mtdt.bsl)[15:56])
# phages.bsl <- do.call(rbind.data.frame, phages.bsl)

phages.bsl.cor <- cor.associations(mtdt.bsl, phages.rar.clr, names(mtdt.bsl)[15:56])
phages.bsl.cor <- do.call(rbind.data.frame, phages.bsl.cor)

heat.sort <- phages.bsl.cor[phages.bsl.cor$variable == 'BMI',]
heat.sort <- heat.sort[order(heat.sort$rho),]
heat.sort <- factor(heat.sort$metabolite, levels=heat.sort$metabolite)

phages.bsl.cor %>% 
  filter(fdr<.1) %>% 
  mutate('association' = gsub('\\.[0-9]+', '', row.names(.))) %>% 
  mutate('association'= factor(association, levels=unique(association))) %>% 
  mutate('association' = plyr::revalue(association,
                                       c('Height' = 'Height (cm)',
                                                       'Weight.kg' = 'Weight (kg)',
                                                       'Waist.cm' = 'Waist (cm)', 
                                                       'Hip.cm' =  'Hip (cm)',
                                                       'Waist.Hip' = 'Waist:Hip',
                                                       'BMI' = 'Body Mass Index (kg/m2)',
                                                       'Fasting.Glucose' = 'Fasting P. Glucose (mmol/L)',
                                                       'Basal.Glucose' = 'Basal P. Glucose (mmol/L)',
                                                       'Mean.Glucose' = 'Mean OGTT P. Glucose (mmol/L)',
                                                       'glucoseh' = 'OGTT P. Glucose 2-hours (mmol/L)',
                                                       'Fast.Insulin' = 'Fasting P. Insulin (pmol/L)',
                                                       'Basal.Insulin' = 'Basal P. Insulin (pmol/L)',
                                                       'Mean.Insulin' = 'OGTT P. Insulin (pmol/L)',
                                                       'Basal.InsulinSecretion' = 'Basal insulin secretion (pmol min-1m-2)',
                                                       'Total.InsulinSecretion' = 'OGTT insulin secretion (nmol m-2)',
                                                       'Glucose.Sensitivity' = 'Glucose sensitivity (pmol min-1 m-2 (mmol/L)-1)',
                                                       'OGISh'= 'OGTT insulin sensitivity (ml min-1m-2)',
                                                       'Stumvoll' = "Stumvoll's index (umol min-1kg-1)",
                                                       'Matsuda' = "Matsuda's index (umol min-1kg-1)",
                                                       'Clinsb' = 'Basal insulin clearance (L min-1m-2)',
                                                       'Clins' = 'OGTT insulin clearance (L min-1m-2)',
                                                       'Active.GLP1.Concm' = "Fasting P. GLP1 (pg/ml)",
                                                       'Total.GLP1.Concm' = 'Total P. OGTT GLP1 (pg/ml)',
                                                       'hsCRP' = "Fasting P. hsCRP (mg/L)",
                                                       'Fasting.HDL' = "Fasting P. HDL  (mmol/L)",
                                                       'Fasting.LDL' = "Fasting P. LDL  (mmol/L)",
                                                       'Fasting.TG' = "Fasting P. Triglycerides  (mmol/L)",
                                                       'Fasting.ALT' = "Fasting P. ALT (U/L)",
                                                       'Fasting.AST' = "Fasting P. AST (U/L)",
                                                       'Fasting.Chol' = "Fasting P. Cholesterol  (mmol/L)",
                                                       'Liver.Iron' = "Liver iron (umol/g)",
                                                       'Panc.Iron' = "Pancreatic iron (umol/g)",
                                                       'IAAT' = "IAAT (L)",
                                                       'ASAT' = "ASAT (L)",
                                                       'TAAT' = "TAAT (L)",
                                                       'LiverFat..' = "Liver fat (%)",
                                                       'PancreasFat..' = "Pancreatic fat (%)",
                                                       'BP.Sys' = "Mean Systolic B.P. (mmHg)",
                                                       'BP.D.Mean' = "Mean Diastolic B.P. (mmHg)",
                                                       'Value.vm.hpf.mean.BL' = "Physical activity",
                                                       'HDI' = "Healthy Dietary Index (HDI)"))) %>% 
  mutate('rho' = as.numeric(rho)) %>% 
  mutate('metabolite' = factor(metabolite, levels=heat.sort)) %>% 
  ggplot(aes(association, metabolite, fill=rho)) +
  geom_tile() + scale_fill_gradient2(low='blue', mid='white', high='red') +
    blank_theme() + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5)) + 
    labs(x=NULL, y= NULL, fill = "rho", 
       title = 'Phages vs phenotype')

phages.bsl.cor %>% 
  filter(fdr<.1) %>% 
  select(metabolite) %>% 
  unique() %>% 
  write.table(., '/home/margar/download/phagesSIG.txt')

# phages.bsl %>% 
#   filter(fdr<.1) %>% 
#   mutate('association' = gsub('\\.[0-9]+', '', row.names(.))) %>% 
#   mutate('delta' = as.numeric(delta)) %>% 
#   ggplot(aes(association, mgs, fill=delta)) +
#   geom_tile() + scale_fill_gradient2(low='blue', mid='white', high='red') +
#   theme(panel.border = element_rect(fill=NA), panel.background = element_blank(),
#         panel.grid.major = element_line(color='gray85', linetype='dashed'),
#         axis.text.x = element_text(angle=90, vjust=0.5, size=10),
#         axis.title = element_text(size = 10, face = "bold"), axis.text = element_text(face = "bold"),
#         legend.text = element_text(face = "bold"), legend.title = element_text(face = "bold"),
#         legend.background = element_blank(), legend.box.background = element_blank(), legend.key = element_blank()) +
#   labs(x=NULL, y= NULL, fill = "Cliff's delta Effect Size", 
#        title = 'Phages vs phenotype')

## phages and delta values ----
delta.df <- readRDS('data/DeltaofLogValues.rds')
row.names(delta.df) <- delta.df$studyID
phages.delta <- phages.rar.clr[row.names(phages.rar.clr)%in%row.names(delta.df),]
row.names(phages.delta) == row.names(delta.df)
phages.delta.cor <- cor.associations(delta.df, phages.delta, names(delta.df)[6:24])
phages.delta.cor <- do.call(rbind.data.frame, phages.delta.cor)
phages.delta.cor %>% 
  filter(fdr<.1) %>% 
  mutate('association' = gsub('\\.[0-9]+', '', row.names(.))) %>% 
  mutate('association'= factor(association, levels=unique(association))) %>% 
  mutate('association' = plyr::revalue(association,
                                      c('Fast.Insulin' = 'Fasting P. Insulin (pmol/L)',
                                         'Basal.Insulin' = 'Basal P. Insulin ((pmol/L)',
                                         'Total.InsSecr' = 'OGTT insulin secretion (nmol m-2)',
                                         'Stumvoll' = "Stumvoll's index (umol min-1kg-1)"))) %>% 
  mutate('rho' = as.numeric(rho)) %>% 
  ggplot(aes(association, metabolite, fill=rho)) +
  geom_tile() + scale_fill_gradient2(low='blue', mid='white', high='red') +
  blank_theme() + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5)) + 
  labs(x=NULL, y= NULL, fill = "rho", 
       title = 'Phages vs phenotype progression')

## pages and progression AS CONTINUOUS (fast.glucose followup) ----
mtdt.fw <- read.table('data/curatedmetadata/WP2.1_BlFw-641_MetadataImputed_OnlyBIOCHEM_2.txt', header=T)
fast.glucose <- mtdt.fw %>% 
  filter(timepoint=='Follow-up') %>% 
  select(c('StudyID', 'Fast.Glucose', 'Age', 'CenterID', 'Gender'))

fast.glucose$StudyID == row.names(phages.delta)
fast.gluc.phage <- cor.associations(fast.glucose, phages.delta, 'Fast.Glucose') 
fast.gluc.phage <- as.data.frame(fast.gluc.phage); names(fast.gluc.phage) <- c('var', 'phage', 'rho', 'pval', 'fdr')
fast.gluc.phage %>% 
  filter(fdr<.1)

## phages and progression AS DISCRETE ----
progression <- mtdt.fw %>% 
  filter(timepoint=='Follow-up') %>% 
  select(c('StudyID', 'Gly.Cat', 'Age', 'CenterID', 'Gender'))
progression$StudyID == row.names(phages.delta)
progression$ngr <- plyr::revalue(progression$Gly.Cat, c('T2D'='noNGR', 'IGR'='noNGR'))
progression$t2d <- plyr::revalue(progression$Gly.Cat, c('NGR'='noT2D', 'IGR'='noT2D'))
phages.t2d <- lm.associations(progression, phages.delta, 't2d')
phages.t2d <- as.data.frame(phages.t2d); names(phages.t2d) <- c('phage', 'pval', 'delta', 'delta.low', 'delta.up', 'fdr')
phages.t2d %>% 
  filter(fdr<.1)
phages.ngr <- lm.associations(progression, phages.delta, 'ngr')
phages.ngr <- as.data.frame(phages.ngr); names(phages.ngr) <- c('phage', 'pval', 'delta', 'delta.low', 'delta.up', 'fdr')
phages.ngr %>% 
  filter(fdr<.1)

progression.2 <- progression[progression$Gly.Cat!='IGR',]
phages.delta2 <- phages.delta[row.names(phages.delta)%in%progression.2$StudyID,]
row.names(phages.delta2) == progression.2$StudyID
extreme.pr <- lm.associations(progression.2, phages.delta2, 'Gly.Cat')
extreme.pr <- as.data.frame(extreme.pr); names(extreme.pr) <- c('phage', 'pval', 'delta', 'delta.low', 'delta.up', 'fdr')
extreme.pr %>% 
  filter(fdr<.1)

phages.delta$t2d <- progression$t2d
ggplot(phages.delta, aes(t2d, uvig_17786)) + 
  geom_violin() + stat_compare_means()

ggplot(phages.delta, aes(t2d, uvig_92258)) + 
  geom_violin() + stat_compare_means()

## procrustes phages - QMP ----
qmp <- read.table('data/WP2.1_775_IGR_MGS-QMP_igrHA_Ent_M.csv', header=T, row.names=1)
qmp <- qmp[,7:727]
row.names(qmp) == mtdt.bsl$samplerenamed
row.names(qmp) <- row.names(mtdt.bsl)
qmp.byc <- vegdist(qmp, method='bray')
qmp.byc.pcoa <- ape::pcoa(qmp.byc, scale = T, .center=T)

qmp.phages.proc <- procrustes(qmp.byc.pcoa$vectors, phages.byc.pcoa$vectors)
plot(qmp.phages.proc, kind=1)
protest(qmp.byc.pcoa$vectors, phages.byc.pcoa$vectors)


qmp.phages.proc.df <- data.frame(x1 = qmp.phages.proc$Yrot[,1], y1 = qmp.phages.proc$Yrot[,2],
                     x2 = qmp.phages.proc$X[,1], y2 = qmp.phages.proc$X[,2])


ggplot(qmp.phages.proc.df) +
  geom_point(aes(x1, y1), colour = 'darkolivegreen2', shape = 19, size = 5) +
  geom_point(aes(x2, y2), colour = 'cornflowerblue', shape = 17, size = 5) +
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), arrow=arrow(length=unit(.3, 'cm')), alpha=.15) +
  geom_hline(yintercept = 0, linetype = 'dashed') + geom_vline(xintercept = 0, linetype = 'dashed') +
  blank_theme() +
  labs(x = "Dimension 1", y = "Dimension 2") +
  annotate('text', x=0.4 , y=-0.2, label = 'Correlation 0.582 \n Monte-Carlo p = 0.001', fontface =2)

phages.pca <- prcomp(phages.byc, center = T, scale. = T)
qmp.pca <- prcomp(qmp.byc, center = T, scale.=T)
qmp.phages.proc2 <- procrustes(qmp.pca, phages.pca)
plot(qmp.phages.proc, kind=1)
protest(qmp.pca, phages.pca)


qmp.phages.proc.df <- data.frame(x1 = qmp.phages.proc$Yrot[,1], y1 = qmp.phages.proc$Yrot[,2],
                                 x2 = qmp.phages.proc$X[,1], y2 = qmp.phages.proc$X[,2])


ggplot(qmp.phages.proc.df) +
  geom_point(aes(x1, y1), colour = 'darkolivegreen2', shape = 19, size = 5) +
  geom_point(aes(x2, y2), colour = 'cornflowerblue', shape = 17, size = 5) +
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), arrow=arrow(length=unit(.3, 'cm')), alpha=.15) +
  geom_hline(yintercept = 0, linetype = 'dashed') + geom_vline(xintercept = 0, linetype = 'dashed') +
  theme(plot.subtitle = element_text(vjust = 1), plot.caption = element_text(vjust = 1),
        axis.title = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 14, face = "bold", colour = "black"),
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(fill=NA, color = 'black', size = .5)) +
  labs(x = "Dimension 1", y = "Dimension 2") +
  annotate('text', x=0.4 , y=-0.2, label = 'Correlation 0.582 \n Monte-Carlo p = 0.001', fontface =2)

## phages josef ----
phages.josef <- phages.delta.cor %>% 
  filter(fdr<.1) %>% 
  select('metabolite') 



phages.josef <- unique(c(phages.josef$metabolite, 'uvig_17786', 'uvig_92258'))
phages.rar.clr2 <- readRDS('/home/scratch/margar/data/Phages_CLR.rds')
phages.rar.clr2 <- phages.rar.clr2[,names(phages.rar.clr2)%in%phages.josef]
## control for age, sex and center
phages.rar.clr2 <- cbind(phages.rar.clr2, 'Center' = mtdt.bsl$CenterID)

residuals.df <- list()
for (phage in seq(1:9)){
  model <- lm(phages.rar.clr2[,phage] ~ Center, data = phages.rar.clr2)
  residuals.phage <- residuals(model)
  residuals.df[[phage]]<- residuals.phage
}

residuals.df <- do.call(cbind.data.frame, residuals.df)
names(residuals.df) <- names(phages.rar.clr2)[1:9]
phages.rar.clr2 <- residuals.df
saveRDS(phages.rar.clr2, 'data/PhagesForJosef.rds')


## load phylogenetic tree from clustal omega alignment
library(ggtree)
phages.tree <- ape::read.tree('/home/margar/upload/phagestree.txt')
phages.tree.df <- fortify(phages.tree)
dd <- subset(phages.tree.df, isTip)
heatmap.order <- factor(dd$label, levels=dd$label[order(dd$y, decreasing = T)])

phages.bsl.cor %>% 
  filter(fdr<.1) %>% 
  mutate('association' = gsub('\\.[0-9]+', '', row.names(.))) %>% 
  mutate('association'= factor(association, levels=unique(association))) %>% 
  mutate('association' = plyr::revalue(association,
                                       c('Height' = 'Height (cm)',
                                         'Weight.kg' = 'Weight (kg)',
                                         'Waist.cm' = 'Waist (cm)', 
                                         'Hip.cm' =  'Hip (cm)',
                                         'Waist.Hip' = 'Waist:Hip',
                                         'BMI' = 'Body Mass Index (kg/m2)',
                                         'Fasting.Glucose' = 'Fasting P. Glucose (mmol/L)',
                                         'Basal.Glucose' = 'Basal P. Glucose (mmol/L)',
                                         'Mean.Glucose' = 'Mean OGTT P. Glucose (mmol/L)',
                                         'glucoseh' = 'OGTT P. Glucose 2-hours (mmol/L)',
                                         'Fast.Insulin' = 'Fasting P. Insulin (pmol/L)',
                                         'Basal.Insulin' = 'Basal P. Insulin (pmol/L)',
                                         'Mean.Insulin' = 'OGTT P. Insulin (pmol/L)',
                                         'Basal.InsulinSecretion' = 'Basal insulin secretion (pmol min-1m-2)',
                                         'Total.InsulinSecretion' = 'OGTT insulin secretion (nmol m-2)',
                                         'Glucose.Sensitivity' = 'Glucose sensitivity (pmol min-1 m-2 (mmol/L)-1)',
                                         'OGISh'= 'OGTT insulin sensitivity (ml min-1m-2)',
                                         'Stumvoll' = "Stumvoll's index (umol min-1kg-1)",
                                         'Matsuda' = "Matsuda's index (umol min-1kg-1)",
                                         'Clinsb' = 'Basal insulin clearance (L min-1m-2)',
                                         'Clins' = 'OGTT insulin clearance (L min-1m-2)',
                                         'Active.GLP1.Concm' = "Fasting P. GLP1 (pg/ml)",
                                         'Total.GLP1.Concm' = 'Total P. OGTT GLP1 (pg/ml)',
                                         'hsCRP' = "Fasting P. hsCRP (mg/L)",
                                         'Fasting.HDL' = "Fasting P. HDL  (mmol/L)",
                                         'Fasting.LDL' = "Fasting P. LDL  (mmol/L)",
                                         'Fasting.TG' = "Fasting P. Triglycerides  (mmol/L)",
                                         'Fasting.ALT' = "Fasting P. ALT (U/L)",
                                         'Fasting.AST' = "Fasting P. AST (U/L)",
                                         'Fasting.Chol' = "Fasting P. Cholesterol  (mmol/L)",
                                         'Liver.Iron' = "Liver iron (umol/g)",
                                         'Panc.Iron' = "Pancreatic iron (umol/g)",
                                         'IAAT' = "IAAT (L)",
                                         'ASAT' = "ASAT (L)",
                                         'TAAT' = "TAAT (L)",
                                         'LiverFat..' = "Liver fat (%)",
                                         'PancreasFat..' = "Pancreatic fat (%)",
                                         'BP.Sys' = "Mean Systolic B.P. (mmHg)",
                                         'BP.D.Mean' = "Mean Diastolic B.P. (mmHg)",
                                         'Value.vm.hpf.mean.BL' = "Physical activity",
                                         'HDI' = "Healthy Dietary Index (HDI)"))) %>% 
  mutate('rho' = as.numeric(rho)) %>% 
  mutate('metabolite' = factor(metabolite, levels=heatmap.order)) %>% 
  ggplot(aes(association, metabolite, fill=rho)) +
  geom_tile() + scale_fill_gradient2(low='blue', mid='white', high='red') +
  blank_theme() + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5)) + 
  labs(x=NULL, y= NULL, fill = "rho", 
       title = 'Phages vs phenotype')
pdf('/home/margar/download/PHAGES_tree.pdf', height = 8, width = 8)
ggtree(phages.tree, branch.length = 'none') + geom_tiplab()
dev.off()
