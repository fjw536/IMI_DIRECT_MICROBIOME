######################################################################################################
############################# Metabolomics data preparation ##########################################
######################################################################################################

## directory
setwd('/home/scratch/margar')

## libraries
library(factoextra); library(mdatools); library(tidyverse); library(purrr); library(vroom); library(fs)
library(ggplot2); library(ggpubr); library(ggbeeswarm); library(ROCR); library(corrplot)

## upload (and format) data
mets <- read.table('data/metabolomics/DTU_imput_wp21_kk_22sh.txt', header=T, row.names=1)
names(mets) <- gsub('X', '', names(mets))
annotations <- read.table("/home/scratch/margar/data/metabolomics/Chemical_Annotation_Metabolone_DTU_sh.txt", row.names = 1, header=T)
names(mets) == row.names(annotations)

mtdt.bsl <- read.table('data/curatedmetadata/WP2.1_BL_MetadataWithInfl2.txt', header=T, row.names =1, sep ='\t')
mtdt.all <- read.table('data/curatedmetadata/WP2.1_BlFw-641_MetadataImputed_OnlyBIOCHEM_2.txt', header=T, sep = '\t')
mtdt.all[mtdt.all$timepoint=='Baseline',]$StudyID == mtdt.all[mtdt.all$timepoint!='Baseline',]$StudyID 
mtdt.all$progression <- rep(mtdt.all[mtdt.all$timepoint!='Baseline',]$Gly.Cat, 2)
mtdt.all <- mtdt.all[mtdt.all$timepoint=='Baseline',]

mets[] <- as.data.frame(sapply(mets, log))
mets[] <- as.data.frame(scale(mets))

row.names(mets) == row.names(mtdt.bsl)
reduced.mets <- mets[row.names(mets)%in%mtdt.all$StudyID,]

## PCA metabolites vs metadata strats
mets.pca <- prcomp(mets, center=T, scale=T)
fviz_screeplot(mets.pca)

fviz_pca_biplot(mets.pca)
center.plot <- fviz_pca_ind(mets.pca, geom='point', col.ind = as.factor(mtdt.bsl$CenterID), addEllipses = T) + ggtitle(label='Center PCA')
gender.plot <- fviz_pca_ind(mets.pca, geom='point', col.ind = as.factor(mtdt.bsl$Sex), addEllipses = T) + ggtitle(label='Gender PCA')
age.plot <- fviz_pca_ind(mets.pca, geom='point', col.ind = mtdt.bsl$Age, addEllipses = F) + ggtitle(label='Age PCA')

prediab.plot <- fviz_pca_ind(mets.pca, geom='point', col.ind = as.factor(mtdt.bsl$Gly.Cat2), addEllipses = T) + ggtitle(label='Prediabetic type PCA')
gr.plot <- fviz_pca_ind(mets.pca, geom='point', col.ind = as.factor(mtdt.bsl$DSGARLTertileGroup), addEllipses = T) + ggtitle(label='Gene richness tertiles PCA')
entero.plot <- fviz_pca_ind(mets.pca, geom='point', col.ind = as.factor(mtdt.bsl$enterotype), addEllipses = T) + ggtitle(label='Enterotypes PCA')

cowplot::plot_grid(center.plot, gender.plot, age.plot, prediab.plot, gr.plot, entero.plot, ncol = 3)
loadings <- fviz_pca_var(mets.pca, geom='point', col.var=annotations$SUPER_PATHWAY, addEllipses=F, size.var=10, alpha.var=.6) + 
  scale_shape_manual(values=c(15, 16, 17, 18, 19, 20, 7, 8, 9, 14))


reduced.mets.pca <- prcomp(reduced.mets, center = T, scale=T)
fviz_pca_ind(reduced.mets.pca, geom='point', col.ind = as.factor(mtdt.all$progression), addEllipses = T)

## PLS-DA on progression
reduced.mets$progression <- mtdt.all$progression
set.seed(99)
sample.split <- caTools::sample.split(row.names(reduced.mets), SplitRatio = .8)
train <- subset(reduced.mets, sample.split==T)
test <- subset(reduced.mets, sample.split==F)

progression.plsda <- plsda(train[,1:973], train[,974], 7, scale = TRUE, cv = 1)
summary(progression.plsda)
getConfusionMatrix(progression.plsda$calres) ## it does not seem to be able to identify NGR and T2D...

plotPredictions(progression.plsda)
plotMisclassified(progression.plsda)
plotSensitivity(progression.plsda)
plotSpecificity(progression.plsda)

plotRegcoeffs(progression.plsda, ncomp=1, show.ci=T)

res <- predict(progression.plsda, test[,1:973], test[,974])
summary(res)

## PLS-DA on prediabetic types ----
mets$prediabetics <- mtdt.bsl$Gly.Cat2
set.seed(97)
sample.split <- caTools::sample.split(row.names(mets), SplitRatio = .8)
train <- subset(mets, sample.split==T)
test <- subset(mets, sample.split==F)

prediabet.plsda <- plsda(train[,1:973], train[,974], 7, scale = TRUE, cv = 1)
summary(prediabet.plsda)
getConfusionMatrix(prediabet.plsda$calres) ## it does not seem to be able to identify NGR and T2D...

plotPredictions(prediabet.plsda)
plotMisclassified(prediabet.plsda)
plotSensitivity(prediabet.plsda)
plotSpecificity(prediabet.plsda)

plotRegcoeffs(prediabet.plsda, ncomp=1, show.ci=T)

res <- predict(prediabet.plsda, test[,1:973], test[,974])
summary(res)

## correction for age, sex, center and EGFr ---- 
mets <- cbind(mets, 'Age' = mtdt.bsl$Age, 'Gender' = mtdt.bsl$Gender, 'Center' = mtdt.bsl$CenterID, 'egfr' = mtdt.bsl$egfr)
residuals.df <- list()
for (met in seq(1:973)){
  model <- lm(mets[,met] ~ Age + Center + Gender + egfr, data = mets)
  residuals.met <- residuals(model)
  residuals.df[[met]]<- residuals.met
}

residuals.df <- do.call(cbind.data.frame, residuals.df)
names(residuals.df) <- names(mets)[1:973]
#saveRDS(residuals.df, 'data/metabolomics/MetabolitesResiduals.rds')
## individual metabolite differences to progression ----
red.mets.corr <- residuals.df[row.names(residuals.df)%in%row.names(reduced.mets),]
red.mets.corr$progression <- reduced.mets$progression

red.mets.corr %>% 
  group_by(progression) %>% 
  summarise(across(where(is.numeric), list(mean=mean)))  %>% 
  pivot_longer(ends_with(c('mean', 'sd')),
               names_to = c('name', '.value'),
               names_sep='_') %>% 
  pivot_wider()

anova.results <- sapply(red.mets.corr[,1:973], function(a){
  aov.model <- lm(a~progression, data=red.mets.corr)
  p.val <- anova(aov.model)$`Pr(>F)`[1]
  return(p.val)
})
anova.fdr <- p.adjust(anova.results, method ='fdr')

pairwise.results <- lapply(red.mets.corr[,1:973], function(a){
  aov.model <- aov(a~progression, data=red.mets.corr)
  pairwise.p <- TukeyHSD(aov.model)$progression[,4] # order NGR-IGR, T2D-IGR, T2D-NGR
  return(pairwise.p)
})
pairwise.results.df <- do.call(rbind.data.frame, pairwise.results)
names(pairwise.results.df) <- c('NGRvsIGR', 'T2DvsIGR', 'T2DvsNGR')
row.names(pairwise.results.df) <- names(pairwise.results)
fc.df <- red.mets.corr %>% 
  group_by(progression) %>% 
  summarise(across(where(is.numeric), list(mean=mean))) %>% 
  pivot_longer(ends_with(c('mean'))) %>% 
  pivot_wider(names_from = progression, values_from=value) 

red.mets.corr$progression <- factor(red.mets.corr$progression, levels = c('NGR', 'IGR', 'T2D'))
delta.ngrt2d <- lapply(red.mets.corr[red.mets.corr$progression!='IGR',1:973], function(h){
  deltas <- effsize::cliff.delta(h, droplevels(red.mets.corr[red.mets.corr$progression!='IGR','progression']))
  delta.res <- data.frame('delta.ngrt2d' = deltas$estimate, 
                          'down.ngrt2d' = deltas$conf.int[[1]],
                          'up.ngrt2d' = deltas$conf.int[[2]])
  return(delta.res)
})
delta.ngrt2d <- do.call(rbind.data.frame, delta.ngrt2d)

delta.ngrigr <- lapply(red.mets.corr[red.mets.corr$progression!='T2D',1:973], function(h){
  deltas <- effsize::cliff.delta(h, droplevels(red.mets.corr[red.mets.corr$progression!='T2D','progression']))
  delta.res <- data.frame('delta.ngrigr' = deltas$estimate, 
                          'down.ngrigr' = deltas$conf.int[[1]],
                          'up.ngrigr' = deltas$conf.int[[2]])
  return(delta.res)
})
delta.ngrigr <- do.call(rbind.data.frame, delta.ngrigr)

delta.igrt2d <- lapply(red.mets.corr[red.mets.corr$progression!='NGR',1:973], function(h){
  deltas <- effsize::cliff.delta(h, droplevels(red.mets.corr[red.mets.corr$progression!='NGR','progression']))
  delta.res <- data.frame('delta.igrt2d' = deltas$estimate, 
                          'down.igrt2d' = deltas$conf.int[[1]],
                          'up.igrt2d' = deltas$conf.int[[2]])
  return(delta.res)
})
delta.igrt2d <- do.call(rbind.data.frame, delta.igrt2d)


export.sig <- data.frame(annotations[,4:7], annotations[,14:16], fc.df, delta.ngrt2d, delta.ngrigr, delta.igrt2d, 'anova.fdr' = anova.fdr, pairwise.results.df)
write.table(export.sig, 'data/metabolomics/residualssignificanceprogression.txt', sep = '\t')

# select significant metabolites from ANOVA test
keep.mets <- names(anova.fdr[anova.fdr<.1]) 
red.mets.corr.plot <- red.mets.corr[,names(red.mets.corr)%in%keep.mets]
red.mets.corr.plot$progression <- red.mets.corr$progression

annotations.red <- annotations[row.names(annotations)%in%keep.mets,]
names(red.mets.corr.plot) <- c(make.names(annotations.red$SHORT_NAME), 'progression')

red.mets.corr.plot$progression <- factor(red.mets.corr.plot$progression, levels = c('NGR', 'IGR', 'T2D'))
ggplot(red.mets.corr.plot, aes(progression, AMP)) + 
  geom_boxplot() +
  geom_beeswarm() +
  stat_compare_means(method='anova') +
  stat_compare_means(method='t.test', comparisons = list(c('NGR', 'IGR'),
                                                         c('NGR', 'T2D'),
                                                         c('IGR', 'T2D')))


pdf('ANOVA_univariant_sig_mets.pdf')
for (met in names(red.mets.corr.plot)[1:103]){
  print(ggplot(red.mets.corr.plot, aes_string('progression', met)) + 
          geom_boxplot() +
          geom_beeswarm() +
          stat_compare_means(method='anova') +
          stat_compare_means(method='t.test', comparisons = list(c('NGR', 'IGR'),
                                                                 c('NGR', 'T2D'),
                                                                 c('IGR', 'T2D')))) +
    ggtitle(label=met)
}
dev.off()

red.mets.corr.plot.m <- reshape2::melt(red.mets.corr.plot)
red.mets.corr.plot.m$family <- rep(annotations.red$SUPER_PATHWAY, each = 641)
red.mets.corr.plot.m$pathway <- rep(annotations.red$SUB_PATHWAY, each = 641)
red.mets.corr.plot.m$name <- rep(annotations.red$SHORT_NAME, each = 641)
red.mets.corr.plot.m$name <- gsub('2-hydroxybutyrate/2-hydroxyisobutyrate',  '2-hydroxybutyrate\n/2-hydroxyisobutyrate', red.mets.corr.plot.m$name)
red.mets.corr.plot.m$name <- gsub('palmitoyl-sphingosine-phosphoethanolamine',  'palmitoyl-sphingosine-\nphosphoethanolamine', red.mets.corr.plot.m$name)
red.mets.corr.plot.m$name <- gsub('1-\\(1-enyl-palmitoyl)-2-',  '1-\\(1-enyl-palmitoyl\\)-2-\n', red.mets.corr.plot.m$name)
red.mets.corr.plot.m$name <- gsub('1-oleoyl-2-docosahexaenoyl-',  '1-oleoyl-2-docosahexaenoyl-\n', red.mets.corr.plot.m$name)
red.mets.corr.plot.m$name <- gsub('1-stearoyl-2-docosahexaenoyl-',  '1-stearoyl-2-docosahexaenoyl-\n', red.mets.corr.plot.m$name)
red.mets.corr.plot.m$name <- gsub('glycosyl-N-palmitoyl-',  'glycosyl-N-palmitoyl-\n', red.mets.corr.plot.m$name)
red.mets.corr.plot.m$name <- gsub('lactosyl-N-palmitoyl-',  'lactosyl-N-palmitoyl-\n', red.mets.corr.plot.m$name)
red.mets.corr.plot.m$name <- gsub('linoleoyl-N-arachidnoyl-',  'linoleoyl-N-arachidnoyl-\n', red.mets.corr.plot.m$name)
red.mets.corr.plot.m$name <- gsub('sphingomyelin',  'sphingomyelin\n', red.mets.corr.plot.m$name)
red.mets.corr.plot.m$name <- gsub('glycerol',  'glycerol\n', red.mets.corr.plot.m$name)
red.mets.corr.plot.m$name <- gsub('glycosyl ceramide',  'glycosyl ceramide\n', red.mets.corr.plot.m$name)


red.mets.corr.plot.m %>% 
  filter(family!='Lipid'&family!='Amino Acid') %>% 
  ggplot(aes(name, value, fill=progression)) + 
  geom_boxplot(position='dodge', outlier.size=.4) + facet_wrap(~family, scales = 'free', ncol=4) +
  scale_x_discrete(labels = scales::label_wrap(10)) +
  scale_fill_manual(values=c('darkolivegreen4', 'red', 'darkorchid')) +
  theme(panel.grid.major = element_line(colour = "gray85", linetype = "dashed"), 
        axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90), 
        legend.text = element_text(face = "bold"), legend.title = element_text(face = "bold"), panel.background = element_rect(fill = NA), 
        legend.key = element_rect(fill = NA), legend.background = element_rect(fill = NA), panel.border = element_rect(fill=NA),
        strip.background = element_rect(fill=NA), strip.text.x = element_text(face='bold')) +
  labs(x = NULL, y = "Residuals", fill = "Progression")

bxs.lipids <- red.mets.corr.plot.m %>% 
  filter(family=='Lipid') %>% 
  ggplot(aes(name, value, fill=progression)) + 
  geom_boxplot(position='dodge', outlier.size = .5) + facet_wrap(~pathway, scales = 'free_x', nrow = 1) + 
  scale_fill_manual(values=c('darkolivegreen4', 'red', 'darkorchid')) +
  theme(panel.grid.major = element_line(colour = "gray85", linetype = "dashed"), 
        axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90), 
        legend.text = element_text(face = "bold"), legend.title = element_text(face = "bold"), panel.background = element_rect(fill = NA), 
        legend.key = element_rect(fill = NA), legend.background = element_rect(fill = NA), panel.border = element_rect(fill=NA),
        strip.background = element_rect(fill=NA), strip.text.x = element_text(face='bold'), legend.position = 'none') +
  labs(x = NULL, y = "Residuals", fill = "Progression", title = 'Lipids') + ylim(-6, 5)
gp.lipids <- ggplotGrob(bxs.lipids)
# gtable::gtable_show_layout(gp.lipids)
facet.columns <- gp.lipids$layout$l[grepl("panel", gp.lipids$layout$name)]
x.var <- sapply(ggplot_build(bxs.lipids)$layout$panel_scales_x,
                function(l) length(l$range$range))
gp.lipids$widths[facet.columns] <- gp.lipids$widths[facet.columns] * x.var
grid::grid.draw(gp.lipids)

bxs.amino <- red.mets.corr.plot.m %>% 
  filter(family=='Amino Acid') %>% 
  ggplot(aes(name, value, fill=progression)) + 
  geom_boxplot(position='dodge', outlier.size = .5) + facet_wrap(~pathway, nrow = 1, scales = 'free_x') +
  scale_fill_manual(values=c('darkolivegreen4', 'red', 'darkorchid')) +
  theme(panel.grid.major = element_line(colour = "gray85", linetype = "dashed"), 
        axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90), 
        legend.text = element_text(face = "bold"), legend.title = element_text(face = "bold"), panel.background = element_rect(fill = NA), 
        legend.key = element_rect(fill = NA), legend.background = element_rect(fill = NA), panel.border = element_rect(fill=NA),
        strip.background = element_rect(fill=NA), strip.text.x = element_text(face='bold', margin=margin(b=5)), legend.position = 'none') +
  labs(x = NULL, y = "Residuals", fill = "Progression", title = 'Amino Acids')
gp.amino<- ggplotGrob(bxs.amino)
# gtable::gtable_show_layout(gp.amino)
facet.columns <- gp.amino$layout$l[grepl("panel", gp.amino$layout$name)]
x.var <- sapply(ggplot_build(bxs.amino)$layout$panel_scales_x,
                function(l) length(l$range$range))
gp.amino$widths[facet.columns] <- gp.amino$widths[facet.columns] * x.var

grid::grid.draw(gp.amino)

pdf('boxplotsaminolipids.pdf', width = 18, height = 13)
gridExtra::grid.arrange(gp.lipids, gp.amino)
dev.off()


final.figure.plot <- red.mets.corr.plot.m %>% 
  group_by(progression, variable) %>% 
  summarise(across(where(is.numeric), median))

final.figure.plot <- as.data.frame(final.figure.plot)
final.figure.plot$family <- annotations.red$SUPER_PATHWAY
final.figure.plot.a <- final.figure.plot[final.figure.plot$family=='Amino Acid'|final.figure.plot$family=='Lipid',]
final.figure.plot.b <- final.figure.plot[final.figure.plot$family!='Amino Acid'&final.figure.plot$family!='Lipid',]

ggplot(final.figure.plot.a, aes(value, variable, fill=progression)) + geom_col(position='dodge') +
  facet_wrap(~family, scales = 'free') + 
  scale_fill_manual(values=c('darkolivegreen4', 'red', 'darkorchid'))

ggplot(final.figure.plot.b, aes(value, variable, fill=progression)) + geom_col(position='dodge') +
  facet_wrap(~family, scales = 'free') + 
  scale_fill_manual(values=c('darkolivegreen4', 'red', 'darkorchid'))

# heatmap
annotation.row <- data.frame('progression' = red.mets.corr.plot$progression); row.names(annotation.row) <- row.names(red.mets.corr.plot)
annotation.col <- data.frame('family' = annotations.red$SUPER_PATHWAY); row.names(annotation.col) <- names(red.mets.corr.plot)[1:103]
pheatmap(red.mets.corr.plot[1:103], scale = 'column', annotation_row = annotation.row, annotation_col = annotation.col)

row.color <- plyr::revalue(annotation.row$progression, c('NGR' = 'darkolivegreen4', 'IGR' = 'red', 'T2D' = 'darkorchid'))
made4::heatplot(red.mets.corr.plot[,1:103], RowSideColors = as.character(row.color))

## PLS-DA on progression
set.seed(101)
sample.split <- caTools::sample.split(row.names(red.mets.corr.plot), SplitRatio = .8)
train <- subset(red.mets.corr.plot, sample.split==T)
test <- subset(red.mets.corr.plot, sample.split==F)

progression.plsda.red <- plsda(train[,1:103], train[,104], 7, scale = TRUE, cv = 1)
vip <- vipscores(progression.plsda.red, ncomp=1)
summary(progression.plsda.red)
getConfusionMatrix(progression.plsda.red$calres) ## it does not seem to be able to identify NGR and T2D...

progression.plsda.red2 <- plsda(train[,1:103], train[,104], 7, scale = TRUE, cv = 1, exclcols = (vip[,3]<0.75)) #filter for T2D-VIP inferior a 0.3
summary(progression.plsda.red2)
getConfusionMatrix(progression.plsda.red2$calres) ## it does not seem to be able to identify NGR and T2D...

lm.model.data <- red.mets.corr.plot 
lm.model.data$progression <- plyr::revalue(lm.model.data$progression, c('IGR' = 'noT2D', 'NGR' = 'noT2D'))
lm.model.data$progression <- as.numeric(lm.model.data$progression)

set.seed(121)
sample.split <- caTools::sample.split(row.names(lm.model.data), SplitRatio = .8)
train <- subset(lm.model.data, sample.split==T)
test <- subset(lm.model.data, sample.split==F)

full.model <- lm(progression~., data=train)
nothing <- lm(progression~1, data=train)
mixed <- step(nothing, direction = 'both',  scope = list(lower=formula(nothing), upper=formula(full.model)), trace = 0)
summary(mixed)
p <- predict(mixed, newdata = test, type='response')
pr <- prediction(p, test$progression)
prf <- performance(pr, measure='tpr', x.measure = 'fpr')
auc.perf <- performance(pr, measure = 'auc')
auc <- unlist(slot(auc.perf,"y.values"))
aucprint <- paste(c("AUC ="),round(auc,digits=3),sep="")
print(aucprint)
plot(prf,col="darkgreen")
abline(a=0,b=1,col="red")

## BASELINE - association to phenotypes ----
mets.residuals <- readRDS('data/metabolomics/MetabolitesResiduals.rds')
row.names(mets.residuals) == row.names(mtdt.bsl)
red.mtdt <- data.frame('Age' = mtdt.bsl$Age, 'GR' = mtdt.bsl$DSGeneAbundRichness,
                       'Shannon' = mtdt.bsl$DSMGSAbundShannon, mtdt.bsl[,15:71])

mets.mtdt <- cbind(mets.residuals, red.mtdt)
cor.mets.mtdt <- cor(mets.mtdt)
cor.mets.mtdt.p <- cor.mtest(mets.mtdt)
pdf('mtdtbslmetabolites.pdf', width = 18, height = 14)
corrplot(cor.mets.mtdt[974:1033, 1:973], p.mat = cor.mets.mtdt.p$p[974:1033, 1:973], insig = 'blank', method = 'square', tl.cex=0.25)
dev.off()

preProcValues.bsl <- preProcess(mtdt.bsl %>% 
                                  select(Age, Height.cm, Weight.kg, Waist.cm, Hip.cm, Waist.Hip, BMI,
                                         Fasting.Glucose, Basal.Glucose, Mean.Glucose, glucose.2h,
                                         Fasting.Insulin, Basal.Insulin, Mean.Insulin, Basal.InsulinSecretion,
                                         Total.InsulinSecretion, Glucose.Sensitivity, OGIS.2h, Stumvoll, 
                                         Matsuda, Clinsb, Clins, Active.GLP1.Conc.0m, Total.GLP1.Conc.0m, Total.GLP1.Conc.60m,
                                         hsCRP, Fasting.HDL, Fasting.LDL, Fasting.TG, Fasting.ALT, Fasting.AST, Fasting.Chol,
                                         Liver.Iron, Panc.Iron, IAAT, ASAT, TAAT, LiverFat.., PancreasFat.., BP.S.Mean, BP.D.Mean,
                                         Value.vm.hpf.mean.BL, HDI, egfr),
                                method = c('knnImpute'),
                                k = 25,
                                knnSummary = mean)
mtdt.bsl <- predict(preProcValues.bsl, mtdt.bsl, na.action=na.pass)

procNames.bsl <- data.frame(col = names(preProcValues.bsl$mean), mean = preProcValues.bsl$mean, sd = preProcValues.bsl$std)
for(i in procNames.bsl$col){
  mtdt.bsl[i] <- mtdt.bsl[i]*preProcValues.bsl$std[i]+preProcValues.bsl$mean[i] 
}

mtdt.bsl$Smoking.Status.BL <- plyr::revalue(as.factor(mtdt.bsl$Smoking.Status.BL), c('ex-smoker'='non-smoker', 'never'='non-smoker', 'current-smoker'='smoker'))
red.mtdt <- data.frame('Age' = mtdt.bsl$Age, 'CenterID' = mtdt.bsl$CenterID, 
                       'Gender' = mtdt.bsl$Gender, 'GR' = mtdt.bsl$DSGeneAbundRichness,
                       'Shannon' = mtdt.bsl$DSMGSAbundShannon, mtdt.bsl[,15:71])
## correlation analyses ----
cor.results <- cor.associations(red.mtdt, mets.residuals, names(red.mtdt)[4:47])

results.dfs.sig <- lapply(cor.results, function(x) {
  x <- x[x$fdr < .1, ]
  return(data.frame(
    'variable' = x$variable,
    'metabolite' = x$metabolite,
    'rho' = x$rho,
    'fdr' = x$fdr
  ))
})

mets.associations <-  do.call(rbind.data.frame, results.dfs.sig)
mets.associations$type <- sapply(mets.associations$metabolite, function(a){
  annotations[grep(paste0('^',a,'$'), row.names(annotations)),]$SUPER_PATHWAY
  }
)
mets.associations$path <- sapply(mets.associations$metabolite, function(a){
  annotations[grep(paste0('^',a,'$'), row.names(annotations)),]$SUB_PATHWAY
  }
)
mets.associations$name <- sapply(mets.associations$metabolite, function(a){
  annotations[grep(paste0('^',a,'$'), row.names(annotations)),]$SHORT_NAME
  }
)

mets.associations$rho <- as.numeric(mets.associations$rho)

direction <- c()
for (i in 1:nrow(mets.associations)){
  if(mets.associations[i,]$rho > 0) {
    direction[[i]] <- c("Positive")
  } else {
    if(mets.associations[i,]$rho  == 0) {
      direction[[i]] <- c("Zero")
    } else {
      direction[[i]] <- c("Negative")
    }
  }}
mets.associations$direction <- unlist(direction)
mets.associations$direction <- factor(mets.associations$direction, levels=c('Positive', 'Negative'))

mets.associations2 <- mets.associations[abs(mets.associations$rho)>0.3,]
for (i in unique(mets.associations2$type)){
  df.plot <- mets.associations2[mets.associations2$type == i,]
  print(ggplot(df.plot, aes(variable, name, fill=rho, shape=direction)) +
    geom_point(size=3.9)  +
    scale_shape_manual(values = c(24, 25)) +
    scale_fill_gradient2(low='blue', high = 'red') +
    theme(panel.border = element_rect(fill=NA), panel.background = element_blank(),
          panel.grid.major = element_line(color='gray85', linetype='dashed'),
          axis.text.x = element_text(angle=90, vjust=0.5, face='bold', color='black', size=12),
          axis.text.y = element_text(face='bold.italic'),
          axis.title = element_text(size = 13, face = "bold"), axis.text = element_text(face = "bold"),
          legend.text = element_text(face = "bold"), legend.title = element_text(face = "bold"),
          legend.background = element_blank(), legend.box.background = element_blank(), legend.key = element_blank()) +
    labs(x = NULL, y = NULL, fill = "Cliff's delta", shape = "Direction", title = paste0('Metabolites associations - ', i)))
}

mets.associations2 <- arrange(mets.associations2, type, path)
mets.associations2$variable <- factor(mets.associations2$variable, levels= names(red.mtdt)[4:47])
order.rows <- cbind('name'=mets.associations2$name, 'type' = mets.associations2$type,'path' = mets.associations2$path)
order.rows <- order.rows[!duplicated(order.rows),]
order.rows <- as.data.frame(order.rows)
order.rows <- arrange(order.rows, type, path) 
mets.associations2$name <- factor(mets.associations2$name, levels= order.rows$name)

pdf('baselinephenotype_metabolites.pdf', height = 25, width = 10)
ggplot(mets.associations2, aes(variable, name, fill = rho)) + 
  geom_tile(colour="white",size=0.2) + scale_fill_gradient2(low='blue', high = 'darkred') +
  theme(axis.text.x=element_text(size=10, vjust = 0.5, angle = 90, hjust=1, face = 'bold'), 
        axis.text.y=element_text(vjust=0.2, face = 'bold.italic'),
        plot.background=element_blank(),panel.border=element_rect(fill=NA),
        plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),
        plot.title=element_text(hjust=0,size=14,face="bold"),
        panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank")) +
  labs(x=NULL, y=NULL)
dev.off()
write.table(order.rows, 'baselinephenotype_metabolites_order.txt', sep = '\t')

head(mets.associations2)
mets.associations %>% 
  group_by(path, variable) %>% 
  summarise_at('rho', median) %>% 
  # filter(abs(rho)>0.1) %>%
  mutate(path=factor(path, levels=order.rows[!duplicated(order.rows$path),]$path)) %>% 
  mutate(variable = factor(variable, levels = names(red.mtdt)[4:47])) %>% 
  ggplot(aes(variable, path, fill = rho)) + 
  geom_tile(colour="white",size=0.2) + scale_fill_gradient2(low='blue', high = 'darkred') +
  theme(axis.text.x=element_text(size=10, vjust = 0.5, angle = 90, hjust=1, face = 'bold'), 
        axis.text.y=element_text(vjust=0.2, face = 'bold.italic'),
        plot.background=element_blank(),panel.border=element_rect(fill=NA),
        plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),
        plot.title=element_text(hjust=0,size=14,face="bold"),
        panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank")) +
  labs(x=NULL, y=NULL)

mets.associations %>% 
  group_by(type, variable) %>% 
  summarise_at('rho', median) %>% 
  # filter(abs(rho)>0.3) %>% 
  ggplot(aes(variable, type, fill = rho)) + 
  geom_tile(colour="white",size=0.2) + scale_fill_gradient2(low='blue', high = 'darkred') +
  theme(axis.text.x=element_text(size=10, vjust = 0.5, angle = 90, hjust=1, face = 'bold'), 
        axis.text.y=element_text(vjust=0.2, face = 'bold.italic'),
        plot.background=element_blank(),panel.border=element_rect(fill=NA),
        plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),
        plot.title=element_text(hjust=0,size=14,face="bold"),
        panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank")) +
  labs(x=NULL, y=NULL)

# limit to significant metabolites
mets.associations2.red <- mets.associations2[mets.associations2$metabolite%in%keep.mets,]
mets.associations2.red %>% 
  group_by(type, variable) %>% 
  summarise_at('rho', median) %>% 
  # filter(abs(rho)>0.3) %>% 
  ggplot(aes(variable, type, fill = rho)) + 
  geom_tile(colour="white",size=0.2) + scale_fill_gradient2(low='blue', high = 'darkred') +
  theme(axis.text.x=element_text(size=10, vjust = 0.5, angle = 90, hjust=1, face = 'bold'), 
        axis.text.y=element_text(vjust=0.2, face = 'bold.italic'),
        plot.background=element_blank(),panel.border=element_rect(fill=NA),
        plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),
        plot.title=element_text(hjust=0,size=14,face="bold"),
        panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank")) +
  labs(x=NULL, y=NULL)

mets.associations2.red %>% 
  group_by(path, variable) %>% 
  summarise_at('rho', median) %>% 
  # filter(abs(rho)>0.3) %>% 
  ggplot(aes(variable, path, fill = rho)) + 
  geom_tile(colour="white",size=0.2) + scale_fill_gradient2(low='blue', high = 'darkred') +
  theme(axis.text.x=element_text(size=10, vjust = 0.5, angle = 90, hjust=1, face = 'bold'), 
        axis.text.y=element_text(vjust=0.2, face = 'bold.italic'),
        plot.background=element_blank(),panel.border=element_rect(fill=NA),
        plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),
        plot.title=element_text(hjust=0,size=14,face="bold"),
        panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank")) +
  labs(x=NULL, y=NULL)

## progression correlations ----
delta.values <- readRDS('data/DeltaofLogValues.rds')
mets.residuals.prog <- mets.residuals[row.names(mets.residuals)%in%delta.values$studyID,]
identical(row.names(mets.residuals.prog), delta.values$studyID)

cor.results <- cor.associations(delta.values, mets.residuals.prog, names(delta.values)[6:24])

results.dfs.sig <- lapply(cor.results, function(x) {
  x <- x[x$fdr < .1, ]
  return(data.frame(
    'variable' = x$variable,
    'metabolite' = x$metabolite,
    'rho' = x$rho,
    'fdr' = x$fdr
  ))
})

mets.associations.deltavals <-  do.call(rbind.data.frame, results.dfs.sig)
mets.associations.deltavals$type <- sapply(mets.associations.deltavals$metabolite, function(a){
  annotations[grep(paste0('^',a,'$'), row.names(annotations)),]$SUPER_PATHWAY
}
)
mets.associations.deltavals$path <- sapply(mets.associations.deltavals$metabolite, function(a){
  annotations[grep(paste0('^',a,'$'), row.names(annotations)),]$SUB_PATHWAY
}
)
mets.associations.deltavals$name <- sapply(mets.associations.deltavals$metabolite, function(a){
  annotations[grep(paste0('^',a,'$'), row.names(annotations)),]$SHORT_NAME
}
)

mets.associations.deltavals$rho <- as.numeric(mets.associations.deltavals$rho)

# mets.associations.deltavals <- mets.associations.deltavals[abs(mets.associations.deltavals$rho)>0.3,]

mets.associations.deltavals <- arrange(mets.associations.deltavals, type, path)
mets.associations.deltavals$variable <- factor(mets.associations.deltavals$variable, levels= names(delta.values)[6:24])
order.rows <- cbind('name'=mets.associations.deltavals$name, 'type' = mets.associations.deltavals$type,'path' = mets.associations.deltavals$path)
order.rows <- order.rows[!duplicated(order.rows),]
order.rows <- as.data.frame(order.rows)
order.rows <- arrange(order.rows, type, path) 
mets.associations.deltavals$name <- factor(mets.associations.deltavals$name, levels= order.rows$name)

pdf('deltavalues_metabolites.pdf', height = 25, width = 10)
ggplot(mets.associations.deltavals, aes(variable, name, fill = rho)) + 
  geom_tile(colour="white",size=0.2) + scale_fill_gradient2(low='blue', high = 'darkred') +
  theme(axis.text.x=element_text(size=10, vjust = 0.5, angle = 90, hjust=1, face = 'bold'), 
        axis.text.y=element_text(vjust=0.2, face = 'bold.italic'),
        plot.background=element_blank(),panel.border=element_rect(fill=NA),
        plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),
        plot.title=element_text(hjust=0,size=14,face="bold"),
        panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank")) +
  labs(x=NULL, y=NULL)
dev.off()
write.table(order.rows, 'deltavalues_metabolites_order.txt', sep = '\t')

mets.associations.deltavals %>% 
  group_by(path, variable) %>% 
  summarise_at('rho', median) %>% 
  mutate(path2=factor(path, levels = order.rows[!duplicated(order.rows$path),]$path)) %>% 
  mutate(variable=factor(variable, levels= names(delta.values)[6:24])) %>% 
  ggplot(aes(variable, path2, fill = rho)) + 
  geom_tile(colour="white",size=0.2) + scale_fill_gradient2(low='blue', high = 'darkred') +
  theme(axis.text.x=element_text(size=10, vjust = 0.5, angle = 90, hjust=1, face = 'bold'), 
        axis.text.y=element_text(vjust=0.2, face = 'bold.italic'),
        plot.background=element_blank(),panel.border=element_rect(fill=NA),
        plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),
        plot.title=element_text(hjust=0,size=14,face="bold"),
        panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank")) +
  labs(x=NULL, y=NULL)


## Associations of total or partially bacteria-related metabolites -----
bact.mets <- read.table('explainedvarianceLASSOmodels.txt', header=T) ## this file comes from the LASSO modelling, needs to be generated previously!
bact.mets <- bact.mets[bact.mets$data.type == 'metaG',]
bact.mets <- unique(bact.mets$metabolite)

mets.associations2.bact <- mets.associations2[mets.associations2$metabolite%in%bact.mets,]
mets.associations2.bact %>% 
  group_by(type, variable) %>% 
  summarise_at('rho', median) %>% 
  # filter(abs(rho)>0.3) %>% 
  ggplot(aes(variable, type, fill = rho)) + 
  geom_tile(colour="white",size=0.2) + scale_fill_gradient2(low='blue', high = 'darkred') +
  theme(axis.text.x=element_text(size=10, vjust = 0.5, angle = 90, hjust=1, face = 'bold'), 
        axis.text.y=element_text(vjust=0.2, face = 'bold.italic'),
        plot.background=element_blank(),panel.border=element_rect(fill=NA),
        plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),
        plot.title=element_text(hjust=0,size=14,face="bold"),
        panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank")) +
  labs(x=NULL, y=NULL,title = 'Totally/Partially bacterial metabolites associations to \nmetabolite families clustered into metabolite families')

mets.associations2.bact %>% 
  group_by(path, variable) %>% 
  summarise_at('rho', median) %>% 
  # filter(abs(rho)>0.3) %>% 
  ggplot(aes(variable, path, fill = rho)) + 
  geom_tile(colour="white",size=0.2) + scale_fill_gradient2(low='blue', high = 'darkred') + #coord_flip() +
  theme(axis.text.x=element_text(size=10, vjust = 0.5, angle = 90, hjust=1, face = 'bold'), 
        axis.text.y=element_text(vjust=0.2, face = 'bold.italic'),
        plot.background=element_blank(),panel.border=element_rect(fill=NA),
        plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),
        plot.title=element_text(hjust=0,size=14,face="bold"),
        panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank")) +
  labs(x=NULL, y=NULL, title = 'Totally/Partially bacterial metabolites associations to metabolic pathways\nclustered into pathways')

ggplot(mets.associations2.bact, aes(variable, name, fill = rho)) + 
  geom_tile(colour="white",size=0.2) + scale_fill_gradient2(low='blue', high = 'darkred') +
  theme(axis.text.x=element_text(size=10, vjust = 0.5, angle = 90, hjust=1, face = 'bold'), 
        axis.text.y=element_text(vjust=0.2, face = 'bold.italic'),
        plot.background=element_blank(),panel.border=element_rect(fill=NA),
        plot.margin=margin(0.7,0.4,0.1,0.2,"cm"),
        plot.title=element_text(hjust=0,size=14,face="bold"),
        panel.grid.major = element_line(linetype = "blank"), panel.grid.minor = element_line(linetype = "blank")) +
  labs(x=NULL, y=NULL, title = 'Totally/Partially bacterial metabolites associations to biochemical variables')

