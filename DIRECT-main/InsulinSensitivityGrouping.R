########################################################################################################################
##############################  Supplementary: Insulin sensitivity stratifications #####################################
########################################################################################################################

## directory
setwd('/home/scratch/margar')

## libraries
library(tidyverse); library(tidyr); library(dplyr); library(reshape2); library(ggplot2); library(ggbeeswarm); library(ggpubr)
library(factoextra); library(phyloseq); library(effsize); library(cowplot); library(circlize); library(caret); library(grid)

## upload and preprocess data ----
biochem.bsl <- read.table('data/curatedmetadata/WP2.1_BL_MetadataWithInfl.txt', sep = '\t', header=T)

## Stumvoll, Matsuda and OGIS-2h distributions ----
plot_grid(ggplot(biochem.bsl, aes(Stumvoll)) + geom_density(),
          ggplot(biochem.bsl, aes(Matsuda)) + geom_density(),
          ggplot(biochem.bsl, aes(OGIS.2h)) + geom_density(), ncol = 3)

## We use OGIS2h to stratify the samples (need to impute 4 samples!)
preProcValues <- preProcess(biochem.bsl,
                            method = c('knnImpute'),
                            k = 25,
                            knnSummary = mean)
imputed.bsl.info <- predict(preProcValues, biochem.bsl, na.action=na.pass)

procNames <- data.frame(col = names(preProcValues$mean), mean = preProcValues$mean, sd = preProcValues$std)
for(i in procNames$col){
  imputed.bsl.info[i] <- imputed.bsl.info[i]*preProcValues$std[i]+preProcValues$mean[i] 
}

biochem.bsl <- imputed.bsl.info # replace biochem.bsl DF with the imputed dataset
biochem.bsl$OGIS.strat <- factor(ntile(biochem.bsl$OGIS.2h, 3))
biochem.bsl$OGIS.strat <- plyr::revalue(biochem.bsl$OGIS.strat, c('1'='Low.OGIS','2'='Mid.OGIS','3'='High.OGIS')) #stratify in tertiles 

plot_grid(
  ggplot(biochem.bsl, aes(OGIS.strat, OGIS.2h, group = interaction(OGIS.strat))) +
  geom_violin(alpha = .1, aes(fill=OGIS.strat, color=OGIS.strat)) +
  geom_quasirandom(alpha = .6, dodge.width = .9, aes(color = OGIS.strat)) +
  geom_boxplot(width = .1, alpha = .7, outlier.colour = NA, position = position_dodge(.9), fill='white') +
  scale_color_manual(values = c('cadetblue2', 'cornflowerblue', 'blue')) +
  scale_fill_manual(values = c('cadetblue2', 'cornflowerblue', 'blue')) + 
  stat_compare_means() + 
  stat_compare_means(comparisons = list(c('Low.OGIS', 'Mid.OGIS'), 
                                        c('Low.OGIS', 'High.OGIS'), c('Mid.OGIS', 'High.OGIS'))) +
  theme(panel.grid.major = element_line(colour = "gray85",linetype = "dashed"), 
        axis.text = element_text(size = 11, face = "bold", colour = "black"), 
        axis.text.x = element_text(angle = 90, vjust = .5),
        panel.background = element_rect(fill = NA), legend.title = element_text(face='bold'),
        legend.text = element_text(face='bold'), strip.background = element_rect(fill=NA),
        strip.text = element_text(face='bold', size =12), panel.border = element_rect(fill=NA), legend.position='none') +
  labs(x = NULL, y = NULL, title = 'OGIS 2h'),
  ggplot(biochem.bsl, aes(OGIS.strat, Stumvoll, group = interaction(OGIS.strat))) +
    geom_violin(alpha = .1, aes(fill=OGIS.strat, color=OGIS.strat)) +
    geom_quasirandom(alpha = .6, dodge.width = .9, aes(color = OGIS.strat)) +
    geom_boxplot(width = .1, alpha = .7, outlier.colour = NA, position = position_dodge(.9), fill='white') +
    scale_color_manual(values = c('cadetblue2', 'cornflowerblue', 'blue')) +
    scale_fill_manual(values = c('cadetblue2', 'cornflowerblue', 'blue')) + 
    stat_compare_means() + 
    stat_compare_means(comparisons = list(c('Low.OGIS', 'Mid.OGIS'), 
                                          c('Low.OGIS', 'High.OGIS'), c('Mid.OGIS', 'High.OGIS'))) +
    theme(panel.grid.major = element_line(colour = "gray85",linetype = "dashed"), 
          axis.text = element_text(size = 11, face = "bold", colour = "black"), 
          axis.text.x = element_text(angle = 90, vjust = .5),
          panel.background = element_rect(fill = NA), legend.title = element_text(face='bold'),
          legend.text = element_text(face='bold'), strip.background = element_rect(fill=NA),
          strip.text = element_text(face='bold', size =12), panel.border = element_rect(fill=NA), legend.position='none') +
    labs(x = NULL, y = NULL, title = "Stumvoll's IR"),
  ggplot(biochem.bsl, aes(OGIS.strat, Matsuda, group = interaction(OGIS.strat))) +
    geom_violin(alpha = .1, aes(fill=OGIS.strat, color=OGIS.strat)) +
    geom_quasirandom(alpha = .6, dodge.width = .9, aes(color = OGIS.strat)) +
    geom_boxplot(width = .1, alpha = .7, outlier.colour = NA, position = position_dodge(.9), fill='white') +
    scale_color_manual(values = c('cadetblue2', 'cornflowerblue', 'blue')) +
    scale_fill_manual(values = c('cadetblue2', 'cornflowerblue', 'blue')) + 
    stat_compare_means() + 
    stat_compare_means(comparisons = list(c('Low.OGIS', 'Mid.OGIS'), 
                                          c('Low.OGIS', 'High.OGIS'), c('Mid.OGIS', 'High.OGIS'))) +
    theme(panel.grid.major = element_line(colour = "gray85",linetype = "dashed"), 
          axis.text = element_text(size = 11, face = "bold", colour = "black"), 
          axis.text.x = element_text(angle = 90, vjust = .5),
          panel.background = element_rect(fill = NA), legend.title = element_text(face='bold'),
          legend.text = element_text(face='bold'), strip.background = element_rect(fill=NA),
          strip.text = element_text(face='bold', size =12), panel.border = element_rect(fill=NA), legend.position='none') +
    labs(x = NULL, y = NULL, title = "Matsuda's IR"), ncol = 3)

## differences in gene richness stratifications ----
tertiles.d <- biochem.bsl %>% 
  group_by(OGIS.strat, DSGARLTertileGroup) %>% 
  tally() %>% 
  pivot_wider(names_from=DSGARLTertileGroup, values_from=n) ## count number of individuals in each tertile per each IGR type
tertiles.d <- as.data.frame(t(tertiles.d))
names(tertiles.d) <- tertiles.d[1,]
tertiles.d <- tertiles.d[-1,]
#chisq test on counts#
tertiles.d <- sapply(tertiles.d, as.numeric) #chisquare test for the 2 groups
chisq.test(tertiles.d) # X2 26.729, pvalue 2.255e-05, clearly associated
#
tertiles.d <- data.frame(tertiles.d, 'GR' = c('EGC', 'HGC', 'LGC'))
tertiles.d$GR <- factor(tertiles.d$GR, levels = c('LGC', 'HGC', 'EGC'))
tertiles.d.m <- melt(tertiles.d)
ggplot(tertiles.d.m, aes(GR, value, fill = GR)) + geom_col() +
  facet_wrap(~variable) + scale_fill_manual(values=c('yellow', 'gray75', 'blue'))

## differences in enterotypes stratifications ----
entero.d <- biochem.bsl %>% 
  group_by(OGIS.strat, enterotype) %>% 
  tally() %>% 
  pivot_wider(names_from=enterotype, values_from=n) ## count number of individuals in each tertile per each IGR type
entero.d <- as.data.frame(t(entero.d))
names(entero.d) <- entero.d[1,]
entero.d <- entero.d[-1,]
#chisq test on counts#
entero.d <- sapply(entero.d, as.numeric) #chisquare test for the 2 groups
chisq.test(entero.d) # X2 30.184 pvalue 3.627e-05, clearly associated
entero.d <- as.data.frame(entero.d)
entero.d$entero <- c('Bact1', 'Bact2', 'Prev', 'Rum')
entero.d.m <- melt(entero.d)
ggplot(entero.d.m, aes(entero, value, fill=entero)) + geom_col() + facet_wrap(~variable, scales='free')

## gene richness (continous) and bacterial diversity ----
richness <- biochem.bsl %>% 
  select(c('OGIS.strat', 'DSGeneAbundRichness', 'DSMGSAbundShannon'))
plot_grid(
  ggplot(richness, aes(OGIS.strat, DSGeneAbundRichness)) + geom_boxplot() + stat_compare_means(),
  ggplot(richness, aes(OGIS.strat, DSMGSAbundShannon)) + geom_boxplot() + stat_compare_means())

## compare biochemical variables between Insulin stratifications -----
biochem.variables <- data.frame('OGIS.strat'= biochem.bsl$OGIS.strat, biochem.bsl[,16:57]) # extract only biochemical variables

biochem.stats <- biochem.variables %>% 
  group_by(OGIS.strat) %>% 
  summarise_all(list(mean=mean, median=median, sd=sd), na.rm=T) %>% 
  pivot_longer(ends_with(c("_mean", "_sd", "_median")),
               names_to = c("name", ".value"),
               names_sep = "_") %>% 
  pivot_wider(names_from = OGIS.strat,
              names_sep = "-",
              values_from = c(mean, sd, median)) ## compute mean, median and SD per IGR subtype


p.vals <- list()
for (i in seq(from=2, to=43)){
  model <- aov(biochem.variables[,i]~biochem.variables$OGIS.strat, data = biochem.variables)
  anova.p <- summary(model)[[1]][['Pr(>F)']][1]
  tuks <- TukeyHSD(model)
  tuks.p <- tuks$`biochem.variables$OGIS.strat`[,4]
  p.vals[[i]] <- c('anova'=anova.p, tuks.p)
}

biochem.variables.stats.p <- cbind(biochem.stats, do.call(rbind.data.frame, p.vals)) # attach p.values to stats
names(biochem.variables.stats.p) <- c(names(biochem.stats), 
                                    c('anova', 'Mid-Low', 'High-Low', 'High-Mid'))
biochem.variables.stats.p$aov.fdr  <- p.adjust(biochem.variables.stats.p$anova, method='fdr') #multitest adjust

# write.table(biochem.variables.stats.p, 'InsulinSens/ResultsBiochemInsSensComparison.txt', sep ='\t')

biochem.variables.m <- melt(biochem.variables)
ggplot(biochem.variables.m, aes(OGIS.strat, value), group = interaction(OGIS.strat)) + 
  geom_violin(alpha = .1, aes(fill=OGIS.strat, color=OGIS.strat)) +
  geom_quasirandom(alpha = .6, dodge.width = .9, aes(color = OGIS.strat)) +
  geom_boxplot(width = .1, alpha = .7, outlier.colour = NA, position = position_dodge(.9), fill='white')  +
  facet_wrap(~variable, scales = 'free') + 
  stat_compare_means(label = 'p.signif')

diff.biochem <- biochem.variables[,names(biochem.variables)%in%biochem.variables.stats.p[biochem.variables.stats.p$aov.fdr<.1,]$name] ## plot only the significant ones
diff.biochem$OGIS.strat <- biochem.variables$OGIS.strat
diff.biochem.m <- melt(diff.biochem)
ggplot(diff.biochem.m, aes(OGIS.strat, value), group = interaction(OGIS.strat)) + 
  geom_violin(alpha = .1, aes(fill=OGIS.strat, color=OGIS.strat)) +
  geom_quasirandom(alpha = .6, dodge.width = .9, aes(color = OGIS.strat)) +
  geom_boxplot(width = .1, alpha = .7, outlier.colour = NA, position = position_dodge(.9), fill='white')  +
  facet_wrap(~variable, scales = 'free') + 
  stat_compare_means(label = 'p.signif') +
  scale_color_manual(values = c('cadetblue2', 'cornflowerblue', 'blue')) +
  scale_fill_manual(values = c('cadetblue2', 'cornflowerblue', 'blue')) +
  theme(legend.position = 'bottom')

biochem.variables.pca <- na.omit(biochem.variables[,2:43])
biochem.variables.pca <- biochem.variables.pca[, !(names(biochem.variables.pca)%in%c('Stumvoll', 'Matsuda', 'OGIS.2h'))]
biochem.variables.pca <- prcomp(biochem.variables.pca, center=T, scale. = T) #PCA 
fviz_screeplot(biochem.variables.pca)
row.names(get_pca_ind(biochem.variables.pca)$coord)
fviz_pca_biplot(biochem.variables.pca, geom.ind = 'point', 
                col.ind = biochem.variables[row.names(biochem.variables)%in%row.names(get_pca_ind(biochem.variables.pca)$coord),]$OGIS.strat, addEllipses = T)
fviz_pca_ind(biochem.variables.pca, geom.ind = 'point', 
                col.ind = biochem.variables[row.names(biochem.variables)%in%row.names(get_pca_ind(biochem.variables.pca)$coord),]$OGIS.strat, addEllipses = T)
fviz_pca_var(biochem.variables.pca)


## compare inflammation markers between groups ----
inflammatory <- data.frame('OGIS.strat'= biochem.bsl$OGIS.strat, biochem.bsl[,58:72]) ## extract only inflammatory markers
inflammatory.stats <- inflammatory %>% 
  group_by(OGIS.strat) %>% 
  summarise_all(list(mean=mean, median=median, sd=sd), na.rm=T) %>% 
  pivot_longer(ends_with(c("_mean", "_sd", "_median")),
               names_to = c("name", ".value"),
               names_sep = "_") %>% 
  pivot_wider(names_from = OGIS.strat,
              names_sep = "-",
              values_from = c(mean, sd, median)) ## compute mean, median and SD per IGR subtype


p.vals <- list()
for (i in seq(from=2, to=16)){
  model <- aov(inflammatory[,i]~inflammatory$OGIS.strat, data = inflammatory)
  anova.p <- summary(model)[[1]][['Pr(>F)']][1]
  tuks <- TukeyHSD(model)
  tuks.p <- tuks$`inflammatory$OGIS.strat`[,4]
  p.vals[[i]] <- c('anova'=anova.p, tuks.p)
}

inflammatory.stats.p <- cbind(inflammatory.stats, do.call(rbind.data.frame, p.vals)) # attach p.values to stats
names(inflammatory.stats.p) <- c(names(inflammatory.stats), 
                                      c('anova', 'Mid-Low', 'High-Low', 'High-Mid'))
inflammatory.stats.p$aov.fdr  <- p.adjust(inflammatory.stats.p$anova, method='fdr') #multitest adjust

# write.table(inflammatory.stats.p, 'InsulinSens/ResultsInflammatoryMarkersInsulinSecromparison.txt', sep ='\t')

inflammatory.m <- melt(inflammatory)
ggplot(inflammatory.m, aes(OGIS.strat, value), group = interaction(OGIS.strat)) + 
  geom_violin(alpha = .1, aes(fill=OGIS.strat, color=OGIS.strat)) +
  geom_quasirandom(alpha = .6, dodge.width = .9, aes(color = OGIS.strat)) +
  geom_boxplot(width = .1, alpha = .7, outlier.colour = NA, position = position_dodge(.9), fill='white')  +
  facet_wrap(~variable, scales = 'free') + 
  stat_compare_means(label = 'p.signif') +
  scale_color_manual(values = c('cadetblue2', 'cornflowerblue', 'blue')) +
  scale_fill_manual(values = c('cadetblue2', 'cornflowerblue', 'blue')) +
  theme(legend.position = c(.9,.1))

## Metagenome differences at baseline between insulin sensitivity stratifications ----
qmp.bsl <- read.table('data/WP2.1_775_IGR_MGS-QMP_igrHA_Ent_M.csv', header = T, row.names=1) # MGS data
qmp.bsl <- qmp.bsl[,7:727]

tax <- read.table('data/taxonomy.tsv', header = T, sep ='\t', row.names=1) # read taxonomy for Phyloseq
row.names(tax) <- gsub(':', '\\.', row.names(tax))
tax.m <- as.matrix(tax); colnames(tax.m) <- names(tax); rownames(tax.m) <- row.names(tax)

qmp.bsl <- qmp.bsl[row.names(qmp.bsl)%in%biochem.bsl$samplerenamed,]
physeq.bsl <- phyloseq(otu_table(qmp.bsl, taxa_are_rows=F), tax_table(tax.m)) ## cluster QMPs for genus taxonomical levels
physeq.bsl.genus <- tax_glom(physeq.bsl, taxrank = rank_names(physeq.bsl)[6])

# compute and plot beta diversity
jsd.dist.mgs <- distance(physeq.bsl, method = 'jsd') ## Jenssen-shannon
pcoa.jsd.bsl <- ordinate(physeq.bsl, 'PCoA', distance = jsd.dist.mgs)
pcoa.df <- data.frame(pcoa.jsd.bsl$vectors, 'OGIS.strat' = biochem.bsl$OGIS.strat)
ggplot(pcoa.df, aes(Axis.1, Axis.2, color=OGIS.strat)) + geom_point() +
  geom_hline(yintercept = 0, linetype = 'dashed') + geom_vline(xintercept=0, linetype = 'dashed') +
  stat_ellipse() + scale_color_manual(values=c('cadetblue2', 'cornflowerblue', 'blue')) +
  theme(panel.grid.major = element_line(colour = "gray85", linetype = "dashed"), axis.title = element_text(size = 17, face = "bold"), 
        axis.text = element_text(size = 13,  face = "bold"), legend.text = element_text(size = 12, face = "bold"), 
        legend.title = element_text(face = "bold"), panel.background = element_rect(fill = NA), 
        legend.key = element_rect(fill = NA), legend.background = element_rect(fill = NA), panel.border = element_rect(fill=NA)) +
  labs(x = "PCoA 1 (expl. variance 13.4%)", y = "PCoA 2 (expl. variance 10.2%)", colour = "Insulin sensitivity type")
vegan::adonis(jsd.dist.mgs~OGIS.strat, data=biochem.bsl)

byc.dist.mgs <- distance(physeq.bsl, method = 'bray') ## Bray curtis
pcoa.byc.bsl <- ordinate(physeq.bsl, 'PCoA', distance = byc.dist.mgs)
pcoa.df <- data.frame(pcoa.byc.bsl$vectors, 'OGIS.strat' = biochem.bsl$OGIS.strat)
ggplot(pcoa.df, aes(Axis.1, Axis.2, color=OGIS.strat)) + geom_point() +
  geom_hline(yintercept = 0, linetype = 'dashed') + geom_vline(xintercept=0, linetype = 'dashed') +
  stat_ellipse() + scale_color_manual(values=c('cadetblue2', 'cornflowerblue', 'blue')) +
  theme(panel.grid.major = element_line(colour = "gray85", linetype = "dashed"), axis.title = element_text(size = 17, face = "bold"), 
        axis.text = element_text(size = 13,  face = "bold"), legend.text = element_text(size = 12, face = "bold"), 
        legend.title = element_text(face = "bold"), panel.background = element_rect(fill = NA), 
        legend.key = element_rect(fill = NA), legend.background = element_rect(fill = NA), panel.border = element_rect(fill=NA)) +
  labs(x = "PCoA 1 (expl. variance 9.62%)", y = "PCoA 2 (expl. variance 6.78%)", colour = "Insulin sensitivity type")
vegan::adonis(byc.dist.mgs~OGIS.strat, data=biochem.bsl)

# MGS model
biochem.bsl$samplerenamed == row.names(qmp.bsl)
mgs.model <- data.frame('StudyID' = biochem.bsl$StudyID, 'center' = biochem.bsl$CenterID, # generate DF for the modeling step
                        'gender' = biochem.bsl$Gender, 'age'=biochem.bsl$Age, 
                        'OGIS.strat' = biochem.bsl$OGIS.strat, qmp.bsl)
mgs.model <- mgs.model[mgs.model$OGIS.strat!='Mid.OGIS',]
mgs.model$OGIS.strat <- droplevels(mgs.model$OGIS.strat)

result.model <- list()
for (x in seq(from=6, to=726)){
  print(names(mgs.model[x]))
  model <- data.frame('mgs' = mgs.model[,x], mgs.model[,1:5])
  mod.a <- lm(mgs~OGIS.strat+gender+center+age, data = model, na.action = na.omit) # fixed effects: age sex and center
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['OGIS.stratHigh.OGIS', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['OGIS.stratHigh.OGIS', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$OGIS.strat)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <- c(names(mgs.model)[x],p.val, delta, delta.low, delta.up)
  result.model[[x]] <- result
}

mgs.lm <- do.call(rbind.data.frame, result.model)
names(mgs.lm) <- c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
mgs.lm$fdr <- p.adjust(mgs.lm$pval, method='fdr') #pvalue adjustment
mgs.lm[mgs.lm$fdr<.1,] ## at fdr 10% no differences were observed for MGS between IGT and IFG
hist(mgs.lm$fdr)
mgs.lm[mgs.lm$pval<.05,]$mgs
tax[row.names(tax)%in%mgs.lm[mgs.lm$pval<.05,]$mgs,]$species
mgs.lm.pval.sig <- cbind(mgs.lm[mgs.lm$pval<.05,], 'species' = tax[row.names(tax)%in%mgs.lm[mgs.lm$pval<.05,]$mgs,]$species)

# genus level
genera.counts <- as.data.frame(physeq.bsl.genus@otu_table)
biochem.bsl$samplerenamed == row.names(genera.counts)
genus.model <- data.frame('StudyID' = biochem.bsl$StudyID, 'center' = biochem.bsl$CenterID, # generate DF for the modeling step
                          'gender' = biochem.bsl$Gender, 'age'=biochem.bsl$Age, 
                          'OGIS.strat' = biochem.bsl$OGIS.strat, genera.counts)

genus.model <- genus.model[genus.model$OGIS.strat!='Mid.OGIS',]
genus.model$OGIS.strat <- droplevels(genus.model$OGIS.strat)

result.model <- list()
for (x in seq(from=6, to=128)){
  print(names(genus.model[x]))
  model <- data.frame('genus' = genus.model[,x], genus.model[,1:5])
  mod.a <- lm(genus~OGIS.strat+gender+center+age, data = model, na.action = na.omit) # fixed effects: age sex and center
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['OGIS.stratHigh.OGIS', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['OGIS.stratHigh.OGIS', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$genus, model$OGIS.strat)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <- c(names(genus.model)[x],p.val, delta, delta.low, delta.up)
  result.model[[x]] <- result
}

genus.lm <- do.call(rbind.data.frame, result.model)
names(genus.lm) <- c('genus', 'pval', 'delta', 'delta.low', 'delta.up')
genus.lm$fdr <- p.adjust(genus.lm$pval, method='fdr') #pvalue adjustment
genus.lm[genus.lm$fdr<.1,] ## at fdr 10% no differences were observed for genus level between IGT and IFG
hist(genus.lm$fdr)
genus.lm[genus.lm$pval<.05,]$genus
tax[row.names(tax)%in%genus.lm[genus.lm$pval<.05,]$genus,]$genus
genus.lm.pval.sig <- cbind(genus.lm[genus.lm$pval<.05,], 'genus2' = tax[row.names(tax)%in%genus.lm[genus.lm$pval<.05,]$genus,]$genus)

## GMMs differences between Insulin sensitivity stratifications ----
GMMs <- read.table('data/WP2.1_775_IGR_GMM-QMP_igrHA_Ent_M.csv', header = T, row.names = 1)
GMMs <- GMMs[,7:103]
biochem.bsl$samplerenamed == row.names(GMMs)
gmm.model <- data.frame('StudyID' = biochem.bsl$StudyID, 'center' = biochem.bsl$CenterID, # generate DF for the modeling step
                        'gender' = biochem.bsl$Gender, 'age'=biochem.bsl$Age, 
                        'OGIS.strat' = biochem.bsl$OGIS.strat, GMMs)
gmm.model <- gmm.model[gmm.model$OGIS.strat!='Mid.OGIS',]
gmm.model$OGIS.strat <- droplevels(gmm.model$OGIS.strat)

result.model <- list()
for (x in seq(from=6, to=102)){
  print(names(gmm.model[x]))
  model <- data.frame('GMM' = gmm.model[,x], gmm.model[,1:5])
  mod.a <- lm(GMM~OGIS.strat+gender+center+age, data = model, na.action = na.omit) # fixed effects: age sex and center
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['OGIS.stratHigh.OGIS', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['OGIS.stratHigh.OGIS', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$GMM, model$OGIS.strat)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <- c(names(gmm.model)[x],p.val, delta, delta.low, delta.up)
  result.model[[x]] <- result
}

GMM.lm <- do.call(rbind.data.frame, result.model)
names(GMM.lm) <- c('GMM', 'pval', 'delta', 'delta.low', 'delta.up')
GMM.lm$fdr <- p.adjust(GMM.lm$pval, method='fdr') #pvalue adjustment
GMM.lm[GMM.lm$fdr<.1,] ## at fdr 10% no differences were observed for GMMs between IGT and IFG
hist(GMM.lm$fdr)

GMM.lm.pval.sig <- GMM.lm[GMM.lm$pval<.05,]

# save results
write.table(plyr::rbind.fill(mgs.lm.pval.sig, genus.lm.pval.sig, GMM.lm.pval.sig),
            'InsulinSens/ReducedDatasetMGSGenusGMM_linearmodelsInsulinSecr.txt', sep = '\t')

## structural variants ----
# data upload
deletions <- read.table('SV-followup/dsgv_baseline.csv', header = T, sep = ',', row.names = 1)
row.names(deletions) == biochem.bsl$SampleID
deletions <- deletions[match(biochem.bsl$SampleID, row.names(deletions)),] # match order between DFs

variations <- read.table('SV-followup/vsgv_baseline.csv', header = T, sep = ',', row.names = 1)
row.names(deletions) == row.names(variations)
variations <- variations[match(row.names(deletions), row.names(variations)),] # match order between DFs

SVs.per.sample <- data.frame('deletions' = rowSums(!is.na(deletions)&deletions==0), 'variations' = rowSums(!is.na(variations))) # total number of deletions and variable regions per sample
SVs.per.sample$OGIS.strat <- biochem.bsl$OGIS.strat

p1<-ggplot(SVs.per.sample, aes(deletions, variations, color=OGIS.strat)) + geom_point() + geom_smooth(method='lm') + # correlations between deletions and variable regions per sample
  scale_color_manual(values = c('cadetblue2', 'cornflowerblue', 'blue')) + 
  scale_fill_manual(values = c('cadetblue2', 'cornflowerblue', 'blue')) +
  stat_cor(label.x = 400,) +
  theme(axis.line = element_line(size = 0.6), axis.title = element_text(size = 13, face = "bold"), axis.text = element_text(face = "bold"), 
        legend.text = element_text(face = "bold"), legend.title = element_text(face = "bold"), panel.background = element_rect(fill = NA), 
        panel.border = element_rect(fill=NA, colour = 'black', size=.5), legend.key = element_rect(fill = NA), 
        legend.background = element_rect(fill = NA), 
        legend.position = c(0.94, 0.15)) + labs(x = "Number of deletions", y = "Number of variable SVs", color=NULL)


p1 <- ggExtra::ggMarginal(  # correlations
  p = p1,
  type = 'density',
  margins = 'both',
  size = 5,
  groupColour = TRUE,
  groupFill = FALSE,
  alpha=.3
)

p2 <- ggplot(SVs.per.sample, aes(OGIS.strat, deletions, color=OGIS.strat)) +  # boxplots for the deletions
  geom_boxplot(outlier.colour = NA, lwd = 1.4) +
  geom_quasirandom(alpha=.4) + 
  stat_compare_means()+
  scale_color_manual(values = c('cadetblue2', 'cornflowerblue', 'blue')) + 
  scale_fill_manual(values = c('cadetblue2', 'cornflowerblue', 'blue')) + 
  theme(panel.grid.major = element_line(colour = "gray85", linetype = "dashed"), 
        axis.title = element_text(size = 18, face = "bold"), 
        axis.text = element_text(size = 15, face = "bold"), 
        panel.background = element_rect(fill = NA), legend.position = "none", 
        panel.border = element_rect(fill=NA)) +
  labs(x = NULL, y = "Number of deletions", colour = NULL) 

p3 <- ggplot(SVs.per.sample, aes(OGIS.strat, variations, color=OGIS.strat)) + # boxplots for the variable regions
  geom_boxplot(outlier.colour = NA, lwd = 1.4) +
  geom_quasirandom(alpha=.4) + 
  stat_compare_means()+
  scale_color_manual(values = c('cadetblue2', 'cornflowerblue', 'blue')) + 
  scale_fill_manual(values = c('cadetblue2', 'cornflowerblue', 'blue')) + 
  theme(panel.grid.major = element_line(colour = "gray85", linetype = "dashed"), 
        axis.title = element_text(size = 18, face = "bold"), 
        axis.text = element_text(size = 15, face = "bold"), 
        panel.background = element_rect(fill = NA), legend.position = "none", 
        panel.border = element_rect(fill=NA)) +
  labs(x = NULL, y = "Number of variable regions", colour = NULL) 

plot_grid(p1, plot_grid(p2, p3, ncol =1), ncol =2)  # put everything in one plot


# deletions associations to phenotype
biochem.bsl$SampleID == row.names(deletions)
deletions.model <- data.frame('StudyID' = biochem.bsl$StudyID, 'center' = biochem.bsl$CenterID, # generate DF for the modeling step
                              'gender' = biochem.bsl$Gender, 'age'=biochem.bsl$Age, 
                              'OGIS.strat' = biochem.bsl$OGIS.strat, deletions)
deletions.model <- deletions.model[deletions.model$OGIS.strat!='Mid.OGIS',]
deletions.model$OGIS.strat <- droplevels(deletions.model$OGIS.strat)

result.model <- list()
for (x in seq(from=6, to=7720)){
  print(names(deletions.model[x]))
  model <- data.frame('Deletion' = deletions.model[,x], deletions.model[,1:5])
  mod.a <- lm(Deletion~OGIS.strat+gender+center+age, data = model, na.action = na.omit) # fixed effects: age sex and center
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['OGIS.stratHigh.OGIS', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['OGIS.stratHigh.OGIS', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$Deletion, model$OGIS.strat)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <- c(names(deletions.model)[x],p.val, delta, delta.low, delta.up)
  result.model[[x]] <- result
}

deletions.lm <- do.call(rbind.data.frame, result.model)
names(deletions.lm) <- c('Deletion', 'pval', 'delta', 'delta.low', 'delta.up')
deletions.lm$fdr <- p.adjust(deletions.lm$pval, method='fdr') #pvalue adjustment
deletions.lm[deletions.lm$fdr<.1,] ## at fdr 10% no differences were observed for the deleted regions between IGT and IFG

# variations associations to phenotype
biochem.bsl$SampleID == row.names(variations)
variations.model <- data.frame('StudyID' = biochem.bsl$StudyID, 'center' = biochem.bsl$CenterID, # generate DF for the modeling step
                               'gender' = biochem.bsl$Gender, 'age'=biochem.bsl$Age, 
                               'OGIS.strat' = biochem.bsl$OGIS.strat, variations)
variations.model <- variations.model[variations.model$OGIS.strat!='Mid.OGIS',]
variations.model$OGIS.strat <- droplevels(variations.model$OGIS.strat)

result.model <- list()
for (x in seq(from=6, to=3104)){
  print(names(variations.model[x]))
  model <- data.frame('Variation' = variations.model[,x], variations.model[,1:5])
  mod.a <- lm(Variation~OGIS.strat+gender+center+age, data = model, na.action = na.omit) # fixed effects: age sex and center
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['OGIS.stratHigh.OGIS', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['OGIS.stratHigh.OGIS', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$Variation, model$OGIS.strat)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <- c(names(variations.model)[x],p.val, delta, delta.low, delta.up)
  result.model[[x]] <- result
}

variations.lm <- do.call(rbind.data.frame, result.model)
names(variations.lm) <- c('Variation', 'pval', 'delta', 'delta.low', 'delta.up')
variations.lm$fdr <- p.adjust(variations.lm$pval, method='fdr') #pvalue adjustment
variations.lm[variations.lm$fdr<.1,] ## at fdr 10% no differences were observed for the deleted regions between IGT and IFG  
write.table(variations.lm[variations.lm$fdr<.1,], 'InsulinSens/SignificantVariationAssocitions.txt', sep='\t')

## Progression and regression differences in high-low insulins individuals ----
biochem.mix <- read.table('data/curatedmetadata/WP2.1_BlFw-641_MetadataImputed_OnlyBIOCHEM.txt', header = T) # read the mix of baseline and followup biochemical data
biochem.mix$Gly.Cat.fw <- c(biochem.mix[biochem.mix$timepoint=='Follow-up', 'Gly.Cat.2'],
                            biochem.mix[biochem.mix$timepoint=='Follow-up', 'Gly.Cat.2']) # add column with only the glycaemic status at followup

ab <- biochem.bsl[biochem.bsl$samplerenamed%in%biochem.mix[biochem.mix$timepoint=='Baseline',]$sample.renamed,]
biochem.mix$OGIS.strat <- rep(ab$OGIS.strat, 2)
progression.summ <- biochem.mix %>% 
  filter(timepoint == 'Baseline') %>% 
  select(c('OGIS.strat', 'Gly.Cat.fw')) %>% 
  group_by(OGIS.strat, Gly.Cat.fw) %>%
  tally()                                     # count the number of individuals progressing to each status

progression.all <- progression.summ %>% 
  pivot_wider(names_from=Gly.Cat.fw, values_from=n) # prepare data for chisquare test
progression.all <- as.data.frame(t(progression.all))
names(progression.all) <- progression.all[1,]
progression.all <- progression.all[-1,]
#chisq test on progression#
progression.all <- sapply(progression.all, as.numeric) #chisquare test for the 4 groups
chisq.test(progression.all) # pvalue 3.252e-15, X-squared 80.198

progression.summ.red <- progression.summ %>%              # limit chisquare test to only extreme phenotypes
  filter(Gly.Cat.fw == 'NGR' | Gly.Cat.fw == 'T2D') %>% 
  pivot_wider(names_from=Gly.Cat.fw, values_from=n)

progression.summ.red <- as.data.frame(t(progression.summ.red))
names(progression.summ.red) <- progression.summ.red[1,]
progression.summ.red <- progression.summ.red[-1,]
#chisq test on progression#
progression.summ.red <- sapply(progression.summ.red, as.numeric) #chisquare test for the 2 groups
chisq.test(progression.summ.red) # pvalue 5.218e-11, X-squared 47.353

chord.progression <- data.frame('Baseline' = c(rep('Low.OGIS', 4), rep('Mid.OGIS', 4), rep('High.OGIS', 4)),
                                'Follow-up' = c(rep(c('IFG', 'IGT', 'NGR', 'T2D'),3)),
                                'n' = c(101, 55, 18, 46, 126, 34, 35, 24, 132, 12, 52, 6))

chord.progression$percentage <- c(chord.progression[chord.progression$Baseline=='Low.OGIS',]$n/220, # make in percentage for better analyses
                                  chord.progression[chord.progression$Baseline=='Mid.OGIS',]$n/219,
                                  chord.progression[chord.progression$Baseline=='High.OGIS',]$n/202)
chord.progression$Baseline <- factor(chord.progression$Baseline, levels = c('High.OGIS', 'Mid.OGIS', 'Low.OGIS'))
chord.progression$Follow.up <- factor(chord.progression$Follow.up, levels = c('NGR', 'IFG', 'IGT', 'T2D'))

alluvial::alluvial(chord.progression[,1:2], freq = chord.progression$n,
                   col = c('brown1','darkred', 'darkolivegreen4', 'darkorchid'), cex = 1.5, cex.axis = 1.5, blocks = T)
alluvial::alluvial(chord.progression[,1:2], freq = chord.progression$percentage,
                   col = c('brown1','darkred', 'darkolivegreen4', 'darkorchid'), cex = 1.5, cex.axis = 1.5, blocks = T)

chord.progression %>% 
  select(-c('percentage')) %>% 
  chordDiagram(grid.col = c(Low.OGIS= 'cadetblue2', Mid.OGIS = 'cornflowerblue', High.OGIS='blue',
                            IFG= 'brown1', IGT= 'darkred', 
                            NGR = 'darkolivegreen3', T2D = 'darkorchid'))

chord.progression %>% 
  select(-c('n')) %>% 
  chordDiagram(grid.col = c(Low.OGIS= 'cadetblue2', Mid.OGIS = 'cornflowerblue', High.OGIS='blue',
                            IFG= 'brown1', IGT= 'darkred', 
                            NGR = 'darkolivegreen3', T2D = 'darkorchid'))

progression.all <- progression.summ %>% 
  pivot_wider(names_from=OGIS.strat.fw, values_from=n) # prepare data for chisquare test
progression.all <- as.data.frame(t(progression.all))
names(progression.all) <- progression.all[1,]
progression.all <- progression.all[-1,]
#chisq test on progression#
progression.all <- sapply(progression.all, as.numeric) #chisquare test for the 4 groups
chisq.test(progression.all) # pvalue <2.2e-16, X-squared 229.1

## stability of OGIS stratifications ----
OGIS.strat.fw <- factor(ntile(biochem.mix[biochem.mix$timepoint=='Follow-up',]$OGIS.2h, 3))
OGIS.strat.fw <- plyr::revalue(OGIS.strat.fw, c('1'='Low.OGIS','2'='Mid.OGIS','3'='High.OGIS')) #stratify in tertiles 

ab <- biochem.bsl[biochem.bsl$samplerenamed%in%biochem.mix$sample.renamed,]
biochem.mix2 <- biochem.mix
biochem.mix2$OGIS.strat <- c(ab$OGIS.strat, OGIS.strat.fw)
biochem.mix2$OGIS.strat.fw <- factor(rep(biochem.mix2[biochem.mix2$timepoint=='Follow-up',]$OGIS.strat, 2))
biochem.mix2$OGIS.strat.fw  <- plyr::revalue(biochem.mix2$OGIS.strat.fw , c('1'='Low.OGIS','2'='Mid.OGIS','3'='High.OGIS')) #stratify in tertiles 
biochem.mix2$OGIS.strat  <- plyr::revalue(factor(biochem.mix2$OGIS.strat) , c('1'='Low.OGIS','2'='Mid.OGIS','3'='High.OGIS')) #stratify in tertiles 

progression.summ <- biochem.mix2 %>% 
  filter(timepoint == 'Baseline') %>% 
  select(c('OGIS.strat', 'OGIS.strat.fw')) %>% 
  group_by(OGIS.strat, OGIS.strat.fw) %>%
  tally()   

chord.progression <- data.frame('Baseline' = c(rep('Low.OGIS', 3), rep('Mid.OGIS', 3), rep('High.OGIS', 3)),
                                'Follow-up' = c(rep(c('Low.OGIS', 'Mid.OGIS', 'High.OGIS'),3)),
                                'n' = c(149, 48, 23, 55, 97, 67, 10, 69, 123))
chord.progression$percentage <- c(chord.progression[chord.progression$Baseline=='Low.OGIS',]$n/220, # make in percentage for better analyses
                                 chord.progression[chord.progression$Baseline=='Mid.OGIS',]$n/219,
                                 chord.progression[chord.progression$Baseline=='High.OGIS',]$n/202)
chord.progression$Baseline <- factor(chord.progression$Baseline, levels=c('Low.OGIS', 'Mid.OGIS', 'High.OGIS'))
chord.progression$Follow.up <- factor(chord.progression$Follow.up, levels=c('Low.OGIS', 'Mid.OGIS', 'High.OGIS'))

alluvial::alluvial(chord.progression[,1:2], freq = chord.progression$n,
                   col = c('cadetblue2','cornflowerblue', 'blue'), cex = 1.5, cex.axis = 1.5, blocks = T)
                  
alluvial::alluvial(chord.progression[,1:2], freq = chord.progression$percentage,
                   col = c('cadetblue2','cornflowerblue', 'blue'), cex = 1.5, cex.axis = 1.5, blocks = T)

## Biochemical progression at follow-up per OGIS strat and glycaemic status at follow-up ----
low.O.individuals <- biochem.mix %>% 
  filter(timepoint=='Baseline'&OGIS.strat=='Low.OGIS') %>% 
  pull(sample.renamed)
mid.O.individuals<- biochem.mix %>% 
  filter(timepoint=='Baseline'&OGIS.strat=='Mid.OGIS') %>% 
  pull(sample.renamed)
high.O.individuals<- biochem.mix %>% 
  filter(timepoint=='Baseline'&OGIS.strat=='High.OGIS') %>% 
  pull(sample.renamed)

biochem.mix.lowO <- biochem.mix[biochem.mix$sample.renamed%in%low.O.individuals,] # generate DF only for low OGIS 
biochem.mix.midO <- biochem.mix[biochem.mix$sample.renamed%in%mid.O.individuals,] # generate DF only for mid OGIS
biochem.mix.highO <- biochem.mix[biochem.mix$sample.renamed%in%high.O.individuals,] # generate DF only for high OGIS

# Low OGIS statistical divergences
biochem.mix.lowO.stats <- biochem.mix.lowO %>% 
  filter(timepoint=='Follow-up') %>% 
  group_by(Gly.Cat.fw) %>% 
  select(-c(names(biochem.mix.lowO)[1:11], 'OGIS.strat')) %>% 
  summarise_all(list(mean=mean, median=median, sd=sd), na.rm=T) %>% 
  pivot_longer(ends_with(c("_mean", "_sd", "_median")),
               names_to = c("name", ".value"),
               names_sep = "_") %>% 
  pivot_wider(names_from = Gly.Cat.fw,
              names_sep = "-",
              values_from = c(mean, sd, median))  # mean, sd and median at follow-up

biochem.mix.lowO <- biochem.mix.lowO[biochem.mix.lowO$timepoint=='Follow-up',] # p-value of stats only for follow-up
p.vals <- list()
for (i in seq(from=12, to=30)){
  model <- aov(biochem.mix.lowO[,i]~biochem.mix.lowO$Gly.Cat.fw, data = biochem.mix.lowO)
  anova.p <- summary(model)[[1]][['Pr(>F)']][1]
  tuks <- TukeyHSD(model)
  tuks.p <- tuks$`biochem.mix.lowO$Gly.Cat.fw`[,4]
  p.vals[[i]] <- c('anova'=anova.p, tuks.p)
}

biochem.mix.lowO.stats.p <- cbind(biochem.mix.lowO.stats, do.call(rbind.data.frame, p.vals)) # attach p.values to stats
names(biochem.mix.lowO.stats.p) <- c(names(biochem.mix.lowO.stats), 
                                    c('anova', 'IGT-IFG', 'NGR-IFG', 'T2D-IFG', 'NGR-IGT', 'T2D-IGT', 'T2D-NGR'))
write.table(biochem.mix.lowO.stats.p, 'InsulinSens/ResultsProgressionLowOGIS.txt', sep ='\t')

# Mid OGIS statistical divergences
biochem.mix.midO.stats <- biochem.mix.midO %>% 
  filter(timepoint=='Follow-up') %>% 
  group_by(Gly.Cat.fw) %>% 
  select(-c(names(biochem.mix.midO)[1:11], 'OGIS.strat')) %>% 
  summarise_all(list(mean=mean, median=median, sd=sd), na.rm=T) %>% 
  pivot_longer(ends_with(c("_mean", "_sd", "_median")),
               names_to = c("name", ".value"),
               names_sep = "_") %>% 
  pivot_wider(names_from = Gly.Cat.fw,
              names_sep = "-",
              values_from = c(mean, sd, median))  # mean, sd and median at follow-up

biochem.mix.midO <- biochem.mix.midO[biochem.mix.midO$timepoint=='Follow-up',] # p-value of stats only for follow-up
p.vals <- list()
for (i in seq(from=12, to=30)){
  model <- aov(biochem.mix.midO[,i]~biochem.mix.midO$Gly.Cat.fw, data = biochem.mix.midO)
  anova.p <- summary(model)[[1]][['Pr(>F)']][1]
  tuks <- TukeyHSD(model)
  tuks.p <- tuks$`biochem.mix.midO$Gly.Cat.fw`[,4]
  p.vals[[i]] <- c('anova'=anova.p, tuks.p)
}

biochem.mix.midO.stats.p <- cbind(biochem.mix.midO.stats, do.call(rbind.data.frame, p.vals)) # attach p.values to stats
names(biochem.mix.midO.stats.p) <- c(names(biochem.mix.midO.stats), 
                                     c('anova', 'IGT-IFG', 'NGR-IFG', 'T2D-IFG', 'NGR-IGT', 'T2D-IGT', 'T2D-NGR'))
write.table(biochem.mix.midO.stats.p, 'InsulinSens/ResultsProgressionMiddleOGIS.txt', sep ='\t')

# High OGIS statistical divergences
biochem.mix.highO.stats <- biochem.mix.highO %>% 
  filter(timepoint=='Follow-up') %>% 
  group_by(Gly.Cat.fw) %>% 
  select(-c(names(biochem.mix.highO)[1:11], 'OGIS.strat')) %>% 
  summarise_all(list(mean=mean, median=median, sd=sd), na.rm=T) %>% 
  pivot_longer(ends_with(c("_mean", "_sd", "_median")),
               names_to = c("name", ".value"),
               names_sep = "_") %>% 
  pivot_wider(names_from = Gly.Cat.fw,
              names_sep = "-",
              values_from = c(mean, sd, median))  # mean, sd and median at follow-up

biochem.mix.highO <- biochem.mix.highO[biochem.mix.highO$timepoint=='Follow-up',] # p-value of stats only for follow-up
p.vals <- list()
for (i in seq(from=12, to=30)){
  model <- aov(biochem.mix.highO[,i]~biochem.mix.highO$Gly.Cat.fw, data = biochem.mix.highO)
  anova.p <- summary(model)[[1]][['Pr(>F)']][1]
  tuks <- TukeyHSD(model)
  tuks.p <- tuks$`biochem.mix.highO$Gly.Cat.fw`[,4]
  p.vals[[i]] <- c('anova'=anova.p, tuks.p)
}

biochem.mix.highO.stats.p <- cbind(biochem.mix.highO.stats, do.call(rbind.data.frame, p.vals)) # attach p.values to stats
names(biochem.mix.highO.stats.p) <- c(names(biochem.mix.highO.stats), 
                                     c('anova', 'IGT-IFG', 'NGR-IFG', 'T2D-IFG', 'NGR-IGT', 'T2D-IGT', 'T2D-NGR'))
write.table(biochem.mix.highO.stats.p, 'InsulinSens/ResultsProgressionHighOGIS.txt', sep ='\t')

# boxplots
biochem.mix.lowO <- biochem.mix[biochem.mix$sample.renamed%in%low.O.individuals,] # re-generate DF only for low OGIS 
biochem.mix.midO <- biochem.mix[biochem.mix$sample.renamed%in%mid.O.individuals,] # re-generate DF only for mid OGIS
biochem.mix.highO <- biochem.mix[biochem.mix$sample.renamed%in%high.O.individuals,] # re-generate DF only for high OGIS

biochem.mix.lowO.plot <- cbind(biochem.mix.lowO[,10:31], 'timepoint'=biochem.mix.lowO$timepoint)
biochem.mix.lowO.plot$Gly.Cat.fw <- factor(biochem.mix.lowO.plot$Gly.Cat.fw, levels=c('NGR', 'IFG', 'IGT', 'T2D'))
biochem.mix.lowO.plot$progression <- paste0(biochem.mix.lowO.plot$timepoint, '_', c(rep('Low.OGIS', 220), as.character(biochem.mix.lowO.plot[biochem.mix.lowO.plot$timepoint=='Follow-up',]$Gly.Cat.fw)))
biochem.mix.lowO.plot$progression <- factor(biochem.mix.lowO.plot$progression, 
                                           levels = c('Baseline_Low.OGIS', 'Follow-up_NGR', 'Follow-up_IFG',
                                                      'Follow-up_IGT', 'Follow-up_T2D'))
biochem.mix.lowO.plot <- melt(biochem.mix.lowO.plot)
ggplot(biochem.mix.lowO.plot, aes(timepoint, value, group = interaction(progression))) + 
  geom_violin(alpha = .1, aes(fill=progression, color=progression)) +
  geom_quasirandom(alpha = .6, dodge.width = .9, aes(color = progression)) +
  geom_boxplot(width = .1, alpha = .7, outlier.colour = NA, position = position_dodge(.9), fill='white') +
  facet_wrap(~variable, scales = 'free') +
  scale_color_manual(values=c('cadetblue2', 'darkolivegreen3', 'brown1', 'darkred', 'darkorchid')) +
  scale_fill_manual(values=c('cadetblue2', 'darkolivegreen3', 'brown1', 'darkred', 'darkorchid')) +
  theme(panel.grid.major = element_line(colour = "gray85",linetype = "dashed"), 
        axis.text = element_text(size = 11, face = "bold", colour = "black"), 
        panel.background = element_rect(fill = NA), legend.position = c(0.9, 0.1), legend.title = element_text(face='bold'),
        legend.text = element_text(face='bold'), strip.background = element_rect(fill=NA),
        strip.text = element_text(face='bold', size =12), panel.border = element_rect(fill=NA)) +
  labs(x = NULL, y = NULL)

biochem.mix.midO.plot <- cbind(biochem.mix.midO[,10:31], 'timepoint'=biochem.mix.midO$timepoint)
biochem.mix.midO.plot$Gly.Cat.fw <- factor(biochem.mix.midO.plot$Gly.Cat.fw, levels=c('NGR', 'IFG', 'IGT', 'T2D'))
biochem.mix.midO.plot$progression <- paste0(biochem.mix.midO.plot$timepoint, '_', c(rep('Mid.OGIS', 219), as.character(biochem.mix.midO.plot[biochem.mix.midO.plot$timepoint=='Follow-up',]$Gly.Cat.fw)))
biochem.mix.midO.plot$progression <- factor(biochem.mix.midO.plot$progression, 
                                            levels = c('Baseline_Mid.OGIS', 'Follow-up_NGR', 'Follow-up_IFG',
                                                       'Follow-up_IGT', 'Follow-up_T2D'))
biochem.mix.midO.plot <- melt(biochem.mix.midO.plot)
ggplot(biochem.mix.midO.plot, aes(timepoint, value, group = interaction(progression))) + 
  geom_violin(alpha = .1, aes(fill=progression, color=progression)) +
  geom_quasirandom(alpha = .6, dodge.width = .9, aes(color = progression)) +
  geom_boxplot(width = .1, alpha = .7, outlier.colour = NA, position = position_dodge(.9), fill='white') +
  facet_wrap(~variable, scales = 'free') +
  scale_color_manual(values=c('cornflowerblue', 'darkolivegreen3', 'brown1', 'darkred', 'darkorchid')) +
  scale_fill_manual(values=c('cornflowerblue', 'darkolivegreen3', 'brown1', 'darkred', 'darkorchid')) +
  theme(panel.grid.major = element_line(colour = "gray85",linetype = "dashed"), 
        axis.text = element_text(size = 11, face = "bold", colour = "black"), 
        panel.background = element_rect(fill = NA), legend.position = c(0.9, 0.1), legend.title = element_text(face='bold'),
        legend.text = element_text(face='bold'), strip.background = element_rect(fill=NA),
        strip.text = element_text(face='bold', size =12), panel.border = element_rect(fill=NA)) +
  labs(x = NULL, y = NULL)

biochem.mix.highO.plot <- cbind(biochem.mix.highO[,10:31], 'timepoint'=biochem.mix.highO$timepoint)
biochem.mix.highO.plot$Gly.Cat.fw <- factor(biochem.mix.highO.plot$Gly.Cat.fw, levels=c('NGR', 'IFG', 'IGT', 'T2D'))
biochem.mix.highO.plot$progression <- paste0(biochem.mix.highO.plot$timepoint, '_', c(rep('high.OGIS', 202), as.character(biochem.mix.highO.plot[biochem.mix.highO.plot$timepoint=='Follow-up',]$Gly.Cat.fw)))
biochem.mix.highO.plot$progression <- factor(biochem.mix.highO.plot$progression, 
                                            levels = c('Baseline_high.OGIS', 'Follow-up_NGR', 'Follow-up_IFG',
                                                       'Follow-up_IGT', 'Follow-up_T2D'))
biochem.mix.highO.plot <- melt(biochem.mix.highO.plot)
ggplot(biochem.mix.highO.plot, aes(timepoint, value, group = interaction(progression))) + 
  geom_violin(alpha = .1, aes(fill=progression, color=progression)) +
  geom_quasirandom(alpha = .6, dodge.width = .9, aes(color = progression)) +
  geom_boxplot(width = .1, alpha = .7, outlier.colour = NA, position = position_dodge(.9), fill='white') +
  facet_wrap(~variable, scales = 'free') +
  scale_color_manual(values=c('blue', 'darkolivegreen3', 'brown1', 'darkred', 'darkorchid')) +
  scale_fill_manual(values=c('blue', 'darkolivegreen3', 'brown1', 'darkred', 'darkorchid')) +
  theme(panel.grid.major = element_line(colour = "gray85",linetype = "dashed"), 
        axis.text = element_text(size = 11, face = "bold", colour = "black"), 
        panel.background = element_rect(fill = NA), legend.position = c(0.9, 0.1), legend.title = element_text(face='bold'),
        legend.text = element_text(face='bold'), strip.background = element_rect(fill=NA),
        strip.text = element_text(face='bold', size =12), panel.border = element_rect(fill=NA)) +
  labs(x = NULL, y = NULL)


## Association between low OGIS/ high OGIS categorization and delta values (continuous manner) ----
delta.values <- readRDS('DeltaofLogValues.rds')
delta.values <- delta.values[delta.values$sample%in%c(low.O.individuals, high.O.individuals),]
biochem.mix.xtr <- biochem.mix[biochem.mix$OGIS.strat!='Mid.OGIS',]
delta.values$sample == biochem.mix.xtr[biochem.mix.xtr$timepoint=='Baseline',]$sample.renamed
delta.values$OGIS.strat <- biochem.mix.xtr[biochem.mix.xtr$timepoint=='Baseline',]$OGIS.strat 
delta.values$OGIS.strat <- droplevels(delta.values$OGIS.strat)

result.model <- list()
for (x in seq(from=6, to=23)){
  print(names(delta.values[x]))
  model <- data.frame('Variable' = delta.values[,x], delta.values[,1:5], 'OGIS' = delta.values$OGIS.strat)
  mod.a <- lm(Variable~OGIS+Gender+CenterID+Age, data = model, na.action = na.omit) # fixed effects: age sex and center
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['OGISHigh.OGIS', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['OGISHigh.OGIS', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$Variable, model$OGIS)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <- c(names(delta.values)[x],p.val, delta, delta.low, delta.up)
  result.model[[x]] <- result
}

delta.lm <- do.call(rbind.data.frame, result.model)
names(delta.lm) <- c('Variable', 'pval', 'delta', 'delta.low', 'delta.up')
delta.lm$fdr <- p.adjust(delta.lm$pval, method='fdr') #pvalue adjustment

delta.lm[,2:6] <- sapply(delta.lm[,2:6], as.numeric) 
delta.lm$Variable <- factor(delta.lm$Variable, levels=delta.lm$Variable)
ggplot(delta.lm, aes(delta, Variable, xmin=delta-abs(delta.low), xmax=delta+abs(delta.up))) + 
  geom_point(size =3) +
  geom_errorbarh(height=.2) +
  geom_vline(xintercept = 0, linetype = 'dashed') + xlim(-1,1) +
  annotation_custom(grobTree(textGrob('Low OGIS 2h', x = .1, y = .5, gp=gpar(col = 'gray', fontsize = 35, fontface = 'bold.italic', alpha =.6), rot = 90))) +
  annotation_custom(grobTree(textGrob('High OGIS 2h', x = .9, y = .5, gp=gpar(col = 'gray', fontsize = 35, fontface = 'bold.italic', alpha =.6), rot = 90))) +
  theme(panel.grid.major = element_line(colour = "gray85",linetype = "dashed"), 
        axis.text = element_text(size = 11, face = "bold", colour = "black"), 
        panel.background = element_rect(fill = NA), legend.position = c(0.9, 0.1), legend.title = element_text(face='bold'),
        legend.text = element_text(face='bold'), strip.background = element_rect(fill=NA),
        strip.text = element_text(face='bold', size =12), panel.border = element_rect(fill=NA),
        axis.title.x = element_text(face='bold', size = 15))+
  labs(x = "Cliff's Delta Estimate", y = NULL, title = 'Biochemical measures progression at 4 years follow-up vs extreme OGIS 2h stratifications')
