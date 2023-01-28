########################################################################################################################
##########################################  Supplementary: IFG vs IGT ##################################################
########################################################################################################################

## directory
setwd('/home/scratch/margar')

## libraries
library(tidyverse); library(tidyr); library(dplyr); library(reshape2); library(ggplot2); library(ggbeeswarm); library(ggpubr)
library(factoextra); library(phyloseq); library(effsize); library(cowplot); library(circlize)

## upload and preprocess data ----
biochem.bsl <- read.table('data/curatedmetadata/WP2.1_BL_MetadataWithInfl.txt', sep = '\t', header=T)

## differences in gene richness stratifications ----
tertiles.d <- biochem.bsl %>% 
  group_by(Gly.Cat2, DSGARLTertileGroup) %>% 
  tally() %>% 
  pivot_wider(names_from=DSGARLTertileGroup, values_from=n) ## count number of individuals in each tertile per each IGR type
tertiles.d <- as.data.frame(t(tertiles.d))
names(tertiles.d) <- tertiles.d[1,]
tertiles.d <- tertiles.d[-1,]
#chisq test on counts#
tertiles.d <- sapply(tertiles.d, as.numeric) #chisquare test for the 2 groups
chisq.test(tertiles.d) # pvalue 0.7538, non-related

## differences in enterotypes stratifications ----
entero.d <- biochem.bsl %>% 
  group_by(Gly.Cat2, enterotype) %>% 
  tally() %>% 
  pivot_wider(names_from=enterotype, values_from=n) ## count number of individuals in each tertile per each IGR type
entero.d <- as.data.frame(t(entero.d))
names(entero.d) <- entero.d[1,]
entero.d <- entero.d[-1,]
#chisq test on counts#
entero.d <- sapply(entero.d, as.numeric) #chisquare test for the 2 groups
chisq.test(entero.d) # pvalue 0.1587, non-related
entero.d <- as.data.frame(entero.d)
entero.d$entero <- c('Bact1', 'Bact2', 'Prev', 'Rum')
entero.d.m <- melt(entero.d)
ggplot(entero.d.m, aes(entero, value, fill=entero)) + geom_col() + facet_wrap(~variable, scales='free')

## gene richness (continous) and bacterial diversity ----
richness <- biochem.bsl %>% 
  select(c('Gly.Cat2', 'DSGeneAbundRichness', 'DSMGSAbundShannon'))

ggplot(richness, aes(Gly.Cat2, DSGeneAbundRichness)) + geom_boxplot() + stat_compare_means()
ggplot(richness, aes(Gly.Cat2, DSMGSAbundShannon)) + geom_boxplot() + stat_compare_means()

## compare biochemical variables between IFG and IGT -----
biochem.variables <- data.frame('Gly.Cat2'= biochem.bsl$Gly.Cat2, biochem.bsl[,16:57]) # extract only biochemical variables

biochem.stats <- biochem.variables %>% 
  group_by(Gly.Cat2) %>% 
  summarise_all(list(mean=mean, median=median, sd=sd), na.rm=T) %>% 
  pivot_longer(ends_with(c("_mean", "_sd", "_median")),
                    names_to = c("name", ".value"),
                    names_sep = "_") %>% 
  pivot_wider(names_from = Gly.Cat2,
              names_sep = "-",
              values_from = c(mean, sd, median)) ## compute mean, median and SD per IGR subtype


p.vals <- sapply(biochem.variables[,2:43], function(x){ #compute pvalues
  t.test(x~Gly.Cat2, data=biochem.variables)$p.val
})
p.vals  <- p.adjust(p.vals, method='fdr') #multitest adjust

biochem.stats$name == names(p.vals)
biochem.stats$FDR <- p.vals # generate final table

# write.table(biochem.stats, 'IFGvsIGT/ResultsBiochemPrediabeticsComparison.txt', sep ='\t')

biochem.variables.m <- melt(biochem.variables)
ggplot(biochem.variables.m, aes(Gly.Cat2, value), group = interaction(Gly.Cat2)) + 
  geom_violin(alpha = .1, aes(fill=Gly.Cat2, color=Gly.Cat2)) +
  geom_quasirandom(alpha = .6, dodge.width = .9, aes(color = Gly.Cat2)) +
  geom_boxplot(width = .1, alpha = .7, outlier.colour = NA, position = position_dodge(.9), fill='white')  +
  facet_wrap(~variable, scales = 'free') + 
  stat_compare_means(label = 'p.signif')

diff.biochem <- biochem.variables[,names(biochem.variables)%in%biochem.stats[biochem.stats$FDR<.1,]$name] ## plot only the significant ones
diff.biochem$Gly.Cat2 <- biochem.variables$Gly.Cat2
diff.biochem.m <- melt(diff.biochem)
ggplot(diff.biochem.m, aes(Gly.Cat2, value), group = interaction(Gly.Cat2)) + 
  geom_violin(alpha = .1, aes(fill=Gly.Cat2, color=Gly.Cat2)) +
  geom_quasirandom(alpha = .6, dodge.width = .9, aes(color = Gly.Cat2)) +
  geom_boxplot(width = .1, alpha = .7, outlier.colour = NA, position = position_dodge(.9), fill='white')  +
  facet_wrap(~variable, scales = 'free') + 
  stat_compare_means(label = 'p.signif') +
  scale_color_manual(values=c('brown1', 'darkred')) +
  scale_fill_manual(values=c('brown1', 'darkred')) +
  theme(legend.position = c(.9,.1))

biochem.variables.pca <- prcomp(na.omit(biochem.variables[,2:43]), center=T, scale. = T) #PCA 
fviz_screeplot(biochem.variables.pca)
row.names(get_pca_ind(biochem.variables.pca)$coord)
fviz_pca_biplot(biochem.variables.pca, geom.ind = 'point', col.ind = biochem.variables[row.names(biochem.variables)%in%row.names(get_pca_ind(biochem.variables.pca)$coord)
,]$Gly.Cat2, addEllipses = T)


## compare inflammation markers between groups ----
inflammatory <- data.frame('Gly.Cat2'= biochem.bsl$Gly.Cat2, biochem.bsl[,58:72]) ## extract only inflammatory markers
inflammatory.stats <- inflammatory %>% 
  group_by(Gly.Cat2) %>% 
  summarise_all(list(mean=mean, median=median, sd=sd), na.rm=T) %>% 
  pivot_longer(ends_with(c("_mean", "_sd", "_median")),
               names_to = c("name", ".value"),
               names_sep = "_") %>% 
  pivot_wider(names_from = Gly.Cat2,
              names_sep = "-",
              values_from = c(mean, sd, median)) ## compute mean, median and SD per IGR subtype


p.vals <- sapply(inflammatory[,2:16], function(x){ #compute pvalues
  t.test(x~Gly.Cat2, data=inflammatory)$p.val
})
p.vals  <- p.adjust(p.vals, method='fdr') #multitest adjust

inflammatory.stats$name == names(p.vals)
inflammatory.stats$FDR <- p.vals # generate final table
# write.table(inflammatory.stats, 'IFGvsIGT/ResultsInflammatoryMarkersPrediabeticsComparison.txt', sep ='\t')

inflammatory.m <- melt(inflammatory)
ggplot(inflammatory.m, aes(Gly.Cat2, value), group = interaction(Gly.Cat2)) + 
  geom_violin(alpha = .1, aes(fill=Gly.Cat2, color=Gly.Cat2)) +
  geom_quasirandom(alpha = .6, dodge.width = .9, aes(color = Gly.Cat2)) +
  geom_boxplot(width = .1, alpha = .7, outlier.colour = NA, position = position_dodge(.9), fill='white')  +
  facet_wrap(~variable, scales = 'free') + 
  stat_compare_means(label = 'p.signif') +
  scale_color_manual(values=c('brown1', 'darkred')) +
  scale_fill_manual(values=c('brown1', 'darkred')) +
  theme(legend.position = c(.9,.1))

## Metagenome differences at baseline between IGT and IFG ----
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
pcoa.df <- data.frame(pcoa.jsd.bsl$vectors, 'Gly.Cat2' = biochem.bsl$Gly.Cat2)
ggplot(pcoa.df, aes(Axis.1, Axis.2, color=Gly.Cat2)) + geom_point() +
  geom_hline(yintercept = 0, linetype = 'dashed') + geom_vline(xintercept=0, linetype = 'dashed') +
  stat_ellipse() + scale_color_manual(values=c('brown1', 'darkred')) +
  theme(panel.grid.major = element_line(colour = "gray85", linetype = "dashed"), axis.title = element_text(size = 17, face = "bold"), 
        axis.text = element_text(size = 13,  face = "bold"), legend.text = element_text(size = 12, face = "bold"), 
        legend.title = element_text(face = "bold"), panel.background = element_rect(fill = NA), 
        legend.key = element_rect(fill = NA), legend.background = element_rect(fill = NA), panel.border = element_rect(fill=NA)) +
  labs(x = "PCoA 1 (expl. variance 13.4%)", y = "PCoA 2 (expl. variance 10.2%)", colour = "Prediabetic type")
vegan::adonis(jsd.dist.mgs~Gly.Cat2, data=biochem.bsl)

byc.dist.mgs <- distance(physeq.bsl, method = 'bray') ## Bray curtis
pcoa.byc.bsl <- ordinate(physeq.bsl, 'PCoA', distance = byc.dist.mgs)
pcoa.df <- data.frame(pcoa.byc.bsl$vectors, 'Gly.Cat2' = biochem.bsl$Gly.Cat2)
ggplot(pcoa.df, aes(Axis.1, Axis.2, color=Gly.Cat2)) + geom_point() +
  geom_hline(yintercept = 0, linetype = 'dashed') + geom_vline(xintercept=0, linetype = 'dashed') +
  stat_ellipse() + scale_color_manual(values=c('brown1', 'darkred')) +
  theme(panel.grid.major = element_line(colour = "gray85", linetype = "dashed"), axis.title = element_text(size = 17, face = "bold"), 
        axis.text = element_text(size = 13,  face = "bold"), legend.text = element_text(size = 12, face = "bold"), 
        legend.title = element_text(face = "bold"), panel.background = element_rect(fill = NA), 
        legend.key = element_rect(fill = NA), legend.background = element_rect(fill = NA), panel.border = element_rect(fill=NA)) +
  labs(x = "PCoA 1 (expl. variance 9.62%)", y = "PCoA 2 (expl. variance 6.78%)", colour = "Prediabetic type")
vegan::adonis(byc.dist.mgs~Gly.Cat2, data=biochem.bsl)

# MGS level
biochem.bsl$samplerenamed == row.names(qmp.bsl)
mgs.model <- data.frame('StudyID' = biochem.bsl$StudyID, 'center' = biochem.bsl$CenterID, # generate DF for the modeling step
                        'gender' = biochem.bsl$Gender, 'age'=biochem.bsl$Age, 
                        'Gly.Cat2' = biochem.bsl$Gly.Cat2, qmp.bsl)

mgs.model.igt <- mgs.model[mgs.model$Gly.Cat2=='IGT',]
mgs.model.ifg <- mgs.model[mgs.model$Gly.Cat2=='IFG',]
mgs.model <- as.data.frame(rbind(mgs.model.igt, sample_n(mgs.model.ifg, nrow(mgs.model.igt))))

result.model <- list()
for (x in seq(from=6, to=726)){
  print(names(mgs.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <- data.frame('mgs' = mgs.model[,x], mgs.model[,1:5])
  mod.a <- lm(mgs~Gly.Cat2+gender+center+age, data = model, na.action = na.omit) # fixed effects: age sex and center
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Gly.Cat2IGT', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Gly.Cat2IGT', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Gly.Cat2)
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
                        'Gly.Cat2' = biochem.bsl$Gly.Cat2, genera.counts)

genus.model.igt <- genus.model[genus.model$Gly.Cat2=='IGT',]
genus.model.ifg <- genus.model[genus.model$Gly.Cat2=='IFG',]
genus.model <- as.data.frame(rbind(genus.model.igt, sample_n(genus.model.ifg, nrow(genus.model.igt))))

result.model <- list()
for (x in seq(from=6, to=128)){
  print(names(genus.model[x]))
  model <- data.frame('genus' = genus.model[,x], genus.model[,1:5])
  mod.a <- lm(genus~Gly.Cat2+gender+center+age, data = model, na.action = na.omit) # fixed effects: age sex and center
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Gly.Cat2IGT', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Gly.Cat2IGT', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$genus, model$Gly.Cat2)
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

## GMMs differences between IGT and IFG ----
GMMs <- read.table('data/WP2.1_775_IGR_GMM-QMP_igrHA_Ent_M.csv', header = T, row.names = 1)
GMMs <- GMMs[,7:103]
biochem.bsl$samplerenamed == row.names(GMMs)
gmm.model <- data.frame('StudyID' = biochem.bsl$StudyID, 'center' = biochem.bsl$CenterID, # generate DF for the modeling step
                        'gender' = biochem.bsl$Gender, 'age'=biochem.bsl$Age, 
                        'Gly.Cat2' = biochem.bsl$Gly.Cat2, GMMs)
gmm.model.igt <- gmm.model[gmm.model$Gly.Cat2=='IGT',]
gmm.model.ifg <- gmm.model[gmm.model$Gly.Cat2=='IFG',]
gmm.model <- as.data.frame(rbind(gmm.model.igt, sample_n(gmm.model.ifg, nrow(gmm.model.igt))))

result.model <- list()
for (x in seq(from=6, to=102)){
  print(names(gmm.model[x]))
  model <- data.frame('GMM' = gmm.model[,x], gmm.model[,1:5])
  mod.a <- lm(GMM~Gly.Cat2+gender+center+age, data = model, na.action = na.omit) # fixed effects: age sex and center
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Gly.Cat2IGT', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Gly.Cat2IGT', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$GMM, model$Gly.Cat2)
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
            'IFGvsIGT/ReducedDatasetMGSGenusGMM_linearmodelsIFGvsIGT.txt', sep = '\t')

## structural variants ----
# data upload
deletions <- read.table('SV-followup/dsgv_baseline.csv', header = T, sep = ',', row.names = 1)
row.names(deletions) == biochem.bsl$SampleID
deletions <- deletions[match(biochem.bsl$SampleID, row.names(deletions)),] # match order between DFs

variations <- read.table('SV-followup/vsgv_baseline.csv', header = T, sep = ',', row.names = 1)
row.names(deletions) == row.names(variations)
variations <- variations[match(row.names(deletions), row.names(variations)),] # match order between DFs

SVs.per.sample <- data.frame('deletions' = rowSums(!is.na(deletions)&deletions==0), 'variations' = rowSums(!is.na(variations))) # total number of deletions and variable regions per sample
SVs.per.sample$Gly.Cat2 <- biochem.bsl$Gly.Cat2

p1<-ggplot(SVs.per.sample, aes(deletions, variations, color=Gly.Cat2)) + geom_point() + geom_smooth(method='lm') + # correlations between deletions and variable regions per sample
  scale_color_manual(values = c('brown1','darkred')) + 
  scale_fill_manual(values = c('brown1','darkred')) +
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

p2 <- ggplot(SVs.per.sample, aes(Gly.Cat2, deletions, color=Gly.Cat2)) +  # boxplots for the deletions
  geom_boxplot(outlier.colour = NA, lwd = 1.4) +
  geom_quasirandom(alpha=.4) + 
  stat_compare_means()+
  scale_color_manual(values = c('brown1','darkred')) + 
  scale_fill_manual(values = c('brown1','darkred')) + 
  theme(panel.grid.major = element_line(colour = "gray85", linetype = "dashed"), 
        axis.title = element_text(size = 18, face = "bold"), 
        axis.text = element_text(size = 15, face = "bold"), 
        panel.background = element_rect(fill = NA), legend.position = "none", 
        panel.border = element_rect(fill=NA)) +
  labs(x = NULL, y = "Number of deletions", colour = NULL) 

p3 <- ggplot(SVs.per.sample, aes(Gly.Cat2, variations, color=Gly.Cat2)) + # boxplots for the variable regions
  geom_boxplot(outlier.colour = NA, lwd = 1.4) +
  geom_quasirandom(alpha=.4) + 
  stat_compare_means()+
  scale_color_manual(values = c('brown1','darkred')) + 
  scale_fill_manual(values = c('brown1','darkred')) + 
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
                        'Gly.Cat2' = biochem.bsl$Gly.Cat2, deletions)

result.model <- list()
for (x in seq(from=6, to=7720)){
  print(names(deletions.model[x]))
  model <- data.frame('Deletion' = deletions.model[,x], deletions.model[,1:5])
  mod.a <- lm(Deletion~Gly.Cat2+gender+center+age, data = model, na.action = na.omit) # fixed effects: age sex and center
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Gly.Cat2IGT', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Gly.Cat2IGT', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$Deletion, model$Gly.Cat2)
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
                              'Gly.Cat2' = biochem.bsl$Gly.Cat2, variations)

result.model <- list()
for (x in seq(from=6, to=3104)){
  print(names(variations.model[x]))
  model <- data.frame('Variation' = variations.model[,x], variations.model[,1:5])
  mod.a <- lm(Variation~Gly.Cat2+gender+center+age, data = model, na.action = na.omit) # fixed effects: age sex and center
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Gly.Cat2IGT', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Gly.Cat2IGT', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$Variation, model$Gly.Cat2)
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
  
## Progression and regression differences in IFG/IGT individuals ----
biochem.mix <- read.table('data/curatedmetadata/WP2.1_BlFw-641_MetadataImputed_OnlyBIOCHEM.txt', header = T) # read the mix of baseline and followup biochemical data
biochem.mix$Gly.Cat.fw <- c(biochem.mix[biochem.mix$timepoint=='Follow-up', 'Gly.Cat.2'],
                                biochem.mix[biochem.mix$timepoint=='Follow-up', 'Gly.Cat.2']) # add column with only the glycaemic status at followup
progression.summ <- biochem.mix %>% 
  filter(timepoint == 'Baseline') %>% 
  select(c('Gly.Cat.2', 'Gly.Cat.fw')) %>% 
  group_by(Gly.Cat.2, Gly.Cat.fw) %>%
  tally()                                     # count the number of individuals progressing to each status

progression.all <- progression.summ %>% 
  pivot_wider(names_from=Gly.Cat.fw, values_from=n) # prepare data for chisquare test
progression.all <- as.data.frame(t(progression.all))
names(progression.all) <- progression.all[1,]
progression.all <- progression.all[-1,]
#chisq test on progression#
progression.all <- sapply(progression.all, as.numeric) #chisquare test for the 4 groups
chisq.test(progression.all) # pvalue 3.928e-14, X-squared 65.496

progression.summ.red <- progression.summ %>%              # limit chisquare test to only extreme phenotypes
  filter(Gly.Cat.fw == 'NGR' | Gly.Cat.fw == 'T2D') %>% 
  pivot_wider(names_from=Gly.Cat.fw, values_from=n)

progression.summ.red <- as.data.frame(t(progression.summ.red))
names(progression.summ.red) <- progression.summ.red[1,]
progression.summ.red <- progression.summ.red[-1,]
#chisq test on progression#
progression.summ.red <- sapply(progression.summ.red, as.numeric) #chisquare test for the 2 groups
chisq.test(progression.summ.red) # pvalue 0.01552, X-squared 5.8563

chord.progression <- data.frame('Baseline' = c(rep('IFG', 4), rep('IGT', 4)),
                                'Follow-up' = c(rep(c('IFG', 'IGT', 'NGR', 'T2D'),2)),
                                'n' = c(327, 60, 90, 53, 32, 41, 15, 23))

chord.progression$percentage <- c(chord.progression[chord.progression$Baseline=='IFG',]$n/530, # make in percentage for better analyses
                                  chord.progression[chord.progression$Baseline=='IGT',]$n/111)
alluvial::alluvial(chord.progression[,1:2], freq = chord.progression$n,
                   col = c('brown1','darkred', 'darkolivegreen4', 'darkorchid'), cex = 1.5, cex.axis = 1.5, blocks = T)
alluvial::alluvial(chord.progression[,1:2], freq = chord.progression$percentage,
                   col = c('brown1','darkred', 'darkolivegreen4', 'darkorchid'), cex = 1.5, cex.axis = 1.5, blocks = T)

chord.progression$Baseline <- paste0(chord.progression$Baseline, '_BSL')
chord.progression$Follow.up <- paste0(chord.progression$Follow.up, '_FW')
chord.progression %>% 
  select(-c('percentage')) %>% 
  chordDiagram(grid.col = c(IFG_BSL= 'brown1', IGT_BSL = 'darkred', 
                            IFG_FW= 'brown1', IGT_FW = 'darkred', 
                            NGR_FW = 'darkolivegreen3', T2D_FW = 'darkorchid'))

chord.progression %>% 
  select(-c('n')) %>% 
  chordDiagram(grid.col = c(IFG_BSL= 'brown1', IGT_BSL = 'darkred', 
                            IFG_FW= 'brown1', IGT_FW = 'darkred', 
                            NGR_FW = 'darkolivegreen3', T2D_FW = 'darkorchid'))

## Biochemical progression at follow-up per prediabetic group and glycaemic status at follow-up ----
ifg.individuals <- biochem.mix %>% 
  filter(timepoint=='Baseline'&Gly.Cat.2=='IFG') %>% 
  pull(sample.renamed)
igt.individuals<- biochem.mix %>% 
  filter(timepoint=='Baseline'&Gly.Cat.2=='IGT') %>% 
  pull(sample.renamed)
  
biochem.mix.ifg <- biochem.mix[biochem.mix$sample.renamed%in%ifg.individuals,] # generate DF only for IFG 
biochem.mix.igt <- biochem.mix[biochem.mix$sample.renamed%in%igt.individuals,] # generate DF only for IGT 

# IFG statistical divergences
biochem.mix.ifg.stats <- biochem.mix.ifg %>% 
  filter(timepoint=='Follow-up') %>% 
  group_by(Gly.Cat.fw) %>% 
  select(-names(biochem.mix.ifg)[1:11]) %>% 
  summarise_all(list(mean=mean, median=median, sd=sd), na.rm=T) %>% 
  pivot_longer(ends_with(c("_mean", "_sd", "_median")),
               names_to = c("name", ".value"),
               names_sep = "_") %>% 
  pivot_wider(names_from = Gly.Cat.fw,
              names_sep = "-",
              values_from = c(mean, sd, median))  # mean, sd and median at follow-up

biochem.mix.ifg <- biochem.mix.ifg[biochem.mix.ifg$timepoint=='Follow-up',] # p-value of stats only for follow-up
p.vals <- list()
for (i in seq(from=12, to=30)){
  model <- aov(biochem.mix.ifg[,i]~biochem.mix.ifg$Gly.Cat.fw, data = biochem.mix.ifg)
  anova.p <- summary(model)[[1]][['Pr(>F)']][1]
  tuks <- TukeyHSD(model)
  tuks.p <- tuks$`biochem.mix.ifg$Gly.Cat.fw`[,4]
  p.vals[[i]] <- c('anova'=anova.p, tuks.p)
}

biochem.mix.ifg.stats.p <- cbind(biochem.mix.ifg.stats, do.call(rbind.data.frame, p.vals)) # attach p.values to stats
names(biochem.mix.ifg.stats.p) <- c(names(biochem.mix.ifg.stats), 
                                    c('anova', 'IGT-IFG', 'NGR-IFG', 'T2D-IFG', 'NGR-IGT', 'T2D-IGT', 'T2D-NGR'))
write.table(biochem.mix.ifg.stats.p, 'IFGvsIGT/ResultsProgressionIFG.txt', sep ='\t')

# IGT statistical divergences
biochem.mix.igt.stats <- biochem.mix.igt %>% 
  filter(timepoint=='Follow-up') %>% 
  group_by(Gly.Cat.fw) %>% 
  select(-names(biochem.mix.igt)[1:11]) %>% 
  summarise_all(list(mean=mean, median=median, sd=sd), na.rm=T) %>% 
  pivot_longer(ends_with(c("_mean", "_sd", "_median")),
               names_to = c("name", ".value"),
               names_sep = "_") %>% 
  pivot_wider(names_from = Gly.Cat.fw,
              names_sep = "-",
              values_from = c(mean, sd, median))  # mean, sd and median at follow-up

biochem.mix.igt <- biochem.mix.igt[biochem.mix.igt$timepoint=='Follow-up',] # p-value of stats only for follow-up
p.vals <- list()
for (i in seq(from=12, to=30)){
  model <- aov(biochem.mix.igt[,i]~biochem.mix.igt$Gly.Cat.fw, data = biochem.mix.igt)
  anova.p <- summary(model)[[1]][['Pr(>F)']][1]
  tuks <- TukeyHSD(model)
  tuks.p <- tuks$`biochem.mix.igt$Gly.Cat.fw`[,4]
  p.vals[[i]] <- c('anova'=anova.p, tuks.p)
}

biochem.mix.igt.stats.p <- cbind(biochem.mix.igt.stats, do.call(rbind.data.frame, p.vals)) # attach p.values to stats
names(biochem.mix.igt.stats.p) <- c(names(biochem.mix.igt.stats), 
                                    c('anova', 'IGT-IFG', 'NGR-IFG', 'T2D-IFG', 'NGR-IGT', 'T2D-IGT', 'T2D-NGR'))
write.table(biochem.mix.igt.stats.p, 'IFGvsIGT/ResultsProgressionIGT.txt', sep ='\t')

# boxplots
biochem.mix.ifg <- biochem.mix[biochem.mix$sample.renamed%in%ifg.individuals,] # re-generate DF only for IFG 
biochem.mix.igt <- biochem.mix[biochem.mix$sample.renamed%in%igt.individuals,] # re-generate DF only for IGT 

biochem.mix.ifg.plot <- cbind(biochem.mix.ifg[,10:31], 'timepoint'=biochem.mix.ifg$timepoint)
biochem.mix.ifg.plot$Gly.Cat.2 <- factor(biochem.mix.ifg.plot$Gly.Cat.2, levels=c('NGR', 'IFG', 'IGT', 'T2D'))
biochem.mix.ifg.plot$progression <- paste0(biochem.mix.ifg.plot$timepoint, '_', biochem.mix.ifg.plot$Gly.Cat.2)
biochem.mix.ifg.plot$progression <- factor(biochem.mix.ifg.plot$progression, 
                                           levels = c('Baseline_IFG', 'Follow-up_NGR', 'Follow-up_IFG',
                                                      'Follow-up_IGT', 'Follow-up_T2D'))
biochem.mix.ifg.plot <- melt(biochem.mix.ifg.plot)
ggplot(biochem.mix.ifg.plot, aes(timepoint, value, group = interaction(progression))) + 
  geom_violin(alpha = .1, aes(fill=progression, color=progression)) +
  geom_quasirandom(alpha = .6, dodge.width = .9, aes(color = progression)) +
  geom_boxplot(width = .1, alpha = .7, outlier.colour = NA, position = position_dodge(.9), fill='white') +
  facet_wrap(~variable, scales = 'free') +
  scale_color_manual(values=c('brown1', 'darkolivegreen3', 'brown1', 'darkred', 'darkorchid')) +
  scale_fill_manual(values=c('brown1', 'darkolivegreen3', 'brown1', 'darkred', 'darkorchid')) +
  theme(panel.grid.major = element_line(colour = "gray85",linetype = "dashed"), 
        axis.text = element_text(size = 11, face = "bold", colour = "black"), 
        panel.background = element_rect(fill = NA), legend.position = c(0.9, 0.1), legend.title = element_text(face='bold'),
        legend.text = element_text(face='bold'), strip.background = element_rect(fill=NA),
        strip.text = element_text(face='bold', size =12), panel.border = element_rect(fill=NA)) +
  labs(x = NULL, y = NULL)

biochem.mix.igt.plot <- cbind(biochem.mix.igt[,10:31], 'timepoint'=biochem.mix.igt$timepoint)
biochem.mix.igt.plot$Gly.Cat.2 <- factor(biochem.mix.igt.plot$Gly.Cat.2, levels=c('NGR', 'IFG', 'IGT', 'T2D'))
biochem.mix.igt.plot$progression <- paste0(biochem.mix.igt.plot$timepoint, '_', biochem.mix.igt.plot$Gly.Cat.2)
biochem.mix.igt.plot$progression <- factor(biochem.mix.igt.plot$progression, 
                                           levels = c('Baseline_IGT', 'Follow-up_NGR', 'Follow-up_IFG',
                                                      'Follow-up_IGT', 'Follow-up_T2D'))
biochem.mix.igt.plot <- melt(biochem.mix.igt.plot)
ggplot(biochem.mix.igt.plot, aes(timepoint, value, group = interaction(progression))) + 
  geom_violin(alpha = .1, aes(fill=progression, color=progression)) +
  geom_quasirandom(alpha = .6, dodge.width = .9, aes(color = progression)) +
  geom_boxplot(width = .1, alpha = .7, outlier.colour = NA, position = position_dodge(.9), fill='white') +
  facet_wrap(~variable, scales = 'free') +
  scale_color_manual(values=c('darkred', 'darkolivegreen3', 'brown1', 'darkred', 'darkorchid')) +
  scale_fill_manual(values=c('darkred', 'darkolivegreen3', 'brown1', 'darkred', 'darkorchid')) +
  theme(panel.grid.major = element_line(colour = "gray85",linetype = "dashed"), 
        axis.text = element_text(size = 11, face = "bold", colour = "black"), 
        panel.background = element_rect(fill = NA), legend.position = c(0.9, 0.1), legend.title = element_text(face='bold'),
        legend.text = element_text(face='bold'), strip.background = element_rect(fill=NA),
        strip.text = element_text(face='bold', size =12), panel.border = element_rect(fill=NA)) +
  labs(x = NULL, y = NULL)

## Association between IFG/IGT categorization and delta values (continuous manner) ----
delta.values <- readRDS('DeltaofLogValues.rds')
delta.values$sample == biochem.mix[biochem.mix$timepoint=='Baseline',]$sample.renamed
delta.values$Prediab.type <- biochem.mix[biochem.mix$timepoint=='Baseline',]$Gly.Cat.2 

result.model <- list()
for (x in seq(from=6, to=23)){
  print(names(delta.values[x]))
  model <- data.frame('Variable' = delta.values[,x], delta.values[,1:5], 'prediab' = delta.values$Prediab.type)
  mod.a <- lm(Variable~prediab+Gender+CenterID+Age, data = model, na.action = na.omit) # fixed effects: age sex and center
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['prediabIGT', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['prediabIGT', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$Variable, model$prediab)
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
  annotation_custom(grobTree(textGrob('Impaired Fasting Glucose (IFG)', x = .1, y = .5, gp=gpar(col = 'gray', fontsize = 35, fontface = 'bold.italic', alpha =.6), rot = 90))) +
  annotation_custom(grobTree(textGrob('Impaired Glucose Tolerance (IGT)', x = .9, y = .5, gp=gpar(col = 'gray', fontsize = 35, fontface = 'bold.italic', alpha =.6), rot = 90))) +
  theme(panel.grid.major = element_line(colour = "gray85",linetype = "dashed"), 
        axis.text = element_text(size = 11, face = "bold", colour = "black"), 
        panel.background = element_rect(fill = NA), legend.position = c(0.9, 0.1), legend.title = element_text(face='bold'),
        legend.text = element_text(face='bold'), strip.background = element_rect(fill=NA),
        strip.text = element_text(face='bold', size =12), panel.border = element_rect(fill=NA),
        axis.title.x = element_text(face='bold', size = 15))+
  labs(x = "Cliff's Delta Estimate", y = NULL, title = 'Biochemical measures progression at 4 years follow-up vs prediabetic classification')