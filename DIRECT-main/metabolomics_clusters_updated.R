#####################################################################################################
################################## extended metabolomics ############################################
#####################################################################################################

setwd('/home/scratch/margar')

## libraries 
library (xlsx) ### Saving to spreadsheet
library (data.table) ### Fast read of large files into R
library (WGCNA) ### -	Clustering software. Previously reported work done using v1.34
library (flashClust) ### Clustering software
library (ppcor) ### Partial Spearman correlations, for confounder analysis. Previously reported work done using v1.0
library (gplots) ### Plotting
library (cowplot) ### Plotting; to arrange several plots on the same page
library (ggplot2) ### Plotting
library (plyr) ### Data transformations
library (tidyverse)

## settings for WGCNA
cor_method          = "spearman" ### for association with clinical parameters
corFun_tmp          = "bicor"
cluster_method      = "average"
corOptions_list     = list (use = 'pairwise.complete.obs') 
corOptions_str      = "use = 'pairwise.complete.obs'"
BH_pval_asso_cutoff = 0.05
NetworkType         = "signed" ### Signed-network (as the PC1 and median profile does not make sense as a summary measure of a cluster with anticorrelated metabolites.)

## load data
mets <- readRDS('data/metabolomics/MetabolitesResiduals.rds')
mtdt <- read.table('data/curatedmetadata/WP2.1_BL_MetadataWithInfl2.txt', header =T, row.names=1, sep='\t')
annot <- read.table('data/metabolomics/Chemical_Annotation_Metabolone_DTU_sh.txt', header=T, row.names=1)
met.clusters <- readRDS('data/metabolomics/metabolitemodules.rds')
cluster_mapping_file <- read.table('data/metabolomics/moduleannotationscurated2.txt', header=T, sep='\t')
module.content <- read.table('data/metabolomics/modulemappingfile.txt') ## annotate clusters manually

mets.to.be.added <- c(row.names(module.content[module.content$module=='grey',]),
                      row.names(module.content[module.content$module=='turquoise',]))

identical(row.names(mets),row.names(met.clusters))
met.clusters2 <- cbind(met.clusters, mets[,names(mets)%in%mets.to.be.added])
met.clusters2$M38 <- NULL
met.clusters2$M107 <- NULL
# saveRDS(met.clusters2, 'data/metabolomics/mergedMetClusters.rds')

## repeat analyses ----
annotation.pheat <- data.frame('GR'=mtdt$DSGARLTertileGroup, 'entero' = mtdt$enterotype); row.names(annotation.pheat) <- row.names(mtdt)
pheatmap::pheatmap(met.clusters2, 
                   annotation_row = annotation.pheat, scale='row')#, cutree_cols = c(2), cutree_rows = c(2))

met.clusters2.micro <- met.clusters2
met.clusters2.micro$entero <- plyr::revalue(mtdt$enterotype, c('Bact1'='noBact2', 'Prev'='noBact2', 'Rum'= 'noBact2'))
met.clusters2.micro.ent <- sapply(met.clusters2.micro[,1:369], function(x)(
  t.test(x~met.clusters2.micro$entero)$p.val
))
met.clusters2.micro.ent <- p.adjust(met.clusters2.micro.ent, method = 'fdr')
met.clusters.ent <- met.clusters2.micro %>% 
  group_by(entero) %>% 
  summarise(across(where(is.numeric), mean)) %>% 
  t() %>% 
  as.data.frame()
names(met.clusters.ent) <- met.clusters.ent[1,] ; met.clusters.ent <- met.clusters.ent[-1,]
row.names(met.clusters.ent) == names(met.clusters2.micro.ent)
met.clusters.ent$fdr <- met.clusters2.micro.ent

annotation.pheat <- data.frame('GR'=mtdt$DSGARLTertileGroup, 'entero' = met.clusters2.micro$entero); row.names(annotation.pheat) <- row.names(mtdt)
pheatmap::pheatmap(met.clusters2.micro[,names(met.clusters2.micro)%in%names(met.clusters2.micro.ent[met.clusters2.micro.ent<.1])], 
                   annotation_row = annotation.pheat, scale='row')

met.clusters2.micro$GR <- factor(mtdt$DSGARLTertileGroup)
met.clusters2.micro.gr <- sapply(met.clusters2.micro[met.clusters2.micro$GR!='HGC',1:369], function(x)(
  t.test(x~droplevels(met.clusters2.micro[met.clusters2.micro$GR!="HGC",]$GR))$p.val
))
met.clusters2.micro.gr <- p.adjust(met.clusters2.micro.gr, method = 'fdr')
met.clusters.gr <- met.clusters2.micro %>% 
  group_by(GR) %>% 
  summarise(across(where(is.numeric), mean)) %>% 
  t() %>% 
  as.data.frame()
names(met.clusters.gr) <- met.clusters.gr[1,] ; met.clusters.gr <- met.clusters.gr[-1,]
row.names(met.clusters.gr) == names(met.clusters2.micro.ent)
met.clusters.gr$fdr <- met.clusters2.micro.gr
met.clusters.gr$module.name <- sapply(row.names(met.clusters.gr), function(a){
  paste0(a, ': ', cluster_mapping_file[which(cluster_mapping_file$New_Name==a),]$Description)
})


write.table(cbind(met.clusters.ent, met.clusters.gr), '/home/margar/download/export2_MEtClustersMean_GRandENtero.txt', sep='\t')
pheatmap::pheatmap(met.clusters2.micro[met.clusters2.micro$GR!='HGC',names(met.clusters2.micro)%in%names(met.clusters2.micro.ent[met.clusters2.micro.gr<.1])], 
                   annotation_row = annotation.pheat, scale='row')

keep.modules.micro <- unique(names(met.clusters2.micro.gr[met.clusters2.micro.gr<.1]), names(met.clusters2.micro.ent[met.clusters2.micro.ent<.1]))
pheatmap::pheatmap(met.clusters2[, names(met.clusters2)%in%keep.modules.micro], 
                   annotation_row = annotation.pheat, scale='column')

## do clusters associate to progression'? ----
mtdt.all <- read.table('data/curatedmetadata/WP2.1_BlFw-641_MetadataImputed_OnlyBIOCHEM_2.txt', header=T, sep = '\t')
mtdt.all[mtdt.all$timepoint=='Baseline',]$StudyID == mtdt.all[mtdt.all$timepoint!='Baseline',]$StudyID 
mtdt.all$progression <- rep(mtdt.all[mtdt.all$timepoint!='Baseline',]$Gly.Cat, 2)
mtdt.all <- mtdt.all[mtdt.all$timepoint=='Baseline',]
met.clusters2.prog <- met.clusters2[row.names(met.clusters2)%in%mtdt.all$StudyID,]
identical(row.names(met.clusters2.prog), mtdt.all$StudyID)
progression.annotation <- data.frame(mtdt.all$progression); row.names(progression.annotation) <- mtdt.all$StudyID
pheatmap::pheatmap(met.clusters2.prog, annotation_row = progression.annotation, scale = 'column')

## progression differential metabolites ----
met.clusters2.prog$t2d <- plyr::revalue(mtdt.all$progression, c('IGR'='noT2D', 'NGR'='noT2D'))
met.clusters2.prog$ngr <- plyr::revalue(mtdt.all$progression, c('IGR'='noNGR', 'T2D'='noNGR'))

# sapply(met.clusters2.prog[,1:369], function(x)(
#   print(hist(x))
# ))

met.clusters2.prog.t2d <- sapply(met.clusters2.prog[,1:369], function(x)(
  t.test(x~met.clusters2.prog$t2d)$p.val
))
met.clusters2.prog.t2d <- p.adjust(met.clusters2.prog.t2d, method = 'fdr')

met.clusters2.prog.ngr <- sapply(met.clusters2.prog[,1:369], function(x)(
  t.test(x~met.clusters2.prog$ngr)$p.val
))
met.clusters2.prog.ngr <- p.adjust(met.clusters2.prog.ngr, method = 'fdr')

library(tibble)
modules.stats <- data.frame(met.clusters2.prog %>% 
                              group_by(t2d) %>% 
                              dplyr::summarise(across(where(is.numeric), list(mean=mean))) %>% 
                              t() %>% 
                              as.data.frame() %>% 
                              rownames_to_column() %>%
                              `colnames<-`(.[1,]) %>%
                              .[-1,],
                            met.clusters2.prog.t2d,
                            met.clusters2.prog %>% 
                              group_by(ngr) %>% 
                              dplyr::summarise(across(where(is.numeric), list(mean=mean))) %>% 
                              t() %>% 
                              as.data.frame() %>% 
                              rownames_to_column() %>%
                              `colnames<-`(.[1,]) %>%
                              .[-1,],
                            met.clusters2.prog.ngr)
names(modules.stats) <- c('t2d', 'noT2D', 'T2D', 'fdr.t2d', 'ngr', 'NGR', 'noNGR', 'fdr.ngr')
modules.stats$module <- modules.stats$t2d; modules.stats$module <- gsub('_mean', '', modules.stats$module)
modules.stats$module.name <- sapply(modules.stats$module, function(a){
  paste0(a, ': ', cluster_mapping_file[which(cluster_mapping_file$New_Name==a),]$Description)
})

write.table(modules.stats, '/home/margar/download/export2_Metaboiteclusters_progressionT2DNGR.txt', sep='\t')

keep.modules <- unique(c(names(met.clusters2.prog.ngr[met.clusters2.prog.ngr<.1]), names(met.clusters2.prog.t2d[met.clusters2.prog.t2d<.1])))
met.clusters2.prog.sig <- met.clusters2.prog[, names(met.clusters2.prog)%in%keep.modules]

pheatmap::pheatmap(met.clusters2.prog.sig, annotation_row = progression.annotation, scale = 'row')
met.clusters2.prog.sig.m <- reshape2::melt(data.frame(met.clusters2.prog.sig, 'progression' = mtdt.all$progression))
met.clusters2.prog.sig.m$progression <- factor(met.clusters2.prog.sig.m$progression, levels = c('NGR', 'IGR', 'T2D'))
ggplot(met.clusters2.prog.sig.m, aes(progression, value)) + 
  geom_boxplot() + facet_wrap(~variable)

progression.annotation2 <- data.frame('extreme.progression'=mtdt.all[mtdt.all$progression!='IGR',]$progression); row.names(progression.annotation2) <- mtdt.all[mtdt.all$progression!='IGR',]$StudyID
names.sig.modules <- cluster_mapping_file[cluster_mapping_file$New_Name%in%names(met.clusters2.prog.sig),]
identical(names.sig.modules$New_Name, names(met.clusters2.prog.sig))
names.sig.modules <- names.sig.modules[match(names(met.clusters2.prog.sig), names.sig.modules$New_Name),]
names(met.clusters2.prog.sig) <- paste0(names.sig.modules$New_Name, ': ', names.sig.modules$Description)
pheatmap::pheatmap(met.clusters2.prog.sig[row.names(met.clusters2.prog.sig)%in%row.names(progression.annotation2),], 
                   annotation_row = progression.annotation2, scale='row', cutree_cols = c(2), cutree_rows = c(2))

modules.stats <- data.frame(met.clusters2.prog %>% 
                              dplyr::select(c(keep.modules, 't2d', 'ngr')) %>% 
                              group_by(t2d) %>% 
                              dplyr::summarise(across(where(is.numeric), list(mean=mean))) %>% 
                              t() %>% 
                              as.data.frame() %>% 
                              rownames_to_column() %>%
                              `colnames<-`(.[1,]) %>%
                              .[-1,] %>%
                              `rownames<-`(NULL),
                            met.clusters2.prog.t2d[names(met.clusters2.prog.t2d)%in%keep.modules],
                            met.clusters2.prog %>% 
                              dplyr::select(c(keep.modules, 't2d', 'ngr')) %>% 
                              group_by(ngr) %>% 
                              dplyr::summarise(across(where(is.numeric), list(mean=mean))) %>% 
                              t() %>% 
                              as.data.frame() %>% 
                              rownames_to_column() %>%
                              `colnames<-`(.[1,]) %>%
                              .[-1,] %>%
                              `rownames<-`(NULL),
                            met.clusters2.prog.ngr[names(met.clusters2.prog.ngr)%in%keep.modules]) 
names(modules.stats) <- c('t2d', 'noT2D', 'T2D', 'fdr.t2d', 'ngr', 'NGR', 'noNGR', 'fdr.ngr')
identical(modules.stats$t2d, modules.stats$ngr)
row.names(modules.stats) <- gsub('_mean', '', modules.stats$t2d)
modules.stats$t2d <- NULL; modules.stats$ngr <- NULL
row.names(modules.stats) <- sapply(row.names(modules.stats), function(z){
  paste0(z, ': ', cluster_mapping_file[grep(paste0('^',z,'$'), cluster_mapping_file$New_Name),]$Description)
})
# write.table(modules.stats, 'export_progressionassociatedmodules.txt', sep ='\t')

## metabolite cluster associations ----
# phenotype at baseline ----
mtdt.bsl.red <- data.frame('Age' = mtdt$Age, 'GR' = mtdt$DSGeneAbundRichness, 'Shannon' = mtdt$DSMGSAbundShannon,
                           mtdt[,15:56])
mtdt.bsl.red$Stumvoll[mtdt.bsl.red$Stumvoll<0] <- 0
mtdt.bsl.red[,4:45] <- log(mtdt.bsl.red[,4:45]+1)

bch.met.corr <- matrix (NA, nrow = ncol (mtdt.bsl.red), ncol = ncol (met.clusters2))
rownames (bch.met.corr) = names (mtdt.bsl.red)
colnames (bch.met.corr) = colnames (met.clusters2)
bch.met.corr.p <- matrix (NA, nrow = ncol (mtdt.bsl.red), ncol = ncol (met.clusters2))
rownames (bch.met.corr.p) = names (mtdt.bsl.red)
colnames (bch.met.corr.p) = colnames (met.clusters2)

for (m in colnames (met.clusters2)) {
  bch.met.corr [ , m] = apply (mtdt.bsl.red, MARGIN = 2, FUN = function (x) 
    cor (x, met.clusters2 [, m],
         method = cor_method, use = "pairwise.complete.obs"))
  bch.met.corr.p [ , m] = apply (mtdt.bsl.red, MARGIN = 2, FUN = function (x)
    cor.test (x, met.clusters2 [, m],
              method = cor_method, use = "pairwise.complete.obs")$p.val)
}

bch.met.corr <- as.data.frame(bch.met.corr)
bch.met.corr$bch <- row.names(bch.met.corr)
bch.met.corr.p <- as.data.frame(bch.met.corr.p)
bch.met.corr.p$bch <- row.names(bch.met.corr.p)

bch.met.corr.m <- reshape2::melt(bch.met.corr)
bch.met.corr.p.m <- reshape2::melt(bch.met.corr.p)
identical(bch.met.corr.m$bch, bch.met.corr.p.m$bch)
identical(bch.met.corr.m$variable, bch.met.corr.p.m$variable)

bch.met.corr.m$pval <- bch.met.corr.p.m$value
bch.met.corr.m$fdr <- p.adjust(bch.met.corr.m$pval, method ='fdr')
ggplot(bch.met.corr.m, aes(bch, variable, fill=value)) +
  geom_tile() + 
  scale_fill_gradient2(low='blue', high='red')

bch.met.corr.m.sig <- bch.met.corr.m[bch.met.corr.m$fdr <= .1,]

bch.met.corr.m.sig$module.name <- sapply(bch.met.corr.m.sig$variable, function(a){
  paste0(a, ': ', cluster_mapping_file[which(cluster_mapping_file$New_Name==a),]$Description)
})
bch.met.corr.m.sig$bch <- factor(bch.met.corr.m, levels=unique(bch.met.corr.m$bch))
ggplot(bch.met.corr.m.sig, aes(bch, module.name, fill=value)) +
  geom_tile() + 
  scale_fill_gradient2(low='blue', high='red') +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))

bch.met.corr.m$module.name <- sapply(bch.met.corr.m$variable, function(a){
  paste0(a, ': ', cluster_mapping_file[which(cluster_mapping_file$New_Name==a),]$Description)
})
write.table(bch.met.corr.m, '/home/margar/download/export2_correlationsmetclusters_BCHMBASELINE.txt', sep = '\t')

# phenotype progression (delta values) ----
delta.df <- readRDS('data/DeltaofLogValues.rds')
row.names(delta.df) <- delta.df$studyID
delta.df <- delta.df[,6:24]
identical(row.names(delta.df), row.names(met.clusters2.prog))

met.clusters2.prog <- met.clusters2.prog[,1:369]
delta.met.corr <- matrix (NA, nrow = ncol (delta.df), ncol = ncol (met.clusters2.prog))
rownames (delta.met.corr) = names (delta.df)
colnames (delta.met.corr) = colnames (met.clusters2.prog)
delta.met.corr.p <- matrix (NA, nrow = ncol (delta.df), ncol = ncol (met.clusters2.prog))
rownames (delta.met.corr.p) = names (delta.df)
colnames (delta.met.corr.p) = colnames (met.clusters2.prog)

for (m in colnames (met.clusters2.prog)) {
  delta.met.corr [ , m] = apply (delta.df, MARGIN = 2, FUN = function (x) 
    cor (x, met.clusters2.prog [, m],
         method = cor_method, use = "pairwise.complete.obs"))
  delta.met.corr.p [ , m] = apply (delta.df, MARGIN = 2, FUN = function (x)
    cor.test (x, met.clusters2.prog [, m],
              method = cor_method, use = "pairwise.complete.obs")$p.val)
}

delta.met.corr <- as.data.frame(delta.met.corr)
delta.met.corr$delta <- row.names(delta.met.corr)
delta.met.corr.p <- as.data.frame(delta.met.corr.p)
delta.met.corr.p$delta <- row.names(delta.met.corr.p)

delta.met.corr.m <- reshape2::melt(delta.met.corr)
delta.met.corr.p.m <- reshape2::melt(delta.met.corr.p)
identical(delta.met.corr.m$delta, delta.met.corr.p.m$delta)
identical(delta.met.corr.m$variable, delta.met.corr.p.m$variable)

delta.met.corr.m$pval <- delta.met.corr.p.m$value
delta.met.corr.m$fdr <- p.adjust(delta.met.corr.m$pval, method ='fdr')

delta.met.corr.m$module.name <- sapply(delta.met.corr.m$variable, function(a){
  paste0(a, ': ', cluster_mapping_file[which(cluster_mapping_file$New_Name==a),]$Description)
})


ggplot(delta.met.corr.m[delta.met.corr.m$fdr<.1,], aes(delta, module.name, fill=value)) +
  geom_tile() + 
  scale_fill_gradient2(low='blue', high='red') +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))

ggplot(delta.df, aes(Mean.Insulin, fill=mtdt.all$progression)) + geom_density(alpha=.4)
write.table(delta.met.corr.m, '/home/margar/download/export2_MetClusters_deltavalues.txt', sep='\t')

# GMM----
gmm <- read.table('data/WP2.1_775_IGR_GMM-QMP_igrHA_Ent_M.csv', header=T,row.names=1)
gmm <- gmm[,7:103]
row.names(gmm) == mtdt$samplerenamed
row.names(gmm) <- row.names(mtdt)

gmm.met.corr <- matrix (NA, nrow = ncol (gmm), ncol = ncol (met.clusters2))
rownames (gmm.met.corr) = names (gmm)
colnames (gmm.met.corr) = colnames (met.clusters2)
gmm.met.corr.p <- matrix (NA, nrow = ncol (gmm), ncol = ncol (met.clusters2))
rownames (gmm.met.corr.p) = names (gmm)
colnames (gmm.met.corr.p) = colnames (met.clusters2)

for (m in colnames (met.clusters2)) {
  gmm.met.corr [ , m] = apply (gmm, MARGIN = 2, FUN = function (x) 
    cor (x, met.clusters2 [, m],
         method = cor_method, use = "pairwise.complete.obs"))
  gmm.met.corr.p [ , m] = apply (gmm, MARGIN = 2, FUN = function (x)
    cor.test (x, met.clusters2 [, m],
              method = cor_method, use = "pairwise.complete.obs")$p.val)
}

gmm.met.corr <- as.data.frame(gmm.met.corr)
gmm.met.corr$gmm <- row.names(gmm.met.corr)
gmm.met.corr.p <- as.data.frame(gmm.met.corr.p)
gmm.met.corr.p$gmm <- row.names(gmm.met.corr.p)

gmm.met.corr.m <- reshape2::melt(gmm.met.corr)
gmm.met.corr.p.m <- reshape2::melt(gmm.met.corr.p)
identical(gmm.met.corr.m$gmm, gmm.met.corr.p.m$gmm)
identical(gmm.met.corr.m$variable, gmm.met.corr.p.m$variable)

gmm.met.corr.m$pval <- gmm.met.corr.p.m$value
gmm.met.corr.m$fdr <- p.adjust(gmm.met.corr.m$pval, method ='fdr')
ggplot(gmm.met.corr.m, aes(gmm, variable, fill=value)) +
  geom_tile() + 
  scale_fill_gradient2(low='blue', high='red')

gmm.annotations <- read.table('data/GMMs.v1.07.names', header=F, sep='\t')
gmm.met.corr.m$module.name <- sapply(gmm.met.corr.m$variable, function(a){
  paste0(a, ': ', cluster_mapping_file[which(cluster_mapping_file$New_Name==a),]$Description)
})
gmm.met.corr.m$gmm.name <- sapply(gmm.met.corr.m$gmm, function(a){
  gmm.annotations[which(gmm.annotations$V1==a),]$V2
})

write.table(gmm.met.corr.m, '/home/margar/download/export2_MetCluster_GMMs.txt', sep='\t')

gmm.met.corr.m.sig <- gmm.met.corr.m[gmm.met.corr.m$fdr <= .1,]

ggplot(gmm.met.corr.m.sig, aes(gmm.name, module.name, fill=value)) +
  geom_tile() + 
  scale_fill_gradient2(low='blue', high='red') +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))

## bacterial metabolites: we will define potential bacterial metabolites those metabolite clusters that are associated, at some extent, to microbial community
## this includes: gene richness tertiles, enterotype B2 and those resulting from the LASSO modelling done below ----
mgs <- read.table('data/WP2.1_775_IGR_MGS-QMP_igrHA_Ent_M.csv', header=T, row.names=1)
mgs <- mgs[,7:727]
identical(row.names(mgs), mtdt$samplerenamed)
row.names(mgs) <- row.names(mtdt)
identical(row.names(mgs), row.names(met.clusters2))

library(plyr)
library(readr)
library(dplyr)
library(caret)
library(ggplot2)
library(repr)
library(glmnet)

data.model <- cbind(mgs, met.clusters2)

set.seed(100) 
index <- sample(1:nrow(data.model), 0.8*nrow(data.model)) 
train <- data.model[index,] # Create the training data 
test <- data.model[-index,] # Create the test data


lambdas <- 10^seq(2, -3, by = -.1)
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

lasso.metaG <- lapply(met.clusters2, function(x){
  set.seed(100) 
  data.model <- cbind(mgs, 'met'=x)
  index <- sample(1:nrow(data.model), 0.8*nrow(data.model)) 
  train <- data.model[index,] # Create the training data 
  test <- data.model[-index,] # Create the test data
  lasso_reg <- cv.glmnet(as.matrix(train[,1:721]), train$met, alpha = 1, lambda = lambdas, standardize = TRUE, nfolds = 5)
  lasso_model <- glmnet(as.matrix(train[,1:721]), train$met, alpha = 1, lambda=lasso_reg$lambda.min, intercept = F)
  W <- as.matrix(coef(lasso_model))
  W
  keep_X <- rownames(W)[W!=0]
  keep_X <- keep_X[!keep_X == "(Intercept)"]
  train.X <- train[,keep_X]
  model1.data <- data.frame('met'=train$met, train.X)
  model1 <- lm(met~., data = model1.data)
  model1.sm <- summary(model1)
  r2 <- model1.sm$r.squared
  p.val <- if(length(coef(model1))>1){
    lmp(model1)} else {'NA'}
  return(data.frame('var'=r2, 'pval'=p.val, 'MGS' = paste(keep_X, collapse='|')))
})

lasso.metaG.df <- do.call(rbind.data.frame, lasso.metaG)
lasso.metaG.df$fdr <- p.adjust(lasso.metaG.df$pval, method = 'fdr')
lasso.metaG.df.taxa <- lasso.metaG.df[complete.cases(lasso.metaG.df),]
lasso.metaG.df$MGS <- NULL

lasso.metaG.df.taxa2 <- lasso.metaG.df.taxa$MGS
names(lasso.metaG.df.taxa2) <- row.names(lasso.metaG.df.taxa) 
df2<-  list() 
for(i in names(lasso.metaG.df.taxa2)){
  df <- lasso.metaG.df.taxa2[i]
  df <- strsplit(df, '\\|')
  df <- unlist(df)
  df <- data.frame('mgs' = df, 'met' = rep(i, length(df)))
  df2[[i]] <- df
}
metaG.mets.relation <- do.call(rbind.data.frame, df2)
metaG.mets.relation$module.name <- sapply(metaG.mets.relation$met, function(a){
  paste0(a, ': ', cluster_mapping_file[which(cluster_mapping_file$New_Name==a),]$Description)
})
write.table(metaG.mets.relation, '/home/margar/download/export2_MGScoeffsLASSOmetaboliteCLUSTERS.txt', sep ='\t', row.names = F)

biochem <- mtdt[,15:56]
inflam <- mtdt[,57:71]
lifestyle <- mtdt[,c(15,56:57)]

lasso.biochem <- lapply(met.clusters2, function(x){
  set.seed(101) 
  data.model <- cbind(biochem, 'met'=x)
  index <- sample(1:nrow(data.model), 0.8*nrow(data.model)) 
  train <- data.model[index,] # Create the training data 
  test <- data.model[-index,] # Create the test data
  lasso_reg <- cv.glmnet(as.matrix(train[,1:40]), train$met, alpha = 1, lambda = lambdas, standardize = TRUE, nfolds = 5)
  lasso_model <- glmnet(as.matrix(train[,1:40]), train$met, alpha = 1, lambda=lasso_reg$lambda.min, intercept = F)
  W <- as.matrix(coef(lasso_model))
  W
  keep_X <- rownames(W)[W!=0]
  keep_X <- keep_X[!keep_X == "(Intercept)"]
  train.X <- train[,keep_X]
  model1.data <- data.frame('met'=train$met, train.X)
  model1 <- lm(met~., data = model1.data)
  model1.sm <- summary(model1)
  r2 <- model1.sm$r.squared
  p.val <- if(length(coef(model1))>1){
    lmp(model1)} else {'NA'}
  return(data.frame('var'=r2, 'pval'=p.val))
})
lasso.biochem.df <- do.call(rbind.data.frame, lasso.biochem)
lasso.biochem.df$fdr <- p.adjust(lasso.biochem.df$pval, method = 'fdr')

lasso.infl <- lapply(met.clusters2, function(x){
  set.seed(102) 
  data.model <- cbind(inflam, 'met'=x)
  data.model <- data.model[complete.cases(data.model),]
  index <- sample(1:nrow(data.model), 0.8*nrow(data.model)) 
  train <- data.model[index,] # Create the training data 
  test <- data.model[-index,] # Create the test data
  lasso_reg <- cv.glmnet(as.matrix(train[,1:15]), train$met, alpha = 1, lambda = lambdas, standardize = TRUE, nfolds = 5)
  lasso_model <- glmnet(as.matrix(train[,1:15]), train$met, alpha = 1, lambda=lasso_reg$lambda.min, intercept = F)
  W <- as.matrix(coef(lasso_model))
  W
  keep_X <- rownames(W)[W!=0]
  keep_X <- keep_X[!keep_X == "(Intercept)"]
  train.X <- train[,keep_X]
  model1.data <- data.frame('met'=train$met, train.X)
  model1 <- lm(met~., data = model1.data)
  model1.sm <- summary(model1)
  r2 <- model1.sm$r.squared
  p.val <- if(length(coef(model1))>1){
    lmp(model1)} else {'NA'}
  return(data.frame('var'=r2, 'pval'=p.val))
})

lasso.infl.df <- do.call(rbind.data.frame, lasso.infl)
lasso.infl.df$fdr <- p.adjust(lasso.infl.df$pval, method = 'fdr')

fullmodel <- data.frame(mgs, biochem, inflam)
lasso.full <- lapply(met.clusters2, function(x){
  set.seed(104) 
  data.model <- cbind(fullmodel, 'met'=x)
  data.model <- data.model[complete.cases(data.model),]
  index <- sample(1:nrow(data.model), 0.8*nrow(data.model)) 
  train <- data.model[index,] # Create the training data 
  test <- data.model[-index,] # Create the test data
  lasso_reg <- cv.glmnet(as.matrix(train[,1:779]), train$met, alpha = 1, lambda = lambdas, standardize = TRUE, nfolds = 5)
  lasso_model <- glmnet(as.matrix(train[,1:779]), train$met, alpha = 1, lambda=lasso_reg$lambda.min, intercept = F)
  W <- as.matrix(coef(lasso_model))
  W
  keep_X <- rownames(W)[W!=0]
  keep_X <- keep_X[!keep_X == "(Intercept)"]
  train.X <- train[,keep_X]
  model1.data <- data.frame('met'=train$met, train.X)
  model1 <- lm(met~., data = model1.data)
  model1.sm <- summary(model1)
  r2 <- model1.sm$r.squared
  p.val <- if(length(coef(model1))>1){
    lmp(model1)} else {'NA'}
  return(data.frame('var'=r2, 'pval'=p.val))
})

lasso.full.df <- do.call(rbind.data.frame, lasso.full)
lasso.full.df$fdr <- p.adjust(lasso.full.df$pval, method = 'fdr')

lasso.models <- list('full' = lasso.full.df,
                     'metaG' = lasso.metaG.df,
                     'biochem' = lasso.biochem.df,
                     'inflam' = lasso.infl.df
)

lasso.models <- lapply(lasso.models, function(a){
  a <- a[!is.na(a$fdr),]
  a <- a[a$fdr<.1,]
})
lasso.models <- do.call(rbind.data.frame, lasso.models)
lasso.models$data.type <- unlist(strsplit(row.names(lasso.models), '\\.'))[c(T,F)]
lasso.models$metabolite <- unlist(strsplit(row.names(lasso.models), '\\.'))[c(F,T)]
lasso.models$data.type <- factor(lasso.models$data.type, levels = c('full', 'biochem', 'metaG', 'inflam'))


summary(lasso.models$data.type)
lasso.models$data.type <- factor(lasso.models$data.type, levels=c('full', 'biochem', 'metaG', 'inflam'))
lasso.models$data.type <- plyr::revalue(lasso.models$data.type, 
                                        c('full' = 'Combinedmodel',
                                          'biochem' = 'Bioclinical',
                                          'metaG' = 'Microbiome',
                                          'inflam' = 'Inflammatorymarkers'))
p1 <- lasso.models %>% 
  # filter(data.type!='lifestyle') %>% 
  ggplot(aes(data.type, var, color = data.type)) +
  geom_boxplot(fill=NA, color='black', outlier.colour = NA) +
  ggbeeswarm::geom_quasirandom(alpha=.4) + 
  theme(panel.grid.major = element_line(colour = "gray85",  linetype = "dashed"), 
        panel.background = element_rect(fill = NA),  legend.position = "none",
        panel.border = element_rect(fill=NA), 
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10, angle=90, hjust=1, vjust=.5)) +
  labs(x = NULL, y = "Explained variance (/1)", 
       colour = NULL)

ggExtra::ggMarginal(
  p = p1,
  type = 'density',
  margins = 'both',
  size = 5,
  groupColour = TRUE,
  groupFill = TRUE
)

lasso.models %>% 
  dplyr::mutate('data.type' = factor(data.type)) %>% 
  dplyr::select('data.type') %>% 
  summary()

#stats variance
lasso.models %>% 
  filter(data.type == 'Bioclinical') %>% 
  dplyr::select("var") %>% 
  summary()

lasso.models %>% 
  filter(data.type == 'Microbiome') %>% 
  dplyr::select('var') %>% 
  summary()

lasso.models %>% 
  filter(data.type == 'Inflammatorymarkers') %>% 
  dplyr::select('var') %>% 
  summary()

p2 <- lasso.models %>% 
  filter(data.type!='Combined\nmodel') %>% 
  ggplot(aes(data.type, metabolite, fill = var)) + 
  geom_tile() +
  scale_fill_gradient2(low='red', high='darkred') +
  theme(panel.grid.major = element_line(colour = "gray85",  linetype = "dashed"), 
        panel.background = element_rect(fill = NA),  legend.position = "none",
        panel.border = element_rect(fill=NA), 
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10, angle=90, hjust=1, vjust=.5)) +
  labs(x=NULL, y='Expl. var. of the metabolite clusters')

pdf('/home/margar/download/export2_LASSOplot.pdf', width = 8, height = 5)
cowplot::plot_grid(p1, p2)
dev.off()

lasso.models$metabolite2 <- sapply(lasso.models$metabolite, function(a){
  paste0(a, ': ', cluster_mapping_file[which(cluster_mapping_file$New_Name==a),]$Description)
})
write.table(lasso.models, '/home/margar/download/export2_explainedvarianceLASSOmodelsMETABOLITECLUSTERS.txt', sep ='\t')

## list of bacterial-associated modules
length(lasso.models[lasso.models$data.type=='metaG',]$metabolite)
length(keep.modules.micro)
length(unique(c(keep.modules.micro, lasso.models[lasso.models$data.type=='metaG',]$metabolite)))
bacterial.related.modules <- unique(c(keep.modules.micro, lasso.models[lasso.models$data.type=='metaG',]$metabolite))
write.table(bacterial.related.modules, 'bacterialmodules.txt')

met.clusters2.bacterialmod <- met.clusters2[,names(met.clusters2)%in%bacterial.related.modules]
# saveRDS(met.clusters2.bacterialmod, 'data/bacteriametclusters.rds')
modules.dist <- vegan::vegdist(met.clusters2.bacterialmod, method='euclidean')
pca.modules <- ape::pcoa(modules.dist)

# mgs.modules <- mgs[,names(mgs)%in%metaG.mets.relation$mgs]
# mgs2 <- sapply(mgs, function(x)(x/sum(x)))
mgs.dist <- vegan::vegdist(mgs, method='euclidean')
pca.mgs <- ape::pcoa(mgs.dist)

procr <- vegan::procrustes(pca.modules$vectors, pca.mgs$vectors, symmetric = T)
plot(procr, kind=1)
plot(procr, kind=2)
vegan::protest(pca.modules$vectors, pca.mgs$vectors, symmetric =T)

pro.df <- data.frame(x1 = procr$Yrot[,1], y1 = procr$Yrot[,2],
                     x2 = procr$X[,1], y2 = procr$X[,2], 'sample' = row.names(procr$Yrot))

ggplot(pro.df, aes(label = sample)) +
  geom_point(aes(x1, y1), colour = 'darkolivegreen2', shape = 19, size = 5) +
  geom_point(aes(x2, y2), colour = 'cornflowerblue', shape = 17, size = 5) +
  geom_segment(aes(x=x1, y=y1, xend=x2, yend=y2), arrow=arrow(length=unit(.3, 'cm')), alpha=.25) +
  geom_hline(yintercept = 0, linetype = 'dashed') + geom_vline(xintercept = 0, linetype = 'dashed') +
  theme(plot.subtitle = element_text(vjust = 1), plot.caption = element_text(vjust = 1),
        axis.title = element_text(size = 15, face = "bold"),axis.text = element_text(size = 14, face = "bold", colour = "black"),
        panel.background = element_rect(fill = NA), panel.border = element_rect(fill=NA, color = 'black', size = .5)) +
  labs(x = "Dimension 1", y = "Dimension 2") +
  annotate('text', x = 0.001, y =.001, label = 'Correlation 0.228 \n Monte-Carlo p = 0.001', fontface =2)

