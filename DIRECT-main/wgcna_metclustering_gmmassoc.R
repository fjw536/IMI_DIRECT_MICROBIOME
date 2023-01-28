###########################################################################################
#################### metabolite clustering as by helle's paper ############################
###########################################################################################

library (xlsx) ### Saving to spreadsheet
library (data.table) ### Fast read of large files into R
library (WGCNA) ### -	Clustering software. Previously reported work done using v1.34
library (flashClust) ### Clustering software
library (ppcor) ### Partial Spearman correlations, for confounder analysis. Previously reported work done using v1.0
library (gplots) ### Plotting
library (cowplot) ### Plotting; to arrange several plots on the same page
library (ggplot2) ### Plotting
library (plyr) ### Data transformations

# directory
setwd('/home/scratch/margar')

# data
mets <- readRDS('data/metabolomics/MetabolitesResiduals.rds')
mtdt <- read.table('data/curatedmetadata/WP2.1_BL_MetadataWithInfl2.txt', header =T, row.names=1, sep='\t')
annot <- read.table('data/metabolomics/Chemical_Annotation_Metabolone_DTU_sh.txt', header=T, row.names=1)

## settings for WGCNA
cor_method          = "spearman" ### for association with clinical parameters
corFun_tmp          = "bicor"
cluster_method      = "average"
corOptions_list     = list (use = 'pairwise.complete.obs') 
corOptions_str      = "use = 'pairwise.complete.obs'"
BH_pval_asso_cutoff = 0.05
NetworkType         = "signed" ### Signed-network (as the PC1 and median profile does not make sense as a summary measure of a cluster with anticorrelated metabolites.)

### Specify data and parameters 
identical(row.names(mets), row.names(mtdt)) #check we have same samples, same order

### Settings for WGCNA on polar metabolite measurements
### Once these are established the steps below can be run
RsquareCut_val      = 0.89 ### usually ranges 0.80-0.95 but requires inspecting curves
mergingThresh       = 0.20 ### Maximum dissimilarity of module eigengenes (i.e. 1-correlation) for merging modules.
minModuleSize       = 3 ### minimum number of metabolites constituting a cluster
SoftPower           = 13 ### beta-value, main parameter to optimize


## Calculate weighted adjacency matrix
A <- adjacency(mets, power = SoftPower, type = NetworkType, corFnc = corFun_tmp, corOptions = corOptions_str)
colnames (A) <- rownames (A) <- colnames (mets)
### Define dissimilarity based on topological overlap
dissTOM = TOMdist (A, TOMType = NetworkType)
colnames (dissTOM) = rownames (dissTOM) = colnames (mets)
### Hierarchical clustering
metaTree = flashClust (as.dist (dissTOM), method = cluster_method)
### Define modules by cutting branches
moduleLabels1 = cutreeDynamic (dendro = metaTree, distM = dissTOM, method = "hybrid", deepSplit = 4, pamRespectsDendro = T, minClusterSize = minModuleSize)
moduleLabels1 = labels2colors (moduleLabels1)
### Automatically merge highly correlated modules
merge = mergeCloseModules (mets, moduleLabels1, corFnc = corFun_tmp, corOptions = corOptions_list, cutHeight = mergingThresh)
### Determine resulting merged module colors
moduleLabels2 = merge$colors
### Establish eigengenes of the newly merged modules, used for cluster overall abundances
MEs = merge$newMEs
### Choose final module assignments
moduleColorsMeta = moduleLabels2
names (moduleColorsMeta) = colnames (mets)
MEsMeta = orderMEs (MEs)
rownames (MEsMeta) = rownames (mets)

### Determine relevant descriptive statistics of established clusters
### kIN: within-module connectivity, determined by summing connectivity with all
###      other metabolites in the given cluster.
### kME: bicor-correlation between the metabolite profile and module eigenvector; 
### both measures of intramodular hub-metabolite status.
kIN <-      vector (length = ncol (mets)); names (kIN) = colnames (mets)
kME <-      vector (length = ncol (mets)); names (kME) = colnames (mets)
modules <-  vector (length = ncol (mets)); names (modules) = colnames (mets)

for (module in names (table (moduleColorsMeta))) {   
  
  all.metabolites = names (mets)
  inModule = (moduleColorsMeta == module)
  module.metabolites = names (moduleColorsMeta [inModule])
  modules [module.metabolites] = module 
  kIN [module.metabolites] = sapply (module.metabolites, function (x) sum (A [x, module.metabolites]) - 1)
  datKME = signedKME (mets, MEsMeta, corFnc = corFun_tmp, corOptions = corOptions_str)
  rownames (datKME) = colnames (mets)
  kME [module.metabolites] = datKME [module.metabolites, paste ("kME", module, sep = "")]   
  
}
output <- data.frame("module" = modules, "kME" = kME, "kIN" = kIN)
head(output)
identical(row.names(output), row.names(annot))
output <- data.frame(output, 'name' = annot$CHEMICAL_NAME, 'type' = annot$SUPER_PATHWAY, 'subtype' = annot$SUB_PATHWAY)
# write.table(output, 'data/metabolomics/modulemappingfile.txt') ## annotate clusters manually
cluster_mapping_file <- read.table('data/metabolomics/moduleannotationscurated.txt', header=T, sep='\t')
cluster_mapping_file$Module_Name <- paste0('ME', cluster_mapping_file$Module_Name) 
output = data.frame ("module" = modules, "kME" = kME, "kIN" = kIN, 
                     "cluster_name" = sapply (modules, function (m) cluster_mapping_file [which(cluster_mapping_file$Module_Name==m), "New_Name"]),
                     'cluster_description' =  sapply (modules, function (m) cluster_mapping_file [which(cluster_mapping_file$Module_Name==m), "Description"]))

# save.image(file='metclustersrun.RData')
#rename cluster names from color to module number
colnames(MEsMeta) <- sapply (colnames(MEsMeta), function (m) cluster_mapping_file [which(cluster_mapping_file$Module_Name==m), "New_Name"])
MEsMeta <- MEsMeta[,order(colnames(MEsMeta))]
# saveRDS(MEsMeta, 'data/metabolomics/metabolitemodules.rds')

# how many clusters and how many metabolites per cluster
MEsMeta <- readRDS('data/metabolomics/metabolitemodules.rds')
cluster_mapping_file <- read.table('data/metabolomics/moduleannotationscurated.txt', header=T, sep='\t')
library(tidyverse)

mod.order <- factor(names(sort(summary(factor(output$cluster_name)))), 
                       levels = names(sort(summary(factor(output$cluster_name)))))
output$cluster_name <- factor(output$cluster_name, levels = mod.order)
ggplot(output, aes(cluster_name)) + geom_bar()

## do clusters associate to Gene Richness and/or Bact2?
annotation.pheat <- data.frame('GR'=mtdt$DSGARLTertileGroup, 'entero' = mtdt$enterotype); row.names(annotation.pheat) <- row.names(mtdt)
pheatmap::pheatmap(MEsMeta, 
                   annotation_row = annotation.pheat, scale='row')#, cutree_cols = c(2), cutree_rows = c(2))

MEsMeta.micro <- MEsMeta
MEsMeta.micro$entero <- plyr::revalue(mtdt$enterotype, c('Bact1'='noBact2', 'Prev'='noBact2', 'Rum'= 'noBact2'))
MEsMeta.micro.ent <- sapply(MEsMeta.micro[,1:113], function(x)(
  t.test(x~MEsMeta.micro$entero)$p.val
))
MEsMeta.micro.ent <- p.adjust(MEsMeta.micro.ent, method = 'fdr')
met.clusters.ent <- MEsMeta.micro %>% 
  group_by(entero) %>% 
  summarise(across(where(is.numeric), mean)) %>% 
  t() %>% 
  as.data.frame()
names(met.clusters.ent) <- met.clusters.ent[1,] ; met.clusters.ent <- met.clusters.ent[-1,]
row.names(met.clusters.ent) == names(MEsMeta.micro.ent)
met.clusters.ent$fdr <- MEsMeta.micro.ent

annotation.pheat <- data.frame('GR'=mtdt$DSGARLTertileGroup, 'entero' = MEsMeta.micro$entero); row.names(annotation.pheat) <- row.names(mtdt)
pheatmap::pheatmap(MEsMeta.micro[,names(MEsMeta.micro)%in%names(MEsMeta.micro.ent[MEsMeta.micro.ent<.1])], 
                   annotation_row = annotation.pheat, scale='row')

MEsMeta.micro$GR <- factor(mtdt$DSGARLTertileGroup)
MEsMeta.micro.gr <- sapply(MEsMeta.micro[MEsMeta.micro$GR!='HGC',1:113], function(x)(
  t.test(x~droplevels(MEsMeta.micro[MEsMeta.micro$GR!="HGC",]$GR))$p.val
))
MEsMeta.micro.gr <- p.adjust(MEsMeta.micro.gr, method = 'fdr')
met.clusters.gr <- MEsMeta.micro %>% 
  group_by(GR) %>% 
  summarise(across(where(is.numeric), mean)) %>% 
  t() %>% 
  as.data.frame()
names(met.clusters.gr) <- met.clusters.gr[1,] ; met.clusters.gr <- met.clusters.gr[-1,]
row.names(met.clusters.gr) == names(MEsMeta.micro.ent)
met.clusters.gr$fdr <- MEsMeta.micro.gr

write.table(cbind(met.clusters.ent, met.clusters.gr), '/home/margar/download/export_MEtClustersMean_GRandENtero.txt', sep='\t')
pheatmap::pheatmap(MEsMeta.micro[MEsMeta.micro$GR!='HGC',names(MEsMeta.micro)%in%names(MEsMeta.micro.ent[MEsMeta.micro.gr<.1])], 
                   annotation_row = annotation.pheat, scale='row')

keep.modules.micro <- unique(names(MEsMeta.micro.gr[MEsMeta.micro.gr<.1]), names(MEsMeta.micro.ent[MEsMeta.micro.ent<.1]))
pheatmap::pheatmap(MEsMeta[, names(MEsMeta)%in%keep.modules.micro], 
                   annotation_row = annotation.pheat, scale='row')

## do clusters associate to progression'? ----
mtdt.all <- read.table('data/curatedmetadata/WP2.1_BlFw-641_MetadataImputed_OnlyBIOCHEM_2.txt', header=T, sep = '\t')
mtdt.all[mtdt.all$timepoint=='Baseline',]$StudyID == mtdt.all[mtdt.all$timepoint!='Baseline',]$StudyID 
mtdt.all$progression <- rep(mtdt.all[mtdt.all$timepoint!='Baseline',]$Gly.Cat, 2)
mtdt.all <- mtdt.all[mtdt.all$timepoint=='Baseline',]
MEsMeta.prog <- MEsMeta[row.names(MEsMeta)%in%mtdt.all$StudyID,]
identical(row.names(MEsMeta.prog), mtdt.all$StudyID)
progression.annotation <- data.frame(mtdt.all$progression); row.names(progression.annotation) <- mtdt.all$StudyID
pheatmap::pheatmap(MEsMeta.prog, annotation_row = progression.annotation, scale = 'row')

## progression differential metabolites ----
MEsMeta.prog$t2d <- plyr::revalue(mtdt.all$progression, c('IGR'='noT2D', 'NGR'='noT2D'))
MEsMeta.prog$ngr <- plyr::revalue(mtdt.all$progression, c('IGR'='noNGR', 'T2D'='noNGR'))

sapply(MEsMeta.prog[,1:113], function(x)(
  print(hist(x))
))

MEsMeta.prog.t2d <- sapply(MEsMeta.prog[,1:113], function(x)(
  t.test(x~MEsMeta.prog$t2d)$p.val
))
MEsMeta.prog.t2d <- p.adjust(MEsMeta.prog.t2d, method = 'fdr')

MEsMeta.prog.ngr <- sapply(MEsMeta.prog[,1:113], function(x)(
  t.test(x~MEsMeta.prog$ngr)$p.val
))
MEsMeta.prog.ngr <- p.adjust(MEsMeta.prog.ngr, method = 'fdr')

library(tibble)
modules.stats <- data.frame(MEsMeta.prog %>% 
                              group_by(t2d) %>% 
                              dplyr::summarise(across(where(is.numeric), list(mean=mean))) %>% 
                              t() %>% 
                              as.data.frame() %>% 
                              rownames_to_column() %>%
                              `colnames<-`(.[1,]) %>%
                              .[-1,],
                            MEsMeta.prog.t2d,
                            MEsMeta.prog %>% 
                              group_by(ngr) %>% 
                              dplyr::summarise(across(where(is.numeric), list(mean=mean))) %>% 
                              t() %>% 
                              as.data.frame() %>% 
                              rownames_to_column() %>%
                              `colnames<-`(.[1,]) %>%
                              .[-1,],
                            MEsMeta.prog.ngr)
names(modules.stats) <- c('t2d', 'noT2D', 'T2D', 'fdr.t2d', 'ngr', 'NGR', 'noNGR', 'fdr.ngr')
write.table(modules.stats, '/home/margar/download/export_Metaboiteclusters_progressionT2DNGR.txt', sep='\t')

keep.modules <- unique(c(names(MEsMeta.prog.ngr[MEsMeta.prog.ngr<.1]), names(MEsMeta.prog.t2d[MEsMeta.prog.t2d<.1])))
MEsMeta.prog.sig <- MEsMeta.prog[, names(MEsMeta.prog)%in%keep.modules]

pheatmap::pheatmap(MEsMeta.prog.sig, annotation_row = progression.annotation, scale = 'row')
MEsMeta.prog.sig.m <- reshape2::melt(data.frame(MEsMeta.prog.sig, 'progression' = mtdt.all$progression))
MEsMeta.prog.sig.m$progression <- factor(MEsMeta.prog.sig.m$progression, levels = c('NGR', 'IGR', 'T2D'))
ggplot(MEsMeta.prog.sig.m, aes(progression, value)) + 
  geom_boxplot() + facet_wrap(~variable)

progression.annotation2 <- data.frame('extreme.progression'=mtdt.all[mtdt.all$progression!='IGR',]$progression); row.names(progression.annotation2) <- mtdt.all[mtdt.all$progression!='IGR',]$StudyID
names.sig.modules <- cluster_mapping_file[cluster_mapping_file$New_Name%in%names(MEsMeta.prog.sig),]
identical(names.sig.modules$New_Name, names(MEsMeta.prog.sig))
names.sig.modules <- names.sig.modules[match(names(MEsMeta.prog.sig), names.sig.modules$New_Name),]
names(MEsMeta.prog.sig) <- paste0(names.sig.modules$New_Name, ': ', names.sig.modules$Description)
pheatmap::pheatmap(MEsMeta.prog.sig[row.names(MEsMeta.prog.sig)%in%row.names(progression.annotation2),], 
                   annotation_row = progression.annotation2, scale='row', cutree_cols = c(2), cutree_rows = c(2))

modules.stats <- data.frame(MEsMeta.prog %>% 
    dplyr::select(c(keep.modules, 't2d', 'ngr')) %>% 
    group_by(t2d) %>% 
    dplyr::summarise(across(where(is.numeric), list(mean=mean))) %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column() %>%
    `colnames<-`(.[1,]) %>%
    .[-1,] %>%
    `rownames<-`(NULL),
    MEsMeta.prog.t2d[names(MEsMeta.prog.t2d)%in%keep.modules],
  MEsMeta.prog %>% 
    dplyr::select(c(keep.modules, 't2d', 'ngr')) %>% 
    group_by(ngr) %>% 
    dplyr::summarise(across(where(is.numeric), list(mean=mean))) %>% 
      t() %>% 
    as.data.frame() %>% 
    rownames_to_column() %>%
    `colnames<-`(.[1,]) %>%
    .[-1,] %>%
    `rownames<-`(NULL),
  MEsMeta.prog.ngr[names(MEsMeta.prog.ngr)%in%keep.modules]) 
names(modules.stats) <- c('t2d', 'noT2D', 'T2D', 'fdr.t2d', 'ngr', 'NGR', 'noNGR', 'fdr.ngr')
identical(modules.stats$t2d, modules.stats$ngr)
row.names(modules.stats) <- gsub('_mean', '', modules.stats$t2d)
modules.stats$t2d <- NULL; modules.stats$ngr <- NULL
row.names(modules.stats) <- sapply(row.names(modules.stats), function(z){
  paste0(z, ': ', cluster_mapping_file[grep(paste0('^',z,'$'), cluster_mapping_file$New_Name),]$Description)
})
# write.table(modules.stats, 'export_progressionassociatedmodules.txt', sep ='\t')

## metabolite cluster associations
# phenotype at baseline ----
mtdt.bsl.red <- data.frame('Age' = mtdt$Age, 'GR' = mtdt$DSGeneAbundRichness, 'Shannon' = mtdt$DSMGSAbundShannon,
                          mtdt[,15:56])
mtdt.bsl.red$Stumvoll[mtdt.bsl.red$Stumvoll<0] <- 0
mtdt.bsl.red[,4:45] <- log(mtdt.bsl.red[,4:45]+1)

bch.met.corr <- matrix (NA, nrow = ncol (mtdt.bsl.red), ncol = ncol (MEsMeta))
rownames (bch.met.corr) = names (mtdt.bsl.red)
colnames (bch.met.corr) = colnames (MEsMeta)
bch.met.corr.p <- matrix (NA, nrow = ncol (mtdt.bsl.red), ncol = ncol (MEsMeta))
rownames (bch.met.corr.p) = names (mtdt.bsl.red)
colnames (bch.met.corr.p) = colnames (MEsMeta)

for (m in colnames (MEsMeta)) {
  bch.met.corr [ , m] = apply (mtdt.bsl.red, MARGIN = 2, FUN = function (x) 
    cor (x, MEsMeta [, m],
         method = cor_method, use = "pairwise.complete.obs"))
  bch.met.corr.p [ , m] = apply (mtdt.bsl.red, MARGIN = 2, FUN = function (x)
    cor.test (x, MEsMeta [, m],
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

write.table(bch.met.corr.m, '/home/margar/download/export_correlationsmetclusters_BCHMBASELINE.txt', sep = '\t')

# phenotype progression (delta values) ----
delta.df <- readRDS('data/DeltaofLogValues.rds')
row.names(delta.df) <- delta.df$studyID
delta.df <- delta.df[,6:24]
identical(row.names(delta.df), row.names(MEsMeta.prog))

MEsMeta.prog <- MEsMeta.prog[,1:113]
delta.met.corr <- matrix (NA, nrow = ncol (delta.df), ncol = ncol (MEsMeta.prog))
rownames (delta.met.corr) = names (delta.df)
colnames (delta.met.corr) = colnames (MEsMeta.prog)
delta.met.corr.p <- matrix (NA, nrow = ncol (delta.df), ncol = ncol (MEsMeta.prog))
rownames (delta.met.corr.p) = names (delta.df)
colnames (delta.met.corr.p) = colnames (MEsMeta.prog)

for (m in colnames (MEsMeta.prog)) {
  delta.met.corr [ , m] = apply (delta.df, MARGIN = 2, FUN = function (x) 
    cor (x, MEsMeta.prog [, m],
         method = cor_method, use = "pairwise.complete.obs"))
  delta.met.corr.p [ , m] = apply (delta.df, MARGIN = 2, FUN = function (x)
    cor.test (x, MEsMeta.prog [, m],
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
write.table(delta.met.corr.m, '/home/margar/download/export_MetClusters_deltavalues.txt', sep='\t')

# GMM----
gmm <- read.table('data/WP2.1_775_IGR_GMM-QMP_igrHA_Ent_M.csv', header=T,row.names=1)
gmm <- gmm[,7:103]
row.names(gmm) == mtdt$samplerenamed
row.names(gmm) <- row.names(mtdt)

gmm.met.corr <- matrix (NA, nrow = ncol (gmm), ncol = ncol (MEsMeta))
rownames (gmm.met.corr) = names (gmm)
colnames (gmm.met.corr) = colnames (MEsMeta)
gmm.met.corr.p <- matrix (NA, nrow = ncol (gmm), ncol = ncol (MEsMeta))
rownames (gmm.met.corr.p) = names (gmm)
colnames (gmm.met.corr.p) = colnames (MEsMeta)

for (m in colnames (MEsMeta)) {
  gmm.met.corr [ , m] = apply (gmm, MARGIN = 2, FUN = function (x) 
    cor (x, MEsMeta [, m],
         method = cor_method, use = "pairwise.complete.obs"))
  gmm.met.corr.p [ , m] = apply (gmm, MARGIN = 2, FUN = function (x)
    cor.test (x, MEsMeta [, m],
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

write.table(gmm.met.corr.m, '/home/margar/download/export_MetCluster_GMMs.txt', sep='\t')

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
identical(row.names(mgs), row.names(MEsMeta))

library(plyr)
library(readr)
library(dplyr)
library(caret)
library(ggplot2)
library(repr)
library(glmnet)

data.model <- cbind(mgs, MEsMeta)

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

lasso.metaG <- lapply(MEsMeta, function(x){
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
write.table(metaG.mets.relation, 'MGScoeffsLASSOmetaboliteCLUSTERS.txt', sep ='\t', row.names = F)

biochem <- mtdt[,15:56]
inflam <- mtdt[,57:71]
lifestyle <- mtdt[,c(15,56:57)]

lasso.biochem <- lapply(MEsMeta, function(x){
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

lasso.infl <- lapply(MEsMeta, function(x){
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
lasso.full <- lapply(MEsMeta, function(x){
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
                                        c('full' = 'Combined\nmodel',
                                          'biochem' = 'Bioclinical',
                                          'metaG' = 'Microbiome',
                                          'inflam' = 'Inflammatory\nmarkers'))
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
  filter(data.type == 'biochem') %>% 
  dplyr::select("var") %>% 
  summary()

lasso.models %>% 
  filter(data.type == 'metaG') %>% 
  dplyr::select('var') %>% 
  summary()

lasso.models %>% 
  filter(data.type == 'inflam') %>% 
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

pdf('/home/margar/download/export_LASSOplot.pdf', width = 8, height = 5)
cowplot::plot_grid(p1, p2)
dev.off()

write.table(lasso.models, 'explainedvarianceLASSOmodelsMETABOLITECLUSTERS.txt', sep ='\t')

## list of bacterial-associated modules
length(lasso.models[lasso.models$data.type=='metaG',]$metabolite)
length(keep.modules.micro)
length(unique(c(keep.modules.micro, lasso.models[lasso.models$data.type=='metaG',]$metabolite)))
bacterial.related.modules <- unique(c(keep.modules.micro, lasso.models[lasso.models$data.type=='metaG',]$metabolite))
write.table(bacterial.related.modules, 'bacterialmodules.txt')

MEsMeta.bacterialmod <- MEsMeta[,names(MEsMeta)%in%bacterial.related.modules]
# saveRDS(MEsMeta.bacterialmod, 'data/bacteriametclusters.rds')
modules.dist <- vegan::vegdist(MEsMeta.bacterialmod, method='euclidean')
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

