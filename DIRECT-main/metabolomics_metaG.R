######################################################################################################
############################# Metabolomics-microbiome integration ####################################
######################################################################################################

## directory
setwd('/home/scratch/margar')

## libraries
library(factoextra); library(mdatools); library(tidyverse); library(purrr); library(vroom); library(fs)
library(ggplot2); library(ggpubr); library(ggbeeswarm); library(ROCR); library(corrplot); library(vegan)
library(Maaslin2)

## upload metabolomics data and QMP values
mets <- readRDS('data/metabolomics/MetabolitesResiduals.rds')

mtdt.bsl <- read.table('data/curatedmetadata/WP2.1_BL_MetadataWithInfl.txt', header=T, sep='\t')
qmp <- read.table('data/WP2.1_775_IGR_MGS-QMP_igrHA_Ent_M.csv', header=T)
qmp$ID == mtdt.bsl$samplerenamed
row.names(qmp) <- mtdt.bsl$StudyID
qmp <- qmp[,8:728]
gmm <- read.table('data/WP2.1_775_IGR_GMM-QMP_igrHA_Ent_M.csv', header=T, sep='\t', row.names=1)
gmm <- gmm[,7:103]
row.names(gmm) == mtdt.bsl$samplerenamed
row.names(gmm) <- mtdt.bsl$StudyID
gmm <- data.frame(gmm, 'Age' = mtdt.bsl$Age, 'Center' = mtdt.bsl$CenterID, 'Gender' = mtdt.bsl$Gender)
gmm.residuals <-  list()
for (module in seq(1:97)){
  model <- lm(gmm[,module] ~ Age + Center + Gender, data = gmm)
  residuals.gmm <- residuals(model)
  gmm.residuals[[module]]<- residuals.gmm
}

gmm.residuals <- do.call(cbind.data.frame, gmm.residuals)
names(gmm.residuals) <- names(gmm)[1:97]

## similarity analysis ----
## procrustes
mets.dist <- vegdist(mets, method = 'euclidean')
pca.mets <- prcomp(mets.dist, center=T, scale=T)
qmp.dist <- vegdist(qmp, method = 'euclidean')
pca.qmp <- prcomp(qmp.dist, center=T, scale=T)

proc.res <- procrustes(pca.mets, pca.qmp)
plot(proc.res)
protest(pca.mets, pca.qmp, permutations = 5000)


gmm.dist <- vegdist(gmm.residuals, method='euclidean')
pca.gmm <- prcomp(gmm.dist, center=T, scale =T)
proc.gmm <- procrustes(pca.mets, pca.gmm)
plot(proc.gmm)
protest(pca.mets, pca.gmm)

## pls-R ====
#1.- divide data
sample.split <- caTools::sample.split(row.names(mets), SplitRatio = .8)
mets.train <- subset(mets, sample.split==T)
qmp.train <- subset(qmp, sample.split==T)

mets.test <- subset(mets, sample.split==F)
qmp.test <- subset(qmp, sample.split==F)

pls.r.model <- pls(qmp.train, mets.train, 15, scale=T, center=T, cv = 1, info = 'Metabolites prediction analysis')
model.reduced <- selectCompNum(pls.r.model, 3)
print(model.reduced)
plotRegcoeffs(model.reduced)
predicted <- predict(model.reduced, qmp.test, mets.test)

## correlation analysis (Spearman's) ----
row.names(mets) == row.names(qmp)
mets.qmp <- cbind(mets, qmp)

cor.matrix <- cor(mets.qmp, method = 'spearman') 
cor.matrix.p <- cor.mtest(mets.qmp, method = 'spearman')
pdf('metabolon_microbiome/correlationspearman.pdf', width=25, height = 18)
corrplot(cor.matrix[1:103, 104:1694], p.mat = cor.matrix.p$p[1:103, 104:1694], insig = 'blank', method = 'square')
dev.off()
made4::heatplot(cor.matrix)

cor.matrix2 <- cor.matrix[1:973, 974:1694]
cor.matrix2p <- cor.matrix.p$p[1:973, 974:1694]

cor.matrix.export <- cor.matrix2 %>% as.data.frame %>% tibble::rownames_to_column() %>% 
  tidyr::pivot_longer(-rowname)
cor.matrix.export$pvalue <- as.data.frame(cor.matrix2p %>% as.data.frame %>% tibble::rownames_to_column() %>% 
  tidyr::pivot_longer(-rowname) %>% select(value))$value
cor.matrix.export <- as.data.frame(cor.matrix.export); names(cor.matrix.export) <- c('metabolite', 'bacteria', 'corr', 'pvalue')
cor.matrix.export <- cor.matrix.export[cor.matrix.export$pvalue<.05, ]
write.table(cor.matrix.export, '/home/margar/download/cormatrixDIRECTfull.txt', sep ='\t')
annotation.corrs <- data.frame('var' = names(mets.qmp), 'type' = c(rep('metabolite', 973), rep('bacteria', 721)))
write.table(annotation.corrs, '/home/margar/download/annotationcorrmatrixDIRECTfull.txt', sep ='\t')
tax <- read.table('data/taxonomy.tsv', header=T, sep ='\t')
tax$X <- gsub(':', '\\.', tax$X)
tax <- tax[tax$X%in%annotation.corrs$var,]
write.table(tax, '/home/margar/download/tax.txt', sep ='\t')

row.names(mets) == row.names(gmm.residuals)
mets.gmm <- cbind(mets, gmm.residuals)
cor.gmm <- cor(mets.gmm, method='spearman')
cor.gmm.p <- cor.mtest(mets.gmm, method='spearman')  
made4::heatplot(cor.gmm)

cor.gmm2 <- cor.gmm[1:973, 974:1070]
cor.gmm2p <- cor.gmm.p$p[1:973, 974:1070]
cor.matrix.export <- cor.gmm2 %>% as.data.frame %>% tibble::rownames_to_column() %>% 
  tidyr::pivot_longer(-rowname)
cor.matrix.export$pvalue <- as.data.frame(cor.gmm2p %>% as.data.frame %>% tibble::rownames_to_column() %>% 
                                            tidyr::pivot_longer(-rowname) %>% select(value))$value
cor.matrix.export <- as.data.frame(cor.matrix.export); names(cor.matrix.export) <- c('metabolite', 'GMM', 'corr', 'pvalue')
cor.matrix.export <- cor.matrix.export[cor.matrix.export$pvalue<.05, ]
write.table(cor.matrix.export, '/home/margar/download/cormatrixDIRECT_GMMfull.txt', sep ='\t')
annotation.corrs <- data.frame('var' = names(mets.gmm), 'type' = c(rep('metabolite', 973), rep('GMM', 97)))
write.table(annotation.corrs, '/home/margar/download/annotationcorrmatrixDIRECT_GMMfull.txt', sep ='\t')


#onlysignificant

## MAaAsLin2 ----
# see sublime text script

## Canonical Correlation Analysis
library(PMA)
set.seed(331)
perm.out <- CCA.permute(mets,qmp,typex="standard",typez="standard",nperms=7)
print(perm.out)
plot(perm.out)
out <- CCA(mets, qmp,typex="standard",typez="standard",K=1,
           penaltyx=perm.out$bestpenaltyx,penaltyz=perm.out$bestpenaltyz,
           v=perm.out$v.init)
print(out)
## Co-occurence network -----

