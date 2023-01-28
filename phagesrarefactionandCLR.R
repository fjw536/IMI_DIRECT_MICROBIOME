## phages
setwd('/home/scratch/margar')
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
