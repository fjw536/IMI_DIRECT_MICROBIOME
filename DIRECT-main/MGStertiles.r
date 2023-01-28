######### 
## metabolite associations to various things
  # gene richness in continuous mode
  # alpha diversity 
  # MGS richness in both tertiles and continuous mode

# directory
setwd('/home/scratch/margar')

# libraries
library(ggplot2); library(ggpmisc); library(ggbeeswarm); library(ggpubr); library(tidyverse)

## Is gene richness associated to reduced MGS? ----
mtdt <- read.table('data/curatedmetadata/WP2.1_BL_MetadataWithInfl2.txt', header=T, sep='\t')
qmp.bsl <- read.table('data/WP2.1_775_IGR_MGS-QMP_igrHA_Ent_M.csv', header=T, sep='\t', row.names=1)
qmp.bsl <- qmp.bsl[,7:727]

qmp.bsl <- data.frame(t(qmp.bsl))
mgs.rich <- colSums(qmp.bsl!=0)
names(mgs.rich) == mtdt$samplerenamed

#
data.compare <- data.frame('MGS.rich' = mgs.rich, 'gene.rich' = mtdt$DSGeneAbundRichness, 'shannon' = mtdt$DSMGSAbundShannon,
                           'Age' = mtdt$Age, 'Gender' = mtdt$Gender, 'CenterID' = mtdt$CenterID)
row.names(data.compare) <- mtdt$X 

rm(list=setdiff(ls(), c('data.compare', 'cor.associations')))
met.clusters <- readRDS('data/metabolomics/mergedMetClusters.rds')
cluster_mapping_file <- read.table('data/metabolomics/moduleannotationscurated2.txt', header=T, sep='\t')
# other.stats <- lm.associations(data.compare, met.clusters, names(data.compare)[1:3])
other.stats2 <- cor.associations(data.compare, met.clusters, names(data.compare)[1:3])
other.stats2 <- do.call(rbind.data.frame, other.stats2)
  
other.stats2 %>% 
  filter(fdr<.1) %>% 
  mutate('rho' = as.numeric(rho)) %>% 
  ggplot(aes(variable, metabolite, fill = rho)) + 
  geom_tile() + coord_flip() + 
  scale_fill_gradient2(low='blue', high='red') +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5), panel.border=element_rect(fill=NA),
        panel.grid.major = element_line(linetype='dashed', color = 'gray75')) +
  labs(y = 'Metabolite Clusters', x ='Richness variables')

other.stats2$metabolite2 <- sapply(other.stats2$metabolite, function(a){
  paste0(a, ': ', cluster_mapping_file[which(cluster_mapping_file$New_Name==a),]$Description)
})
write.table(other.stats2, '/home/margar/download/export2_richnessvariables_correlationMEtClusters.txt', sep='\t')

### MGS tertiles analysis - categorical ----
mgs.tertiles <- ntile(data.compare$MGS.rich, n=3)
data.compare$MGS.tertile <- plyr::revalue(factor(mgs.tertiles), c('1' = 'LowMGS', '2' = 'MidMGS', '3' = 'HighMGS'))

output.MGStert <- lapply(met.clusters, function(a){
  df <- data.frame(a, 'mgs.tert' = as.factor(data.compare$MGS.tertile))
  df <- df[complete.cases(df),]
  mod <- lm(a~mgs.tert, data=df)
  mod.aov <- aov(mod)
  tuk <- TukeyHSD(mod.aov)
  pvals <- c(tuk$mgs.tert[,4])
  out <- c(summary(mod.aov)[[1]][['Pr(>F)']][1], as.vector(pvals)) ## order: mid-low, high-low, high-mid
  return(out)
})

output.MGStert <- do.call(rbind.data.frame, output.MGStert); names(output.MGStert) <- c('anova.p', 'mid.low', 'high.low', 'high.mid')
output.MGStert$anova.fdr <- p.adjust(output.MGStert$anova.p, method = 'fdr')
row.names(output.MGStert) <- names(met.clusters)

mean.per.mod <- met.clusters %>% 
  mutate('MGS.tertile' = data.compare$MGS.tertile) %>% 
  group_by(MGS.tertile) %>% 
  summarise(across(where(is.numeric), mean)) %>% 
  t() %>% 
  as.data.frame()

names(mean.per.mod) <- mean.per.mod[1,] ; mean.per.mod <- mean.per.mod[-1,]
row.names(mean.per.mod) == row.names(output.MGStert)

export.MGSrich.metmodules <- cbind(mean.per.mod, output.MGStert)
export.MGSrich.metmodules[export.MGSrich.metmodules$anova.fdr <.1, ]
write.table(export.MGSrich.metmodules, '/home/margar/download/export2_MGStertiles_vsMetClusters.txt', sep='\t')

## richness variables against delta values ----
head(data.compare)
delta.df <- readRDS('data/DeltaofLogValues.rds')
head(delta.df)
row.names(delta.df) <- delta.df$studyID
data.compare.prog <- data.compare[row.names(data.compare)%in%row.names(delta.df),]
dim(data.compare.prog)
identical(row.names(data.compare.prog), row.names(delta.df))
assoc.richvar.deltas <- cor.associations(delta.df, data.compare.prog[,1:3], names(delta.df)[6:24])
assoc.richvar.deltas <- do.call(rbind.data.frame, assoc.richvar.deltas)
names(assoc.richvar.deltas) <- c('biochem.delta', 'richness.type', 'rho', 'pval', 'fdr')
assoc.richvar.deltas[,3:5] <- sapply(assoc.richvar.deltas[, 3:5], as.numeric)
assoc.richvar.deltas %>% 
  filter(fdr<.1) %>% 
  ggplot(aes(biochem.delta, richness.type, fill = rho)) + 
    geom_tile()
write.table(assoc.richvar.deltas, '/home/margar/download/export_richvariables_deltavalues.txt', sep ='\t')
