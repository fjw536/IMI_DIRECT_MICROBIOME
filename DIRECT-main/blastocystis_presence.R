############################################################################################################
#################################### Blastocystis screen-out################################################
############################################################################################################

# directory
setwd('/home/scratch/margar')

# libraries
library(phyloseq); library(tidyverse); library(corrplot); library(ggpmisc); library(ggpubr); library(ggbeeswarm)

## Baseline ----
# data upload
qmp <- read.table('data/WP2.1_775_IGR_MGS-QMP_igrHA_Ent_M.csv', header = T, row.names=1)
qmp <- qmp[,-c(1:6)]
names(qmp) <- gsub('\\.', ':', names(qmp))
  
tax <- read.table('data/taxonomy.tsv', header = T, sep='\t')
blastos <- tax[grep('Blastocystis', tax$genus),]$X

qmp.blast <- qmp[,names(qmp)%in%blastos]
## it's not in the most abundant MGS

# work with the full dataset
load('data/output.1592samples.dsmgscounts.RData')
qmp <- data.frame(output.1592samples.dsmgscounts)
rm(output.1592samples.dsmgscounts)

names(qmp) <- gsub('\\.', ':', names(qmp))
qmp.blast <- qmp[,names(qmp)%in%blastos]
row.names(qmp.blast) <- gsub('X', '', row.names(qmp.blast))
head(qmp.blast)

mtdt.bsl <- read.table('data/curatedmetadata/WP2.1_BL_MetadataWithInfl.txt', header=T, row.names=1, sep='\t')
qmp.blast <- qmp.blast[row.names(qmp.blast)%in%mtdt.bsl$SampleID,]
qmp.blast <- qmp.blast[match(mtdt.bsl$SampleID, row.names(qmp.blast)),]

## stratification on phenotype (IGR type)
qmp.blast <- cbind(qmp.blast, 'GlyCat'=mtdt.bsl$Gly.Cat2)
qmp.blast %>% 
  group_by(GlyCat) %>% 
  summarise_each(funs(sum(.!=0)))

qmp.blast.m <- reshape2::melt(qmp.blast)
ggplot(qmp.blast.m, aes(value)) + geom_histogram() + facet_wrap(~variable, scales='free')

## differences in phenotype
mtdt.bsl.num <- select_if(mtdt.bsl, is.numeric)
qmp.blast <- cbind(qmp.blast, mtdt.bsl.num)
qmp.blast <- qmp.blast[,-c(4:7)]
qmp.blast <- qmp.blast[complete.cases(qmp.blast),]

qmp.blast.cor <- cor(qmp.blast)
qmp.blast.cor.p <- cor.mtest(qmp.blast)
corrplot(qmp.blast.cor[4:63,1:3], p.mat = qmp.blast.cor.p$p[4:63,1:3], insig = 'blank', method = 'square')

qmp.blast.red <- qmp.blast[,names(qmp.blast)%in%c("MGS:igc0001", "MGS:igc0003", "MGS:igc0004", "DSGeneAbundRichness",
                                                  "DSMGSAbundShannon", "Waist.cm", "IAAT", "BP.S.Mean")]
names(qmp.blast.red) <- gsub(':', '\\.', names(qmp.blast.red))
qmp.blast.red.m <- reshape2::melt(qmp.blast.red, id.vars=c("MGS.igc0001", "MGS.igc0003", "MGS.igc0004"))
qmp.blast.red.m <- data.frame('mgs'= c(rep('MGS.igc0001', nrow(qmp.blast.red.m)), rep('MGS.igc0003', nrow(qmp.blast.red.m)), rep('MGS.igc0004', nrow(qmp.blast.red.m))),
                              'counts'= c(qmp.blast.red.m$MGS.igc0001, qmp.blast.red.m$MGS.igc0003, qmp.blast.red.m$MGS.igc0004),
                              'variable' = rep(qmp.blast.red.m$variable, 3),
                              'value'=rep(qmp.blast.red.m$value, 3))

ggplot(qmp.blast.red.m, aes(counts, value, color=mgs)) + geom_point() + 
  geom_smooth(method = 'lm') + facet_wrap(~variable, scales = 'free') +
  stat_fit_glance(method='lm', geom='text', aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")),size = 3, label.x = 1500)

## mtdt progression
mtdt.fw <- read.table('data/curatedmetadata/WP2.1_BlFw-641_MetadataImputed_OnlyBIOCHEM_2.txt', sep ='\t', header=T)
mtdt.fw <- mtdt.fw[mtdt.fw$timepoint=='Follow-up',]
qmp.blast <- qmp.blast[row.names(qmp.blast)%in%mtdt.fw$SampleID_baseline,]

blast.df <- data.frame('Blastocystis' = rowSums(qmp.blast))
mtdt.fw <- mtdt.fw[match(row.names(qmp.blast), mtdt.fw$SampleID_baseline),]
blast.df$GlyCat <- mtdt.fw$Gly.Cat
ggplot(blast.df, aes(GlyCat, Blastocystis)) + geom_boxplot() + geom_quasirandom()

blast.df.g <- split(blast.df, blast.df$GlyCat)
lapply(blast.df.g, function(a){
  colSums(a==0)
})

blast.df.prev <- data.frame('class' = c('NGR', 'IGR', 'T2D'), 
                            'present' = c(31/105, 119/460, 22/76),
                            'non-present' = c(74/105, 341/460, 54/76))
chisq.test(blast.df.prev[,2:3])
