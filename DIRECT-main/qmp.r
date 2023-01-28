
# this script doesn't include copy number correction, a function for copy number correction is included in RDP classifier 2.12 
# this script uses function rarefy_even_depth from phyloseq 1.20.0, it needs package phyloseq to be installed and loaded in order to work.
# with cnv_corrected_abundance_table: a copy number variation corrected abundance table with sample-identifiers as rows, copy number corrected taxa-abundances as columns
# with cell_counts_table: a table with sample-identifiers as rows, cell counts as columns 
library(phyloseq)
rarefy_even_sampling_depth <- function(cnv_corrected_abundance_table, cell_counts_table) 
{
  try(if(all(row.names(cnv_corrected_abundance_table) == row.names(cell_counts_table))==FALSE) stop("Cnv_corrected_abundance_table and cell_counts_table do not have the same sample-names, Please check!"))
  cnv_corrected_abundance_table = ceiling(cnv_corrected_abundance_table) # data values are rounded up in order to make use of integer values during the calculations
  cell_counts_table = t(cell_counts_table[order(row.names(cnv_corrected_abundance_table)),]) # make sure the order of the samples is the same in both files  
  sample_sizes = rowSums(cnv_corrected_abundance_table) # sample size of each sample (total nr of reads)
  sampling_depths = sample_sizes / cell_counts_table # sampling depth of each sample (total nr of reads divided by the cell count)
  minimum_sampling_depth = min(sampling_depths) # minimum of all sampling depths
  rarefy_to = cell_counts_table * minimum_sampling_depth # nr of reads to rarefy in each sample in order to get to an even sampling depth over all samples
  cnv_corrected_abundance_table_phyloseq = otu_table(cnv_corrected_abundance_table, taxa_are_rows = FALSE) # convert to phyloseq otutable
  rarefied_matrix=matrix(nrow = nrow(cnv_corrected_abundance_table_phyloseq), ncol = ncol(cnv_corrected_abundance_table_phyloseq), dimnames = list(rownames(cnv_corrected_abundance_table_phyloseq), colnames(cnv_corrected_abundance_table_phyloseq)))
  for (i in 1:nrow(cnv_corrected_abundance_table_phyloseq))
  {
    x <- rarefy_even_depth(cnv_corrected_abundance_table_phyloseq[i,], sample.size = rarefy_to[i], rngseed = 711, replace = FALSE, trimOTUs = F, verbose = FALSE)
    rarefied_matrix[i,] = x
  }
  normalised_rarefied_matrix = rarefied_matrix/rowSums(rarefied_matrix)
  QMP = normalised_rarefied_matrix*cell_counts_table[1,]
  return(QMP)
}

cell_counts <- read.table('/home/scratch/vogt/KO_annotation/QMP/WP2.1_Cell-Counts-NAR.csv', header=T, row.names = 1, sep='\t')
colnames(cell_counts) <- c('CC')

load('/home/scratch/vogt/KO_annotation/QMP/output.1592samples.mgscounts.RData')
load('/home/scratch/vogt/KO_annotation/QMP/output.1592samples.KOcounts.RData')
output.1592samples.KOcounts<-koCounts #rename KO object to fit the annotation
rm(koCounts)
load('/home/scratch/vogt/KO_annotation/QMP/output.1592samples.KEGGmodulecounts.RData')
output.1592samples.KEGGmodulecounts<-Module_abundance
rm(Module_abundance)


output.1592samples.mgscounts <- output.1592samples.mgscounts[is.element(rownames(output.1592samples.mgscounts), rownames(cell_counts)),]
output.1592samples.mgscounts <- output.1592samples.mgscounts[order(row.names(output.1592samples.mgscounts)),]

output.1592samples.KOcounts <- output.1592samples.KOcounts[is.element(rownames(output.1592samples.KOcounts), rownames(cell_counts)),]
output.1592samples.KOcounts <- output.1592samples.KOcounts[order(row.names(output.1592samples.KOcounts)),]

output.1592samples.KEGGmodulecounts <- output.1592samples.KEGGmodulecounts[is.element(rownames(output.1592samples.KEGGmodulecounts), rownames(cell_counts)),]
output.1592samples.KEGGmodulecounts <- output.1592samples.KEGGmodulecounts[order(row.names(output.1592samples.KEGGmodulecounts)),]


output.1592samples.mgscounts.qmp <- rarefy_even_sampling_depth(output.1592samples.mgscounts, cell_counts)
save(output.1592samples.mgscounts.qmp, file='/home/scratch/vogt/KO_annotation/QMP/output.1592samples.mgscounts.qmp.RData')

output.1592samples.KOcounts.qmp <- rarefy_even_sampling_depth(output.1592samples.KOcounts, cell_counts)
save(output.1592samples.KOcounts.qmp, file='/home/scratch/vogt/KO_annotation/QMP/output.1592samples.KOcounts.qmp.RData')

output.1592samples.KEGGmodulecounts.qmp <- rarefy_even_sampling_depth(output.1592samples.KEGGmodulecounts, cell_counts)
save(output.1592samples.KEGGmodulecounts.qmp, file='/home/scratch/vogt/KO_annotation/QMP/output.1592samples.KEGGmodulecounts.qmp.RData')

#######
##GMMs
#######

library(omixerRpm)
library(data.table)

tmp<-t(output.1592samples.KOcounts.qmp)
tmp<-as.data.frame(tmp)

setDT(tmp, keep.rownames = "entry")
#mods2 <- rpm(tmp, minimum.coverage=0.33, annotation = 1)
mods <- rpm(tmp, minimum.coverage=0.66, annotation = 1)
# Load the default mapping database
db <- loadDefaultDB()
 
# get the name of the first predicted module
getNames(db, mods@annotation[1,])
 
# get the abundance|coverage as a data.frame with module id and description
annotation<-as.data.frame(mods@annotation)
abundance<-as.data.frame(mods@abundance)
 
#binding and ready for normalising:
#cbind
abundancerownames<-cbind(annotation,abundance)
 
abundancerownames2 <- abundancerownames[,-1]
rownames(abundancerownames2) <- abundancerownames[,1]

abundancerownames2 <- t(abundancerownames2)

save(abundancerownames2, file='/home/scratch/vogt/KO_annotation/QMP/output.1592samples.GMMcounts.qmp.RData')





