
load('/home/scratch/vogt/IGC.keggOrtholog2geneIndex.RData')
load('/home/scratch/vogt/ModulDef_boolean_sum.RData')
load('/home/Data/Repository/Microbiome/WP/2018_09_11/output.1592samples.genecounts.RData')

#####
#KO annotation
#####

koCounts <- sapply(IGC.keggOrtholog2geneIndex, function(x){
    rowSums(output.1592samples.genecounts[, x, drop = FALSE])
})
save(koCounts, file = '/home/scratch/vogt/output.1592samples.KOcounts.RData')

#rarefy
library(GUniFrac)
koCounts.rare<-Rarefy(koCounts)$otu.tab.rff
save(koCounts.rare, file = '/home/scratch/vogt/output.1592samples.KOcounts.ds.RData')

##normalize
koCounts.rarenor <- t(apply(t(koCounts.rare), 2, function(i) i/sum(i)))

save(koCounts.rarenor, file = '/home/scratch/vogt/output.1592samples.KOabund.ds.RData')


#####
#KEGG modules
#####
koCounts_remove<-koCounts[, colnames(koCounts) %in% KEGGset]
#koCounts matrix from CM is missing some KOs (0 counts on all samples) -> find missing KOs
diff_set<-setdiff(KEGGset,colnames(koCounts_remove))

#make 0 value matrix with rownames of all samples and colnames missing KOs found previously + cbin (combine both dataframes)
add<-as.data.frame(matrix(0, ncol=length(diff_set), nrow = length(rownames(koCounts_remove))))
rownames(add)<-rownames(koCounts_remove)
colnames(add)<-diff_set
koCounts_all<-cbind(koCounts_remove,add)


#TRUE and FALSE dataframe of samples and if KO counts are above 0 
sampleHasKo <- as.data.frame(koCounts_all > 0)

#ModulDef_boolean_sum contains the Module KO composition inclding alternative pathways, the overall number of TRUE KOs found in TRUE and FALSE dataframe
ModuleStepCounts <- sapply(ModulDef_boolean_sum, eval, envir = sampleHasKo)
#percentage of actual module count and complete module count
ModuleCompletenessPct <- round(100 * t(t(ModuleStepCounts) / CompleteModuleCount))

#filter with 2/3 of KOs present (can be adjusted)
ModuleCompletenessPct_filter<-as.data.frame(ModuleCompletenessPct > 66)

#KOs abundances
Module_abundance <- sapply(ModulDef_boolean_sum, eval, envir = koCounts_all)


rownames(Module_abundance) <- rownames(koCounts)

#find abundance of KOs 
Module_abundance[ModuleCompletenessPct_filter == FALSE] <- 0


save(Module_abundance, file='/home/scratch/vogt/output.1592samples.KEGGmodulecounts.RData')


#rarefy module abundance
Module_abundance.rare<-Rarefy(Module_abundance)$otu.tab.rff

#normalize to module length
Module_abundance.rarelength <- round(t(t(Module_abundance) / CompleteModuleCount))

save(Module_abundance.rarelength, file='/home/scratch/vogt/output.1592samples.KEGGmodulecounts.ds.RData')


#normalize
Module_abundance.rarelength.nor <- t(apply(t(Module_abundance.rarelength), 2, function(i) i/sum(i)))


save(Module_abundance.rarelength.nor, file='/home/scratch/vogt/output.1592samples.KEGGmoduleabund.ds.RData')



######
#GMM
######

library(omixerRpm)
library(data.table)

tmp<-t(koCounts.rare)
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

save(abundancerownames2, file='/home/scratch/vogt/output.1592samples.GMMcounts.ds.RData')

##normalizing
abundancerownames2nor <- apply(abundancerownames2, 2, function(i) i/sum(i))

save(abundancerownames2nor, file='/home/scratch/vogt/output.1592samples.GMMabund.ds.RData')
 
 









