## 
## Blautia MGS genes 
##

# data
load('/home/scratch/UHRPP/shared_results/IGCMGS_July2018.RData')
blautia.mgs <- c('MGS:igc0074', 'MGS:igc0142', 'MGS:igc0160', 'MGS:igc0228', 'MGS:igc0242', 'MGS:igc0481',
                 'MGS:igc0508', 'MGS:igc0756', 'MGS:igc0779', 'MGS:igc1352', 'MGS:igc1419')

# extract genes
blautia.genes <- IGCMGS$sets[blautia.mgs]
Reduce(intersect, blautia.genes)
lapply(blautia.genes, length)
