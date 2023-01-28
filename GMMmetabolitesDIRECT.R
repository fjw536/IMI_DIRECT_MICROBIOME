## GMM metabolites corr associaitons

setwd('Downloads/data (3)/')

assoc <- read.table('cormatrixDIRECT_GMMfull.txt')
annot <- read.table('cormatrixDIRECT_GMMfull.txt default node.csv', header=T, sep=',')
annot[] <- sapply(annot, make.names)

assoc$fdr <- p.adjust(assoc$pvalue, method='fdr')

annot.to.add<- sapply(assoc$metabolite, function(o){
  a <- runif(1)
  ab <- annot[which(annot$name==make.names(o)),]
  return(paste(ab$CHEMICAL_NAME, ab$SUPER_PATHWAY, ab$SUB_PATHWAY, sep ="|"))
})

annot.to.add2 <- sapply(assoc$GMM, function(o){
  ab <- annot[which(annot$name==make.names(o)),]
  return(paste(ab$GMMname, ab$GMMhier1, ab$GMMhier2, sep ="|"))
})

assoc <- cbind(assoc, annot.to.add, annot.to.add2)
assoc$name <- unlist(strsplit(assoc$annot.to.add, '\\|'))[c(T,F,F)]
assoc$type <- unlist(strsplit(assoc$annot.to.add, '\\|'))[c(F,T,F)]
assoc$path <- unlist(strsplit(assoc$annot.to.add, '\\|'))[c(F,F,T)]
assoc$GMMname <- unlist(strsplit(assoc$annot.to.add2, '\\|'))[c(T,F,F)]
assoc$GMMhier1 <- unlist(strsplit(assoc$annot.to.add2, '\\|'))[c(F,T,F)]
assoc$GMMhier2 <- unlist(strsplit(assoc$annot.to.add2, '\\|'))[c(F,F,T)]

assoc$annot.to.add <- NULL
assoc$annot.to.add2 <- NULL

assoc %>% 
  group_by(GMMhier2, type) %>% 
  tally() %>% 
  write.table('GMMhier2mettype_network.txt', sep ='\t', row.names=F, quote = F)


  
ggplot(assoc, aes())

assoc %>% 
  group_by(GMMhier1, path) %>% 
  tally() #%>% 
  write.table('genusmepath_network.txt', sep ='\t', row.names=F, quote =F)
  
library(tidyverse)  
met.order <- arrange(assoc, type, path)
met.order2 <- factor(unique(met.order$path), levels=unique(met.order$path)) 
gmm.order <- arrange(assoc, GMMhier1, GMMhier2)  
gmm.order2 <- factor(unique(gmm.order$GMMname), levels = unique(gmm.order$GMMname))

assoc %>% 
  group_by(GMMname, path) %>% 
  summarise_at('corr', median) %>% 
  mutate(GMMname = factor(GMMname, levels=gmm.order2)) %>% 
  mutate(path = factor(path, levels=met.order2)) %>% 
  ggplot(aes(GMMname, path, fill=corr)) + geom_tile() +
    scale_fill_gradient2(low='blue', high='red') +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5, size =10),
          axis.text.y = element_text(size = 10)) +
  labs(x=NULL, y=NULL)# +
  ggsave(filename = 'HeatmapGMMMetsclusteredbypaths.pdf', width = 45, height = 35, units = 'cm', dpi = 300)
