# deconf <- readRDS('~/Downloads/deconfoundedfeatures')
deconf <- readRDS('~/Downloads/drugdeconf.rds')
  
deconf.drugs <- deconf[,5:23]
deconf.drugs.df <- as.data.frame(sapply(deconf.drugs, function(x){as.factor(x)}))
deconf.drugs[1:10,1:10]
deconf.drugs$type <- c(rep('MGS', 721), rep('GMM', 97), rep('phages', 347),
                       rep('AMR_drugclass', 33), rep('AMR_geneF', 81),
                       rep('SV_var',352), rep('SV_del', 986),
                       rep('metabolites', 369))


row.names(deconf.drugs[grep('_C', deconf.drugs$MGS.rich),])
confounded.feats <- sapply(deconf.drugs, function(a){
  paste0(row.names(deconf.drugs[grep('_C', a),]), '|', deconf.drugs[grep('_C', a), ]$type)
})

confounded.feats <- unique(unlist(confounded.feats))
write.table(confounded.feats, 'confoundeddeaturestotal.txt')


deconf.drugs2 <- deconf.drugs %>% 
  mutate_all(funs(str_replace(., "SD", "1"))) %>% 
  mutate_all(funs(str_replace(., "NS", "0"))) %>% 
  mutate_all(funs(str_replace(., "LD", "0.5"))) %>% 
  as.data.frame()

row.names(deconf.drugs2) <- row.names(deconf.drugs)
deconf.drugs2[,1:29] <- sapply(deconf.drugs2[,1:29], as.numeric)
deconf.drugs2$feature <- row.names(deconf.drugs2)

theme_clean <- function(){
  theme(panel.background = element_blank(), panel.border = element_rect(fill=NA), 
        panel.grid.major = element_line(color='gray90'), 
        axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1))
  }

deconf.drugs2 %>% 
  filter(type=='MGS') %>% 
  reshape2::melt() %>% 
  # mutate('value' = replace_na(value, 3)) %>% 
  filter(value!=0) %>% 
  n_distinct('feature')

deconf.drugs2 %>% 
  filter(type=='phages') %>% 
  reshape2::melt() %>% 
  # mutate('value' = replace_na(value, 3)) %>% 
  filter(value!=0) %>% 
  n_distinct('feature')

pdf('~/Desktop/DrugDeconfounding_summary.pdf', width = 15, height = 8)
deconf.drugs2 %>% 
  filter(type=='MGS') %>% 
  reshape2::melt() %>% 
  # mutate('value' = replace_na(value, 3)) %>% 
  filter(value!=0) %>% 
  mutate('value' = factor(value)) %>% 
  mutate('value' =plyr::revalue(value, c('1' = 'SD', '3' = 'mixed'))) %>% 
  ggplot(aes(variable, feature, fill=value)) + geom_tile() +
    scale_fill_manual(values=c('darkred', 'orange')) +
    coord_flip() +
    theme_clean() +
    labs(y= "MGS", x = 'Drugs', fill = 'Deconfounded by')

deconf.drugs2 %>% 
  filter(type=='GMM') %>% 
  reshape2::melt() %>% 
  # mutate('value' = replace_na(value, 3)) %>% 
  filter(value!=0) %>% 
  mutate('value' = factor(value)) %>% 
  mutate('value' =plyr::revalue(value, c('1' = 'SD', '3' = 'mixed'))) %>% 
  ggplot(aes(variable, feature, fill=value)) + geom_tile() +
    scale_fill_manual(values=c('darkred', 'orange')) +
    coord_flip() +
    theme_clean() +
    labs(y= "GMMs", x = 'Drugs', fill = 'Deconfounded by')

deconf.drugs2 %>% 
  filter(type=='phages') %>% 
  reshape2::melt() %>% 
  # mutate('value' = replace_na(value, 3)) %>% 
  filter(value!=0) %>% 
  mutate('value' = factor(value)) %>% 
  mutate('value' =plyr::revalue(value, c('1' = 'SD', '3' = 'mixed'))) %>% 
  ggplot(aes(variable, feature, fill=value)) + geom_tile() +
    scale_fill_manual(values=c('darkred', 'orange')) +
    coord_flip() +
    theme_clean() +
    labs(y= "Phages", x = 'Drugs', fill = 'Deconfounded by')

deconf.drugs2 %>% 
  filter(type=='AMR_drugclass') %>% 
  reshape2::melt() %>% 
  # mutate('value' = replace_na(value, 3)) %>% 
  filter(value!=0) %>% 
  mutate('value' = factor(value)) %>% 
  mutate('value' =plyr::revalue(value, c('1' = 'SD', '3' = 'mixed'))) %>% 
  ggplot(aes(variable, feature, fill=value)) + geom_tile() +
    scale_fill_manual(values=c('darkred', 'orange')) +
    coord_flip() +
    theme_clean() +
    labs(y= "AMR drug classes", x = 'Drugs', fill = 'Deconfounded by')

deconf.drugs2 %>% 
  filter(type=='AMR_geneF') %>% 
  reshape2::melt() %>% 
  # mutate('value' = replace_na(value, 3)) %>% 
  filter(value!=0) %>% 
  mutate('value' = factor(value)) %>% 
  mutate('value' =plyr::revalue(value, c('1' = 'SD', '3' = 'mixed'))) %>% 
  ggplot(aes(variable, feature, fill=value)) + geom_tile() +
    scale_fill_manual(values=c('darkred', 'orange')) +
    coord_flip() +
    theme_clean() +
    labs(y= "AMR gene families", x= 'Drugs', fill = 'Deconfounded by')

deconf.drugs2 %>% 
  filter(type=='metabolites') %>% 
  reshape2::melt() %>% 
  # mutate('value' = replace_na(value, 3)) %>% 
  filter(value!=0) %>% 
  mutate('value' = factor(value)) %>% 
  mutate('value' =plyr::revalue(value, c('1' = 'SD', '3' = 'mixed'))) %>% 
  ggplot(aes(variable, feature, fill=value)) + geom_tile() +
    scale_fill_manual(values=c('darkred', 'orange')) +
    coord_flip() +
    theme_clean() +
    labs(y= "Metabolite clusters and individual metabolites", x = 'Drugs', fill = 'Deconfounded by')

svs.deconf <- rbind(
  deconf.drugs2 %>% 
  filter(type=='SV_del') %>% 
  reshape2::melt() ,
  deconf.drugs2 %>% 
    filter(type=='SV_var') %>% 
    reshape2::melt())

svs <- strsplit(svs.deconf$feature, '_')
svs <- unlist(lapply(svs, function(x)(x[1])))
svs.deconf$feature <- svs

svs.deconf %>% 
  filter(type=='SV_var') %>% 
  # mutate('value' = replace_na(value, 3)) %>% 
  filter(value!=0) %>% 
  mutate('value' = factor(value)) %>% 
  mutate('value' =plyr::revalue(value, c('1' = 'SD', '3' = 'mixed'))) %>% 
  ggplot(aes(variable, feature, fill=value)) + geom_tile() +
    scale_fill_manual(values=c('darkred', 'orange')) +
    coord_flip() +
    theme_clean() +
    labs(y= "SVs variable regions", x = 'Drugs', fill = 'Deconfounded by')

svs.deconf %>% 
  filter(type=='SV_del') %>% 
  # mutate('value' = replace_na(value, 3)) %>% 
  filter(value!=0) %>% 
  mutate('value' = factor(value)) %>% 
  mutate('value' =plyr::revalue(value, c('1' = 'SD', '3' = 'mixed'))) %>% 
  ggplot(aes(variable, feature, fill=value)) + geom_tile() +
  scale_fill_manual(values=c('darkred', 'orange')) +
  coord_flip() +
  theme_clean() +
  labs(y= "SVs deletions", x = 'Drugs', fill = 'Deconfounded by')
dev.off()


phages.status <- as.data.frame(readRDS('~/Downloads/phages_status.rds'))
phages.status <- phages.status[,24:52]
phages.status2 <- phages.status %>% 
  mutate_all(funs(str_replace(., "SD", "1"))) %>% 
  mutate_all(funs(str_replace(., "NS", "0"))) %>% 
  mutate_all(funs(str_replace(., "LD", "0.5"))) %>% 
  as.data.frame()

row.names(phages.status2) <- row.names(phages.status)
phages.status2[,1:29] <- sapply(phages.status2[,1:29], as.numeric)
phages.status2$feature <- row.names(phages.status2)

theme_clean <- function(){
  theme(panel.background = element_blank(), panel.border = element_rect(fill=NA), 
        panel.grid.major = element_line(color='gray90'), 
        axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1))
}

phages.status2 %>% 
  reshape2::melt() %>% 
  mutate('value' = replace_na(value, 3)) %>% 
  filter(value!=0) %>% 
  mutate('value' = factor(value)) %>% 
  mutate('value' =plyr::revalue(value, c('1' = 'SD', '3' = 'mixed'))) %>% 
  ggplot(aes(variable, feature, fill=value)) + geom_tile() +
  scale_fill_manual(values=c('red', 'darkred', 'orange')) +
  coord_flip() +
  theme_clean() +
  labs(y= "Phages", x = 'Drugs', fill = 'Deconfounded by')


ph.value <- phages.status2 %>% 
  reshape2::melt() %>% 
  mutate('value' = replace_na(value, 3)) %>% 
  filter(value!=0) %>% 
  mutate('value' = factor(value)) %>% 
  select(value) 


rf.model <- read.table('~/Desktop/DIRECT_drugdeconf/listoffeaturesmodel.txt')
rf.model <- make.names(as.vector(rf.model$V1))
deconf.features <- deconf.drugs2 %>% 
  reshape2::melt() %>% 
  filter(value==1) %>% 
  select(feature) %>% 
  as.vector() %>% 
  unique()

deconf.features <- as.vector(deconf.features$feature)
setdiff(deconf.features, rf.model)
intersect(deconf.features, rf.model)
