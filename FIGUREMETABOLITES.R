###########################################################################################################
############################ Metabolites-Microbiome summary figure ########################################
###########################################################################################################

# directory
setwd('~/Desktop')

# libraries
library(xlsx); library(ggplot2); library(tidyverse)

# data
gmm.mets <- read.xlsx('DIRECT_MetClustersFigures2.xlsx', sheetName = 'metclustersgmm', header=T)
gc()
lassocoeffs <- read.xlsx('DIRECT_MetClustersFigures2.xlsx', sheetName = 'lassocoefs', header=T)
gc()
tax <- read.xlsx('DIRECT_MetClustersFigures2.xlsx', sheetName = 'tax', header=T); tax$mgs <- gsub(':', '\\.', tax$mgs)
gc()
metclusters.annot <- read.xlsx('DIRECT_MetClustersFigures2.xlsx', sheetName = 'metclustersannot', header=T)
gc()
bact2 <-  read.xlsx('DIRECT_MetClustersFigures2.xlsx', sheetName = 'b2', header=T)
gc()
gr.tert <-  read.xlsx('DIRECT_MetClustersFigures2.xlsx', sheetName = 'GRtertiles', header=T)
gc()
gr.cont <-  read.xlsx('DIRECT_MetClustersFigures2.xlsx', sheetName = 'grcontinous', header=T)
gc()
variance <-  read.xlsx('DIRECT_MetClustersFigures2.xlsx', sheetName = 'lasso', header=T)
variance <- variance[!is.na(variance$module.name),]
gc()
t2d <-  read.xlsx('DIRECT_MetClustersFigures2.xlsx', sheetName = 't2d', header=T)
gc()
ngr<-  read.xlsx('DIRECT_MetClustersFigures2.xlsx', sheetName = 'ngr', header=T)
gc()

## metabolites explained by bacteria
# lassocoeffs$met.name <- sapply(lassocoeffs$met, function(a){
#   metclusters.annot[grep(paste0('^', a, '$'), metclusters.annot$cluster),]$name
# })
met.order <- sort(table(lassocoeffs$module.name), decreasing = F)
lassocoeffs$module.name <- factor(lassocoeffs$module.name, levels = names(met.order))
lassocoeffs$mgs.name <- sapply(lassocoeffs$mgs, function(a){
  tax[grep(paste0('^', a, '$'), tax$mgs),]$name
})

tax.met <- tax[tax$mgs%in%lassocoeffs$mgs,]
lassocoeffs$mgs <- factor(lassocoeffs$mgs,levels = tax.met$mgs)
lassocoeffs$mgs.name <- factor(lassocoeffs$mgs.name, levels = unique(tax.met$name))
lassocoeffs$mgs.name <-  factor(lassocoeffs$mgs.name, levels = rev(levels(lassocoeffs$mgs.name)))
lassocoeffs2 <- lassocoeffs[lassocoeffs$module.name%in%names(tail(met.order, n=26)),]
pdf('mgsandgmmassociationsmets.pdf', width = 20, height = 15)
cowplot::plot_grid(ggplot(lassocoeffs2, aes(mgs.name, module.name, fill='a')) +
            geom_tile(color = 'black', size =0.1) + 
            scale_fill_manual(values = c('darkred')) +
            guides(fill=F) +# coord_flip() +
            theme(panel.background = element_blank(), panel.border = element_rect(fill=NA), 
                  panel.grid.major = element_line(color='gray90'), 
                  axis.text.x = element_text(angle = 90, vjust = .9, hjust = 1)) +
            labs(y='Metabolite clusters', x = 'Bacterial species (MGS)') ,
          gmm.mets %>% 
            filter(fdr<=.1) %>% 
            mutate('value' = as.numeric(value)) %>% 
            ggplot(aes(gmm.name, module.name, fill = value)) + 
            geom_tile(color = 'black', size =0.1) +
            scale_fill_gradient2(low='darkblue', high='darkred') +
            guides(fill=F) + #coord_flip() +
            theme(panel.background = element_blank(), panel.border = element_rect(fill=NA), 
                  panel.grid.major = element_line(color='gray90'), 
                  axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1)) +
            labs(y='Metabolite clusters', x = 'Gut Metabolic Modules'),
          rel_heights  = c(1.3,1), ncol = 1)
dev.off()

ggplot(lassocoeffs2, aes(mgs.name, module.name, fill='a')) +
  geom_tile(color = 'black', size =0.1) + 
  scale_fill_manual(values = c('darkred')) +
  guides(fill=F) + coord_flip() +
  theme(panel.background = element_blank(), panel.border = element_rect(fill=NA), 
        panel.grid.major = element_line(color='gray90'), 
        axis.text.x = element_text(angle = 90, vjust = .9, hjust = 1)) +
  labs(y='Metabolite clusters', x = 'Bacterial species (MGS)')

gmm.mets %>% 
  filter(fdr<=.1) %>% 
  mutate('value' = as.numeric(value)) %>% 
  ggplot(aes(gmm.name, module.name, fill = value)) + 
  geom_tile(color = 'black', size =0.1) +
  scale_fill_gradient2(low='darkblue', high='darkred') +coord_flip() +
  guides(fill=F) +
  theme(panel.background = element_blank(), panel.border = element_rect(fill=NA), 
        panel.grid.major = element_line(color='gray90'), 
        axis.text.x = element_text(angle = 90, vjust = .9, hjust = 1, size = 11),
        axis.text.y = element_text(size=11),
        axis.title = element_text(size=13)) +
  labs(y='Metabolite clusters', x = 'Gut Metabolic Modules (GMM)')
 
## circular heatmap
# prepare data
# compute difference of means and retrive only significant
bact2$diff.means <- (bact2$noBact2+1) - (bact2$Bact2+1)
out <- c()
for(i in seq(1:nrow(bact2))){
  if (bact2[i,]$fdr <= 0.1){
    out[i] <- bact2[i,]$diff.means
  } else (
    out[i] <- 0
  )
}
bact2$diff.means2 <- out

gr.tert$diff.means <- (gr.tert$EGC+1) - (gr.tert$LGC+1)
out <- c()
for(i in seq(1:nrow(gr.tert))){
  if (gr.tert[i,]$fdr <= 0.1){
    out[i] <- gr.tert[i,]$diff.means
  } else (
    out[i] <- 0
  )
}
gr.tert$diff.means2 <- out

t2d$diff.means <- (t2d$T2D+1) - (t2d$mean.noT2D+1)
out <- c()
for(i in seq(1:nrow(t2d))){
  if (t2d[i,]$FDR <= 0.1){
    out[i] <- t2d[i,]$diff.means
  } else (
    out[i] <- 0
  )
}
t2d$diff.means2 <- out

ngr$diff.means <- (ngr$mean.NGR+1) - (ngr$mean.noNGR+1)
out <- c()
for(i in seq(1:nrow(ngr))){
  if (ngr[i,]$FDR <= 0.1){
    out[i] <- ngr[i,]$diff.means
  } else (
    out[i] <- 0
  )
}
ngr$diff.means2 <- out

# correlations significant
gr.cont$rho2 <- ifelse(gr.cont$fdr<=0.1, gr.cont$rho, 0)
# check row names
bact2$module.name == t2d$annotation
t2d$annotation == ngr$annotation
unique(gr.cont$metabolite2) == t2d$annotation

heatmap.discr <- data.frame(
  't2d.means' = t2d$diff.means2,
  'ngr.means' = ngr$diff.means2,
  'gr.tert.means' = gr.tert$diff.means2,
  'b2.means' = bact2$diff.means2
)
row.names(heatmap.discr) <- bact2$module.name
# row.names(heatmap.discr) <- metclusters.annot$cluster

heatmap.cont <- data.frame(
  'gr.cont' = gr.cont[gr.cont$variable=='gene.rich',]$rho2,
  'mgs.rich' = gr.cont[gr.cont$variable=='MGS.rich',]$rho2,
  'shannon' = gr.cont[gr.cont$variable=='shannon',]$rho2
)
row.names(heatmap.cont) <- t2d$annotation

variance.circos <- variance$var; names(variance.circos) <- variance$module.name

heatmap.all <- cbind(heatmap.cont, heatmap.discr, variance.circos); row.names(heatmap.all) <- row.names(heatmap.discr)
cluster.keep <- row.names(heatmap.all[!rowSums(heatmap.all)==0,])

heatmap.discr <- heatmap.discr[row.names(heatmap.discr)%in%cluster.keep,]
heatmap.cont <- heatmap.cont[row.names(heatmap.cont)%in%cluster.keep,]
names(variance.circos) <- row.names(heatmap.all)
variance.circos <- variance.circos[names(variance.circos)%in%cluster.keep]
t2d.means <- heatmap.discr$t2d.means; names(t2d.means) <- row.names(heatmap.discr); t2d.means <- as.data.frame(t2d.means)

library(circlize)
col_fun.t2d <- colorRamp2(c(-0.05, -0.025, 0, 0.025, 0.05), c("darkgreen", "forestgreen", "white", "darkorchid2", "darkorchid4"))
col_fun.ngr <- colorRamp2(c(-0.05, -0.025, 0, 0.025, 0.05), c("darkorchid4", "darkorchid2", "white",  "forestgreen","darkgreen"))
col_fun.grT <- colorRamp2(c(-0.05, -0.025, -0.01, -0.005, 0, 0.025, 0.05), c("gold3", 'gold2', "yellow2", 'yellow1', "white", "dodgerblue1", "blue1"))
col_fun.B2 <- colorRamp2(c(-0.05, -0.025, 0, 0.025, 0.05), c('red4', "red2", "white", "forestgreen", "darkgreen"))
col_fun.gr.cont <- colorRamp2(c(-0.4, -0.2, 0, 0.2, 0.4), c("blue1", "dodgerblue1", "white", "red2", "red4"))
col_fun.var <- colorRamp2(c(0, 0.1, 0.2, 0.3, 0.4, 0.5), c("white", "red", "red1", "red2", "red3", "darkred"))

pdf('circularheatmap.pdf', height = 30, width = 30)
circos.clear()
circos.par(start.degree = 25, gap.degree = 68)
circos.heatmap(t2d.means, col = col_fun.t2d, rownames.side = 'outside', rownames.cex = 1.5, cluster = T, track.height =.03, bg.border = 'black') #replace cex with 1.5
circos.heatmap(heatmap.discr$ngr.means, col = col_fun.ngr, track.height =.03,  bg.border = 'black')
circos.heatmap(heatmap.discr$b2.means, col = col_fun.B2, track.height =.03,  bg.border = 'black')
circos.heatmap(heatmap.discr$gr.tert.means, col = col_fun.grT, track.height =.03,  bg.border = 'black')
circos.heatmap(heatmap.cont, col = col_fun.gr.cont, track.height =0.09,  bg.border = 'black')
circos.heatmap(variance.circos, col = col_fun.var, track.height =.03,  bg.border = 'black')
dev.off()
