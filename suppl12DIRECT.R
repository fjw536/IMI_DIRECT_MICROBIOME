tb12 <- read.table('~/Desktop/supplfig12forfigure.txt', head=T, sep='\t')

library(tidyverse)

tb12$bioclinical <- factor(tb12$bioclinical, 
                           levels = c("Waist (cm)",
                                      "Hip (cm)",            
                                      "Fasting plasma insulin (pmol/L)" ,
                                      "Basal insulin secretion (nmol m-2)",
                                      "OGTT insulin secretion (nmol m-2)",
                                      "OGTT insulin sensitivity (ml min-1m-2)",
                                      "Glucose sensitivity (pmol min-1 m-2 (mmol/L)-1)",
                                      "Matsuda (umol min-1kg-1)",
                                      "Fasting plasma GLP1 (pg/ml)",
                                      "Total plasma GLP1 (pg/ml)",
                                      "Fasting plasma hsCRP (mg/L)",
                                      "Fasting plasma HDL (mmol/L)",
                                      "Fasting plasma LDL  (mmol/L)",
                                      "Fasting plasma Triglycerides  (mmol/L)",
                                      "Fasting plasma cholesterol (mmol/L)",
                                      "Fasting plasma aspartate aminotransferase (U/L)", 
                                      "Liver iron (umol/g)",
                                      "Pancreatic iron (umol/g)",
                                      "Intra-Abdominal Adipose Tissue (L)",
                                      "Abdominal Subcutaneous Adipose Tissue (L)",
                                      "Total Abdominal Adiposse Tissue (L)",
                                      "Liver fat (%)",
                                      "Pancreas fat (%)",
                                      "Physical activity",
                                      "Smoking",
                                      "Healthy Dietary Index"))
                                                   

tb12 <- tb12[order(tb12$bioclinical, tb12$delta),]
tb12$Annotation <- factor(tb12$Annotation, levels = unique(tb12$Annotation))
pdf('phenoMGS.pdf', width = 25, height=10)
tb12 %>% 
  filter(FDR<=.1) %>% 
  filter(Feature=='MGS') %>%
  mutate('bioclinical' = droplevels(bioclinical)) %>% 
  ggplot(aes(Annotation, bioclinical, fill=delta)) +
    geom_tile() +
    scale_fill_gradient2(low='darkblue', mid='white', high='darkred') + 
    scale_y_discrete(limits = rev(levels(tb12$bioclinical))) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5),
        axis.text.y = element_text(size=8)) +
  labs(x = 'Metagenomic species (MGS) annotation', y = 'Phenotypical feature',
       fill = "Cliff's Delta effect size")
dev.off()

tb12 %>% 
  filter(FDR<=.1) %>% 
  filter(Feature=='Species') %>%
  ggplot(aes(bioclinical, Annotation, fill=delta)) +
  geom_tile() +
  scale_fill_gradient2(low='darkblue', mid='white', high='darkred') + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))

tb12 %>% 
  filter(FDR<=.1) %>% 
  filter(Feature=='Genus') %>%
  ggplot(aes(bioclinical, Annotation, fill=delta)) +
  geom_tile() +
  scale_fill_gradient2(low='darkblue', mid='white', high='darkred') + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5))

pdf('phenoGMM.pdf', height = 10, width=8)
tb12 %>% 
  filter(FDR<=.1) %>% 
  filter(Feature=='GMM') %>%
  ggplot(aes(bioclinical, Annotation, fill=delta)) +
  geom_tile() +
  scale_fill_gradient2(low='darkblue', mid='white', high='darkred') + 
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5)) +
  labs(y = 'Gut Metabolic Module (GMM) annotation', x = 'Phenotypical feature',
       fill = "Cliff's Delta effect size")
dev.off()
