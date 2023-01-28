setwd('~/Downloads/')

bsl.heat <- readRDS('phages_bslcor.rds')
tree <- ape::read.tree('phagestree.txt')

library(ggtree); library(tidyverse)

tree.df <- fortify(tree)
tree.df <- subset(tree.df, isTip)
heat.order <- factor(tree.df$label, levels = tree.df$label[order(tree.df$y, decreasing = F)])

pdf('~/Desktop/DIRECT_phagesheatmap.pdf', height=15, width = 15, family='Times')
bsl.heat.plot <- bsl.heat %>% 
  filter(fdr<.1) %>% 
  mutate('association' = gsub('\\.[0-9]+', '', row.names(.))) %>% 
  mutate('association'= factor(association, levels=unique(association))) %>% 
  mutate('association' = plyr::revalue(association,
                                       c('Height' = 'Height (cm)',
                                         'Weight.kg' = 'Weight (kg)',
                                         'Waist.cm' = 'Waist (cm)', 
                                         'Hip.cm' =  'Hip (cm)',
                                         'Waist.Hip' = 'Waist:Hip',
                                         'BMI' = 'Body Mass Index (kg/m2)',
                                         'Fasting.Glucose' = 'Fasting P. Glucose (mmol/L)',
                                         'Basal.Glucose' = 'Basal P. Glucose (mmol/L)',
                                         'Mean.Glucose' = 'Mean OGTT P. Glucose (mmol/L)',
                                         'glucoseh' = 'OGTT P. Glucose 2-hours (mmol/L)',
                                         'Fast.Insulin' = 'Fasting P. Insulin (pmol/L)',
                                         'Basal.Insulin' = 'Basal P. Insulin (pmol/L)',
                                         'Mean.Insulin' = 'OGTT P. Insulin (pmol/L)',
                                         'Basal.InsulinSecretion' = 'Basal insulin secretion (pmol min-1m-2)',
                                         'Total.InsulinSecretion' = 'OGTT insulin secretion (nmol m-2)',
                                         'Glucose.Sensitivity' = 'Glucose sensitivity (pmol min-1 m-2 (mmol/L)-1)',
                                         'OGISh'= 'OGTT insulin sensitivity (ml min-1m-2)',
                                         'Stumvoll' = "Stumvoll's index (umol min-1kg-1)",
                                         'Matsuda' = "Matsuda's index (umol min-1kg-1)",
                                         'Clinsb' = 'Basal insulin clearance (L min-1m-2)',
                                         'Clins' = 'OGTT insulin clearance (L min-1m-2)',
                                         'Active.GLP1.Concm' = "Fasting P. GLP1 (pg/ml)",
                                         'Total.GLP1.Concm' = 'Total P. OGTT GLP1 (pg/ml)',
                                         'hsCRP' = "Fasting P. hsCRP (mg/L)",
                                         'Fasting.HDL' = "Fasting P. HDL  (mmol/L)",
                                         'Fasting.LDL' = "Fasting P. LDL  (mmol/L)",
                                         'Fasting.TG' = "Fasting P. Triglycerides  (mmol/L)",
                                         'Fasting.ALT' = "Fasting P. ALT (U/L)",
                                         'Fasting.AST' = "Fasting P. AST (U/L)",
                                         'Fasting.Chol' = "Fasting P. Cholesterol  (mmol/L)",
                                         'Liver.Iron' = "Liver iron (umol/g)",
                                         'Panc.Iron' = "Pancreatic iron (umol/g)",
                                         'IAAT' = "IAAT (L)",
                                         'ASAT' = "ASAT (L)",
                                         'TAAT' = "TAAT (L)",
                                         'LiverFat..' = "Liver fat (%)",
                                         'PancreasFat..' = "Pancreatic fat (%)",
                                         'BP.Sys' = "Mean Systolic B.P. (mmHg)",
                                         'BP.D.Mean' = "Mean Diastolic B.P. (mmHg)",
                                         'Value.vm.hpf.mean.BL' = "Physical activity",
                                         'HDI' = "Healthy Dietary Index (HDI)"))) %>% 
  mutate('rho' = as.numeric(rho)) %>% 
  mutate('metabolite' = factor(metabolite, levels=levels(heat.order))) %>%
  as.data.frame()
  
ggplot(bsl.heat.plot, aes(association, metabolite, fill=rho)) +
  geom_tile() + scale_fill_gradient2(low='blue', mid='white', high='red') +
  scale_y_discrete(limits=rev(levels(bsl.heat.plot$metabolite))) +
  blank_theme() + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5)) + 
  labs(x=NULL, y= NULL, fill = "rho")

dev.off()
pdf('~/Desktop/DIRECT_phagestree.pdf', height=15, width = 7, family='Times')
ggtree(tree, branch.length = 'none') + geom_tiplab()
dev.off()
