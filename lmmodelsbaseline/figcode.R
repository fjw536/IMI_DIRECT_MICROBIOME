## models GR and enteros - correct data
setwd('~/Desktop/01_Papers/04_DIRECT/stuff/lmmodelsbaseline/')

library(ggplot2); library(grid)
models <- read.table('LM_GRtertiles.txt', header = T, row.names=1,  sep = '\t')
models$variables <- factor(row.names(models), levels = row.names(models))

ggplot(models, aes(delta, variables, color = ifelse(fdr<=0.1, 'red', 'black'))) + 
  annotation_custom(grobTree(textGrob('High gene richness', x = .5, y = .9, gp=gpar(col = 'gray', fontsize = 25, fontface = 'bold.italic')))) +
  annotation_custom(grobTree(textGrob('Low gene richness', x = .5, y = .1, gp=gpar(col = 'gray', fontsize = 25, fontface = 'bold.italic')))) +
  geom_point() + 
  geom_segment(aes(x=delta, xend =0, y=variables, yend=variables)) + coord_flip() +
  scale_color_manual(values = c('black', 'red'), labels = c('Non-significant', 'Significant')) + xlim(-.5, .3) + geom_vline(xintercept = 0, lty = 'dashed') + 
  theme(axis.line = element_line(linetype = "solid"), panel.grid.major = element_line(colour = "gray78", linetype = "dashed"), 
        axis.title = element_text(size = 15), 
        axis.text.y = element_text(colour = "black"), 
        axis.text.x = element_text(colour = "black", angle=90, hjust=1, vjust=.5), 
        plot.title = element_text(size = 17),
        legend.title = element_blank(), panel.background = element_rect(fill = NA), legend.key = element_rect(fill = NA), 
        legend.background = element_rect(fill = NA), panel.border = element_rect(fill=NA)) +
  labs(title = "Lifestyle factors association with highest and lowest tertile of microbial gene richness", x = "Delta value", y = NULL,caption = "Significant q-value (BH<0.1) is showed in red") 

models$variables <- plyr::revalue(models$variables,
                                  c('Height.cm' = 'Height (cm)',
                                  'Weight.kg' = 'Weight (kg)',
                                  'Waist.cm' = 'Waist (cm)', 
                                  'Hip.cm' =  'Hip (cm)',
                                  'Waist.Hip' = 'Waist:Hip',
                                  'BMI' = 'Body Mass Index (kg/m2)',
                                  'Fasting.Glucose' = 'Fasting P. Glucose (mmol/L)',
                                  'Basal.Glucose' = 'Basal P. Glucose (mmol/L)',
                                  'Mean.Glucose' = 'Mean OGTT P. Glucose (mmol/L)',
                                  'glucose.2h' = 'OGTT P. Glucose 2-hours (mmol/L)',
                                  'Fasting.Insulin' = 'Fasting P. Insulin (pmol/L)',
                                  'Basal.Insulin' = 'Basal P. Insulin ((pmol/L)',
                                  'Mean.Insulin' = 'OGTT P. Insulin (pmol/L)',
                                  'Basal.InsulinSecretion' = 'Basal insulin secretion (pmol min-1m-2)',
                                  'Total.InsulinSecretion' = 'OGTT insulin secretion (nmol m-2)',
                                  'Glucose.Sensitivity' = 'Glucose sensitivity (pmol min-1 m-2 (mmol/L)-1)',
                                  'OGIS.2h'= 'OGTT insulin sensitivity (ml min-1m-2)',
                                  'Stumvoll' = "Stumvoll's index (umol min-1kg-1)",
                                  'Matsuda' = "Matsuda's index (umol min-1kg-1)",
                                  'Clinsb' = 'Basal insulin clearance (L min-1m-2)',
                                  'Clins' = 'OGTT insulin clearance (L min-1m-2)',
                                  'Active.GLP1.Conc.0m' = "Fasting P. GLP1 (pg/ml)",
                                  'Total.GLP1.Conc.0m' = "Total P. GLP1 (pg/ml)",
                                  'Total.GLP1.Conc.60m' = 'Total P. OGTT GLP1 (pg/ml)',
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
                                  'BP.S.Mean' = "Mean Systolic B.P. (mmHg)",
                                  'BP.D.Mean' = "Mean Diastolic B.P. (mmHg)",
                                  'Value.vm.hpf.mean.BL' = "Physical activity",
                                  'HDI' = "Healthy Dietary Index (HDI)",
                                  'lifestyle' = 'Lifestyle score',
                                  'bact.load'= "Bacterial cell count (cells/g feces)"))
  
pdf('vertical.pdf', height = 10, width = 8)
ggplot(models, aes(delta, variables, color = ifelse(fdr<=0.1, 'red', 'black'))) + 
  annotation_custom(grobTree(textGrob('High gene richness', x = .9, y = .5, gp=gpar(col = 'gray', fontsize = 25, fontface = 'bold.italic'), rot = 90))) +
  annotation_custom(grobTree(textGrob('Low gene richness', x = .1, y = .5, gp=gpar(col = 'gray', fontsize = 25, fontface = 'bold.italic'), rot = 90))) +
  geom_point() + 
  geom_segment(aes(x=delta, xend =0, y=variables, yend=variables)) +
  scale_y_discrete(limits = rev(levels(models$variables))) +
  scale_color_manual(values = c('black', 'red'), labels = c('Non-significant', 'Significant')) + xlim(-.5, .3) + geom_vline(xintercept = 0, lty = 'dashed') + 
  theme(axis.line = element_line(linetype = "solid"), panel.grid.major = element_line(colour = "gray78", linetype = "dashed"), 
        axis.title = element_text(size = 15), 
        axis.text.y = element_text(colour = "black"), 
        axis.text.x = element_text(colour = "black"), 
        plot.title = element_text(size = 17), 
        legend.title = element_blank(), panel.background = element_rect(fill = NA), legend.key = element_rect(fill = NA), 
        legend.background = element_rect(fill = NA), panel.border = element_rect(fill=NA)) +
  labs(title = "Lifestyle factors association with highest and lowest tertile of microbial gene richness", x = "Cliff's Delta effect size", y = NULL,caption = "Significant q-value (BH<0.1) is showed in red") 
dev.off()
