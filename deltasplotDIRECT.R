deltas <- read.table('~/Downloads/export2_MetClusters_deltavalues.txt', header=T, sep='\t')

heat.sort <- deltas[deltas$delta == 'Mean.Insulin',]
heat.sort <- heat.sort[order(heat.sort$value),]
heat.sort <- factor(heat.sort$module.name, levels=heat.sort$module.name)

order.bch <- c('Height', 'Waist',
               'Fast.Glucose', 'Basal.Glucose', 'Mean.Glucose', 'Glucose.2h', 'Fast.Insulin', 'Basal.Insulin',
               'Mean.Insulin', 'Basal.InsSecr', 'Total.InsSecr', 'Glucose.Sensitivity', 'OGIS.2h',
               'Stumvoll', 'Matsuda', 'BP.Sys', 'BP.Dias') 

deltas$delta <- plyr::revalue(factor(deltas$delta, levels = order.bch),
                          c('Height' = 'Height (cm)',
                            'Weight' = 'Weight (kg)',
                            'Waist' = 'Waist (cm)', 
                            'Fast.Glucose' = 'Fasting P. Glucose (mmol/L)',
                            'Basal.Glucose' = 'Basal P. Glucose (mmol/L)',
                            'Mean.Glucose' = 'Mean OGTT P. Glucose (mmol/L)',
                            'Glucose.2h' = 'OGTT P. Glucose 2-hours (mmol/L)',
                            'Fast.Insulin' = 'Fasting P. Insulin (pmol/L)',
                            'Basal.Insulin' = 'Basal P. Insulin ((pmol/L)',
                            'Mean.Insulin' = 'OGTT P. Insulin (pmol/L)',
                            'Basal.InsSecr' = 'Basal insulin secretion (pmol min-1m-2)',
                            'Total.InsSecr' = 'OGTT insulin secretion (nmol m-2)',
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
                            'BP.Sys' = "Mean Systolic B.P. (mmHg)",
                            'BP.Dias' = "Mean Diastolic B.P. (mmHg)",
                            'Value.vm.hpf.mean.BL' = "Physical activity",
                            'HDI' = "Healthy Dietary Index (HDI)",
                            'lifestyle' = 'Lifestyle score',
                            'bact.load'= "Bacterial cell count (cells/g feces)",
                            'GR' = 'Gene richness'))

pdf('04_DIRECT/Manuscript_v25052021/supplfigures/SF20.pdf', height = 12, width = 8)
deltas %>% 
  filter(fdr<=.1) %>% 
  mutate('module.name' = factor(module.name, levels=heat.sort)) %>%
  ggplot(aes(delta, module.name, fill=value)) + 
  geom_tile(color = 'black', size =0.1) + 
  scale_fill_gradient2(low='blue', high='red') +
  scale_x_discrete(limits=rev(levels(deltas$delta))) +
  #coord_flip() + #guides(fill=F) +
  theme(panel.background = element_blank(), panel.border = element_rect(fill=NA), 
        panel.grid.major = element_line(color='gray90'), 
        axis.text.x = element_text(angle = 90, vjust = .9, hjust = 1)) +#,
  # legend.position = 'bottom') +
  labs(y='Metabolite clusters', x = 'Phenotypical variables progression (delta values)',
       fill='rho')
dev.off()
