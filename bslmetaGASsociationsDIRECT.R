setwd('~/Desktop/01_Papers/04_DIRECT/stuff/baselinemetaGassoc')

associations.plot <- readRDS('associationsplot.rds')

associations.plot$variable <- factor(associations.plot$variable, levels=c('MGS', 'Species', 'Genera', 'GMM'))
associations.plot$association <- factor(associations.plot$association, 
                                        levels=c('waist', 'Fast.Insulin', 'Basal.InsSecr',
                                                'Total.InsSecr', 'Gluc.Sensitivity',
                                                'Matsuda', 'OGIS.2h', 'activeGLP1', 
                                                'totalGLP1.0', 'hscrp', 'hdl', 'ldl', 'tg', 'ast', 'chol', 
                                                'liver.iron', 'panc.iron', 'IAAT', 'ASAT', 'TAAT', 'liverfat', 'pancreasfat', 
                                                'physical', 'hdi', 'smoking'
                                                ))

associations.plot$association <- plyr::revalue(associations.plot$association,
                                  c('Height.cm' = 'Height (cm)',
                                    'Weight.kg' = 'Weight (kg)',
                                    'waist' = 'Waist (cm)', 
                                    'Hip.cm' =  'Hip (cm)',
                                    'Waist.Hip' = 'Waist:Hip',
                                    'BMI' = 'Body Mass Index (kg/m2)',
                                    'Fasting.Glucose' = 'Fasting P. Glucose (mmol/L)',
                                    'Basal.Glucose' = 'Basal P. Glucose (mmol/L)',
                                    'Mean.Glucose' = 'Mean OGTT P. Glucose (mmol/L)',
                                    'glucose.2h' = 'OGTT P. Glucose 2-hours (mmol/L)',
                                    'Fast.Insulin' = 'Fasting P. Insulin (pmol/L)',
                                    'Basal.Insulin' = 'Basal P. Insulin ((pmol/L)',
                                    'Mean.Insulin' = 'OGTT P. Insulin (pmol/L)',
                                    'Basal.InsSecr' = 'Basal insulin secretion (pmol min-1m-2)',
                                    'Total.InsSecr' = 'OGTT insulin secretion (nmol m-2)',
                                    'Gluc.Sensitivity' = 'Glucose sensitivity (pmol min-1 m-2 (mmol/L)-1)',
                                    'OGIS.2h'= 'OGTT insulin sensitivity (ml min-1m-2)',
                                    'Stumvoll' = "Stumvoll's index (umol min-1kg-1)",
                                    'Matsuda' = "Matsuda's index (umol min-1kg-1)",
                                    'Clinsb' = 'Basal insulin clearance (L min-1m-2)',
                                    'Clins' = 'OGTT insulin clearance (L min-1m-2)',
                                    'activeGLP1' = "Fasting P. GLP1 (pg/ml)",
                                    'totalGLP1.0' = "Total P. GLP1 (pg/ml)",
                                    'Total.GLP1.Conc.60m' = 'Total P. OGTT GLP1 (pg/ml)',
                                    'hscrp' = "Fasting P. hsCRP (mg/L)",
                                    'hdl' = "Fasting P. HDL  (mmol/L)",
                                    'ldl' = "Fasting P. LDL  (mmol/L)",
                                    'tg' = "Fasting P. Triglycerides  (mmol/L)",
                                    'Fasting.ALT' = "Fasting P. ALT (U/L)",
                                    'ast' = "Fasting P. AST (U/L)",
                                    'chol' = "Fasting P. Cholesterol  (mmol/L)",
                                    'liver.iron' = "Liver iron (umol/g)",
                                    'panc.iron' = "Pancreatic iron (umol/g)",
                                    'IAAT' = "IAAT (L)",
                                    'ASAT' = "ASAT (L)",
                                    'TAAT' = "TAAT (L)",
                                    'liverfat' = "Liver fat (%)",
                                    'pancreasfat' = "Pancreatic fat (%)",
                                    'BP.S.Mean' = "Mean Systolic B.P. (mmHg)",
                                    'BP.D.Mean' = "Mean Diastolic B.P. (mmHg)",
                                    'physical' = "Physical activity",
                                    'hdi' = "Healthy Dietary Index (HDI)",
                                    'lifestyle' = 'Lifestyle score',
                                    'smoking'= "Smoking"))
library(ggplot2)
pdf('phenoassoc.pdf', width = 20, height = 9)
ggplot(associations.plot, aes(n.dir, association, size=n, shape = direction, fill = n.dir)) + geom_point() +
  scale_shape_manual(values = c(24, 25)) +
  scale_fill_gradient2(low='blue', high = 'red') +
  scale_x_continuous(breaks = seq(-30,30, by=5)) +
  scale_y_discrete(limits = rev(levels(associations.plot$association))) +
  geom_vline(xintercept = 0, linetype='dashed') +
  facet_wrap(~variable, scales = 'free', ncol=4) +
  theme(panel.border = element_rect(fill=NA), panel.background = element_blank(),
        panel.grid.major = element_line(color='gray85', linetype='dashed'),
        axis.title = element_text(size = 13), 
        legend.text = element_text(face = "bold"), legend.title = element_text(face = "bold"),
        legend.background = element_blank(), legend.box.background = element_blank(), legend.key = element_blank(),
        strip.background = element_blank(), strip.text = element_text(face='bold', size=18)) +
  labs(x = 'Number of associations', y = NULL, fill = "Number of \nassociations",
       shape = "Direction") +
  guides(size=F)
dev.off()
