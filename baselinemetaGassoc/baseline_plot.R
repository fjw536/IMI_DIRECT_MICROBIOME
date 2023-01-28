
associations.plot <- data.frame(rbind(mgs, spps, genus, gmm), 'variable'=c(rep('MGS', 40), rep('Species', 32), 
                                                      rep('Genera', 24), rep('GMM', 24)))
associations.plot$variable <- factor(associations.plot$variable, levels=c('MGS', 'Species', 'Genera', 'GMM'))

setwd('~/Desktop/01_Papers/04_DIRECT/stuff/baselinemetaGassoc')
associations.plot <- readRDS('associationsplot.rds')
associations.plot$variable <- factor(associations.plot$variable, levels=c('MGS', 'Species', 'Genera', 'GMM'))
associations.plot <- associations.plot[associations.plot$variable!='Species',]

associations.plot$association <- plyr::revalue(associations.plot$association,
              c('waist' = 'Waist (cm)', 
                'hip' =  'Hip (cm)',
                'Fast.Insulin' = 'Fasting P. Insulin (pmol/L)',
                'Basal.InsSecr' = 'Basal insulin secretion (pmol min-1m-2)',
                'Total.InsSecr' = 'OGTT insulin secretion (nmol m-2)',
                'Gluc.Sensitivity' = 'Glucose sensitivity (pmol min-1 m-2 (mmol/L)-1)',
                'OGIS.2h'= 'OGTT insulin sensitivity (ml min-1m-2)',
                'Matsuda' = "Matsuda's index (umol min-1kg-1)",
                'activeGLP1' = "Fasting P. GLP1 (pg/ml)",
                'totalGLP1.0' = 'Total P. OGTT GLP1 (pg/ml)',
                'hscrp' = "Fasting P. hsCRP (mg/L)",
                'hdl' = "Fasting P. HDL  (mmol/L)",
                'ldl' = "Fasting P. LDL  (mmol/L)",
                'tg' = "Fasting P. Triglycerides  (mmol/L)",
                'ast' = "Fasting P. AST (U/L)",
                'chol' = "Fasting P. Cholesterol  (mmol/L)",
                'liver.iron' = "Liver iron (umol/g)",
                'panc.iron' = "Pancreatic iron (umol/g)",
                'IAAT' = "IAAT (L)",
                'ASAT' = "ASAT (L)",
                'TAAT' = "TAAT (L)",
                'liverfat' = "Liver fat (%)",
                'pancreasfat' = "Pancreatic fat (%)",
                'hdi' = "Healthy Dietary Index (HDI)", 
                'smoking' = 'Smoking'))

associations.plot$association <- factor(associations.plot$association, 
                           levels = c('Waist (cm)', 
                                      'Hip (cm)',
                                      'Fasting P. Insulin (pmol/L)',
                                      'Basal insulin secretion (pmol min-1m-2)',
                                      'OGTT insulin secretion (nmol m-2)',
                                      'Glucose sensitivity (pmol min-1 m-2 (mmol/L)-1)',
                                      'OGTT insulin sensitivity (ml min-1m-2)',
                                      "Matsuda's index (umol min-1kg-1)",
                                      "Fasting P. GLP1 (pg/ml)",
                                      'Total P. OGTT GLP1 (pg/ml)',
                                      "Fasting P. hsCRP (mg/L)",
                                      "Fasting P. HDL  (mmol/L)",
                                      "Fasting P. LDL  (mmol/L)",
                                      "Fasting P. Triglycerides  (mmol/L)",
                                      "Fasting P. AST (U/L)",
                                      "Fasting P. Cholesterol  (mmol/L)",
                                      "Liver iron (umol/g)",
                                      "Pancreatic iron (umol/g)",
                                      "IAAT (L)",
                                      "ASAT (L)",
                                      "TAAT (L)",
                                      "Liver fat (%)",
                                      "Pancreatic fat (%)",
                                      "Healthy Dietary Index (HDI)", 
                                      'Smoking'))
associations.plot$variable <- plyr::revalue(associations.plot$variable, 
                                            c('MGS' = 'Metagenomic species', 'GMM' = 'Gut Metabolic Modules'))
ggplot(associations.plot, aes(n.dir, association, size=n, shape = direction, fill = n.dir)) + geom_point() +
  scale_shape_manual(values = c(24, 25)) +
  scale_fill_gradient2(low='blue', high = 'red') +
  scale_x_continuous(breaks = seq(-30,30, by=5)) +
  scale_y_discrete(limits = rev(levels(associations.plot$association))) +
  geom_vline(xintercept = 0, linetype='dashed') +
  facet_wrap(~variable, scales = 'free', ncol=4) +
  theme(panel.border = element_rect(fill=NA), panel.background = element_blank(),
      panel.grid.major = element_line(color='gray85', linetype='dashed'),
      axis.text.x = element_text(color='black'),
      legend.background = element_blank(), legend.box.background = element_blank(), legend.key = element_blank(),
      strip.background = element_blank(), strip.text = element_text(size=15)) +
  labs(x = 'Number of associations', y = NULL, fill = "Number of\nassociations",
     shape = "Direction") +
  guides(size=F)
