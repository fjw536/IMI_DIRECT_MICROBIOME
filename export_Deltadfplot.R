

delta.df <- readRDS('Downloads/data (24)/deltadfplot.rds')

delta.df <- delta.df %>%
  select(-c('Age', "CenterID", 'Gly.Cat.2')) %>% 
  group_by(GlyCat3g) %>% 
  dplyr::summarise(across(is.numeric, list(mean=mean, sd=sd, median=median))) %>% 
  pivot_longer(ends_with(c("_mean", "_sd", "_median")), names_to = c("name", ".value"), names_sep = "_") %>% 
  mutate(GlyCat3g = factor(GlyCat3g, levels=c('NGR', 'IGR', 'T2D')))%>% 
  as.data.frame()
delta.df$name <- factor(delta.df$name, levels=unique(delta.df$name))
delta.df$name <- plyr::revalue(delta.df$name ,
                                  c('Height' = 'Height (cm)',
                                    'Weight' = 'Weight (kg)',
                                    'Waist' = 'Waist (cm)', 
                                    'Hip.cm' =  'Hip (cm)',
                                    'Waist.Hip' = 'Waist:Hip',
                                    'BMI' = 'Body Mass Index (kg/m2)',
                                    'Fast.Glucose' = 'Fasting P. Glucose (mmol/L)',
                                    'Basal.Glucose' = 'Basal P. Glucose (mmol/L)',
                                    'Mean.Glucose' = 'Mean OGTT P. Glucose (mmol/L)',
                                    'Glucose.2h' = 'OGTT P. Glucose 2-hours (mmol/L)',
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
                                    'bact.load'= "Bacterial cell count (cells/g feces)"))

ggplot(delta.df, aes(reorder(name, desc(name)), mean, fill=GlyCat3g, color=GlyCat3g)) + geom_col(position = position_dodge2(reverse = TRUE), alpha = .6) +
  scale_color_manual(values = c('darkolivegreen4', 'red', 'purple')) +
  scale_fill_manual(values = c('darkolivegreen4', 'red', 'purple')) +
  coord_flip() + ylab(label = 'Mean delta value') + xlab(NULL) +
  theme(panel.grid.major = element_line(colour = "gray85",linetype = "dashed"), 
        axis.text = element_text(size = 11, colour = "black"), 
        axis.title = element_text(size=12,  color='black'),
        panel.background = element_rect(fill = NA), legend.title = element_text(),
        legend.text = element_text(), strip.background = element_rect(fill=NA),
        strip.text = element_text(face='bold', size =12), panel.border = element_rect(fill=NA)) +
  labs(y='Difference of means between\nbaseline and follow-up', fill = "Glycaemic\nstatus", color = "Glycaemic\nstatus")
