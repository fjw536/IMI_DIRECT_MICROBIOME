############################################################################################################
######################## Association of baseline metaG to baseline biochemical  ############################
############################################################################################################

# directory
setwd('/home/scratch/margar')

# libraries
library(effsize)
library(phyloseq)
library(caret); library(dplyr)

# data
qmp.bsl <-
  read.table(
    'data/WP2.1_775_IGR_MGS-QMP_igrHA_Ent_M.csv',
    header = T,
    row.names = 1
  )
qmp.bsl <- qmp.bsl[, 7:727]

mtdt.bsl <- read.table('data/curatedmetadata/WP2.1_BL_MetadataWithInfl.txt', header=T, sep='\t')

qmp.bsl <- qmp.bsl[row.names(qmp.bsl) %in% mtdt.bsl$samplerenamed, ]
mtdt.bsl$samplerenamed == row.names(qmp.bsl)
names(qmp.bsl) <- gsub('\\.', ':', names(qmp.bsl))
# names(qmp.bsl) <- gsub(':', '\\.', names(qmp.bsl))

tax <-
  read.table(
    'data/taxonomy.tsv',
    header = T,
    sep = '\t',
    row.names = 1
  )
tax.m <-
  as.matrix(tax)
rownames(tax.m) <- row.names(tax)
colnames(tax.m) <- names(tax)

physeq <-
  phyloseq(otu_table(qmp.bsl, taxa_are_rows = F), tax_table(tax.m))

mtdt.bsl$Gly.Cat2 <- NULL
mtdt.bsl$ADA <- NULL
mtdt.bsl$enterotype <- NULL
mtdt.bsl$DSGARLTertileGroup <- NULL
mtdt.bsl$DSGeneAbundRichness <- NULL
mtdt.bsl$DSMGSAbundShannon <- NULL
mtdt.bsl$Gly.Cat <- NULL
mtdt.bsl$Sex <- NULL
mtdt.bsl <- mtdt.bsl[,4:49]
mtdt.bsl$Smoking.Status.BL <- plyr::revalue(as.factor(mtdt.bsl$Smoking.Status.BL), c('ex-smoker'='non-smoker', 'never'='non-smoker', 'current-smoker'='smoker'))

preProcValues.bsl <- preProcess(mtdt.bsl %>% 
                                  select(Age, Height.cm, Weight.kg, Waist.cm, Hip.cm, Waist.Hip, BMI,
                                         Fasting.Glucose, Basal.Glucose, Mean.Glucose, glucose.2h,
                                         Fasting.Insulin, Basal.Insulin, Mean.Insulin, Basal.InsulinSecretion,
                                         Total.InsulinSecretion, Glucose.Sensitivity, OGIS.2h, Stumvoll, 
                                         Matsuda, Clinsb, Clins, Active.GLP1.Conc.0m, Total.GLP1.Conc.0m, Total.GLP1.Conc.60m,
                                         hsCRP, Fasting.HDL, Fasting.LDL, Fasting.TG, Fasting.ALT, Fasting.AST, Fasting.Chol,
                                         Liver.Iron, Panc.Iron, IAAT, ASAT, TAAT, LiverFat.., PancreasFat.., BP.S.Mean, BP.D.Mean,
                                         Value.vm.hpf.mean.BL, HDI),
                                method = c('knnImpute'),
                                k = 25,
                                knnSummary = mean)
mtdt.bsl <- predict(preProcValues.bsl, mtdt.bsl, na.action=na.pass)

procNames.bsl <- data.frame(col = names(preProcValues.bsl$mean), mean = preProcValues.bsl$mean, sd = preProcValues.bsl$std)
for(i in procNames.bsl$col){
  mtdt.bsl[i] <- mtdt.bsl[i]*preProcValues.bsl$std[i]+preProcValues.bsl$mean[i] 
}

# MGS modelling -----
progression.model <- data.frame(mtdt.bsl, qmp.bsl)
# 
# #
physeq.spps <- tax_glom(physeq, taxrank = rank_names(physeq)[7])
spps <- as.data.frame(physeq.spps@otu_table)

progression.model <- data.frame(mtdt.bsl, spps)
# 
# #
physeq.gen <- tax_glom(physeq, taxrank = rank_names(physeq)[6])
genus <- as.data.frame(physeq.gen@otu_table)

progression.model <- data.frame(mtdt.bsl, genus)
#
GMMs <- read.table('data/WP2.1_775_IGR_GMM-QMP_igrHA_Ent_M.csv', header=T, sep='\t', row.names=1)
GMMs <- GMMs[,7:103]
mtdt.bsl.original <- read.table('data/curatedmetadata/WP2.1_BL_MetadataWithInfl.txt', header=T, sep='\t')
GMMs <- GMMs[row.names(GMMs)%in%mtdt.bsl.original$samplerenamed,]

progression.model <- data.frame(mtdt.bsl, GMMs)

# MGS vs smoking
result.smoking <- list()
for (x in seq(from = 48, to = ncol(progression.model))) {
  print(names(progression.model[x]))
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:47])
  mod.a <-
    lm(mgs ~ Smoking.Status.BL + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Smoking.Status.BLnon-smoker', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Smoking.Status.BLnon-smoker', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Smoking.Status.BL)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.smoking[[x]] <- result
}
smoking.lm <- do.call(rbind.data.frame, result.smoking)
names(smoking.lm) <-
  c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
smoking.lm$fdr <- p.adjust(smoking.lm$pval, method = 'fdr')


# MGS vs Height
result.height <- list()
for (x in seq(from = 48, to = ncol(progression.model))) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:47])
  mod.a <-
    lm(mgs ~ Height.cm + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Height.cm', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Height.cm', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Height.cm)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.height[[x]] <- result
}
height.lm <- do.call(rbind.data.frame, result.height)
names(height.lm) <-
  c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
height.lm$fdr <- p.adjust(height.lm$pval, method = 'fdr')

# MGS vs Weight
result.weight <- list()
for (x in seq(from = 48, to = ncol(progression.model))) {
  print(names(progression.model[x]))
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:47])
  mod.a <-
    lm(mgs ~ Weight.kg + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Weight.kg', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Weight.kg', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Weight.kg)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.weight[[x]] <- result
}
weight.lm <- do.call(rbind.data.frame, result.weight)
names(weight.lm) <-
  c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
weight.lm$fdr <- p.adjust(weight.lm$pval, method = 'fdr')

# MGS vs Waist
result.waist <- list()
for (x in seq(from = 48, to = ncol(progression.model))) {
  print(names(progression.model[x]))
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:47])
  mod.a <-
    lm(mgs ~ Waist.cm + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Waist.cm', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Waist.cm', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Waist.cm)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.waist[[x]] <- result
}
waist.lm <- do.call(rbind.data.frame, result.waist)
names(waist.lm) <-
  c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
waist.lm$fdr <- p.adjust(waist.lm$pval, method = 'fdr')

# MGS vs hip
result.hip <- list()
for (x in seq(from = 48, to = ncol(progression.model))) {
  print(names(progression.model[x]))
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:47])
  mod.a <-
    lm(mgs ~ Hip.cm + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Hip.cm', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Hip.cm', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Hip.cm)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.hip[[x]] <- result
}
hip.lm <- do.call(rbind.data.frame, result.hip)
names(hip.lm) <-
  c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
hip.lm$fdr <- p.adjust(hip.lm$pval, method = 'fdr')

# MGS vs Waist:Hip
result.waist.hip <- list()
for (x in seq(from = 48, to = ncol(progression.model))) {
  print(names(progression.model[x]))
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:47])
  mod.a <-
    lm(mgs ~ Waist.Hip + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Waist.Hip', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Waist.Hip', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Waist.Hip)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.waist.hip[[x]] <- result
}
waist.hip.lm <- do.call(rbind.data.frame, result.waist.hip)
names(waist.hip.lm) <-
  c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
waist.hip.lm$fdr <- p.adjust(waist.hip.lm$pval, method = 'fdr')

# MGS vs BMI
result.bmi <- list()
for (x in seq(from = 48, to = ncol(progression.model))) {
  print(names(progression.model[x]))
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:47])
  mod.a <-
    lm(mgs ~ BMI + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['BMI', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['BMI', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$BMI)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.bmi[[x]] <- result
}
bmi.lm <- do.call(rbind.data.frame, result.bmi)
names(bmi.lm) <- c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
bmi.lm$fdr <- p.adjust(bmi.lm$pval, method = 'fdr')

# MGS vs Fast.Glucose
result.Fast.Glucose <- list()
for (x in seq(from = 48, to = ncol(progression.model))) {
  print(names(progression.model[x]))
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:47])
  mod.a <-
    lm(mgs ~ Fasting.Glucose + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Fasting.Glucose', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Fasting.Glucose', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Fasting.Glucose)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.Fast.Glucose[[x]] <- result
}
Fast.Glucose.lm <- do.call(rbind.data.frame, result.Fast.Glucose)
names(Fast.Glucose.lm) <-
  c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
Fast.Glucose.lm$fdr <- p.adjust(Fast.Glucose.lm$pval, method = 'fdr')
Fast.Glucose.lm[Fast.Glucose.lm$fdr < .1, ]

# MGS vs Basal.Glucose
result.Basal.Glucose <- list()
for (x in seq(from = 48, to = ncol(progression.model))) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:47])
  mod.a <-
    lm(mgs ~ Basal.Glucose + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Basal.Glucose', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Basal.Glucose', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Basal.Glucose)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.Basal.Glucose[[x]] <- result
}
Basal.Glucose.lm <- do.call(rbind.data.frame, result.Basal.Glucose)
names(Basal.Glucose.lm) <-
  c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
Basal.Glucose.lm$fdr <-
  p.adjust(Basal.Glucose.lm$pval, method = 'fdr')
Basal.Glucose.lm[Basal.Glucose.lm$fdr < .1, ]

# MGS vs Mean.Glucose
result.Mean.Glucose <- list()
for (x in seq(from = 48, to = ncol(progression.model))) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:47])
  mod.a <-
    lm(mgs ~ Mean.Glucose + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Mean.Glucose', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Mean.Glucose', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Mean.Glucose)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.Mean.Glucose[[x]] <- result
}
Mean.Glucose.lm <- do.call(rbind.data.frame, result.Mean.Glucose)
names(Mean.Glucose.lm) <-
  c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
Mean.Glucose.lm$fdr <- p.adjust(Mean.Glucose.lm$pval, method = 'fdr')
Mean.Glucose.lm[Mean.Glucose.lm$fdr < .1, ]

# MGS vs Glucose.2h
result.Glucose.2h <- list()
for (x in seq(from = 48, to = ncol(progression.model))) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:47])
  mod.a <-
    lm(mgs ~ glucose.2h + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['glucose.2h', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['glucose.2h', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$glucose.2h)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.Glucose.2h[[x]] <- result
}
Glucose.2h.lm <- do.call(rbind.data.frame, result.Glucose.2h)
names(Glucose.2h.lm) <-
  c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
Glucose.2h.lm$fdr <- p.adjust(Glucose.2h.lm$pval, method = 'fdr')

# MGS vs Fast.Insulin
result.Fast.Insulin <- list()
for (x in seq(from = 48, to = ncol(progression.model))) {
  print(names(progression.model[x]))
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:47])
  mod.a <-
    lm(mgs ~ Fasting.Insulin + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Fasting.Insulin', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Fasting.Insulin', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Fasting.Insulin)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.Fast.Insulin[[x]] <- result
}
Fast.Insulin.lm <- do.call(rbind.data.frame, result.Fast.Insulin)
names(Fast.Insulin.lm) <-
  c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
Fast.Insulin.lm$fdr <- p.adjust(Fast.Insulin.lm$pval, method = 'fdr')
Fast.Insulin.lm[Fast.Insulin.lm$fdr < .1, ]

# MGS vs Basal.Insulin
result.Basal.Insulin <- list()
for (x in seq(from = 48, to = ncol(progression.model))) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:47])
  mod.a <-
    lm(mgs ~ Basal.Insulin + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Basal.Insulin', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Basal.Insulin', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Basal.Insulin)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.Basal.Insulin[[x]] <- result
}
Basal.Insulin.lm <- do.call(rbind.data.frame, result.Basal.Insulin)
names(Basal.Insulin.lm) <-
  c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
Basal.Insulin.lm$fdr <-
  p.adjust(Basal.Insulin.lm$pval, method = 'fdr')
Basal.Insulin.lm[Basal.Insulin.lm$fdr < .1, ]

# MGS vs Mean.Insulin
result.Mean.Insulin <- list()
for (x in seq(from = 48, to = ncol(progression.model))) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:47])
  mod.a <-
    lm(mgs ~ Mean.Insulin + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Mean.Insulin', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Mean.Insulin', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Mean.Insulin)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  
  result.Mean.Insulin[[x]] <- result
}
Mean.Insulin.lm <- do.call(rbind.data.frame, result.Mean.Insulin)
names(Mean.Insulin.lm) <-
  c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
Mean.Insulin.lm$fdr <- p.adjust(Mean.Insulin.lm$pval, method = 'fdr')
Mean.Insulin.lm[Mean.Insulin.lm$fdr < .1, ]

# MGS vs Basal.InsSecr
result.Basal.InsSecr <- list()
for (x in seq(from = 48, to = ncol(progression.model))) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:47])
  mod.a <-
    lm(mgs ~ Basal.InsulinSecretion + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Basal.InsulinSecretion', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Basal.InsulinSecretion', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Basal.InsulinSecretion)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.Basal.InsSecr[[x]] <- result
}
Basal.InsSecr.lm <- do.call(rbind.data.frame, result.Basal.InsSecr)
names(Basal.InsSecr.lm) <-
  c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
Basal.InsSecr.lm$fdr <-
  p.adjust(Basal.InsSecr.lm$pval, method = 'fdr')
Basal.InsSecr.lm[Basal.InsSecr.lm$fdr < .1, ]

# MGS vs Total.InsSecr
result.Total.InsSecr <- list()
for (x in seq(from = 48, to = ncol(progression.model))) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:47])
  mod.a <-
    lm(mgs ~ Total.InsulinSecretion + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Total.InsulinSecretion', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Total.InsulinSecretion', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Total.InsulinSecretion)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.Total.InsSecr[[x]] <- result
}
Total.InsSecr.lm <- do.call(rbind.data.frame, result.Total.InsSecr)
names(Total.InsSecr.lm) <-
  c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
Total.InsSecr.lm$fdr <-
  p.adjust(Total.InsSecr.lm$pval, method = 'fdr')
Total.InsSecr.lm[Total.InsSecr.lm$fdr < .1, ]

# MGS vs Gluc.Sensitivity
result.Gluc.Sensitivity <- list()
for (x in seq(from = 48, to = ncol(progression.model))) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:47])
  mod.a <-
    lm(
      mgs ~ Glucose.Sensitivity + Gender + CenterID + Age,
      data = model,
      na.action = na.omit
    )
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Glucose.Sensitivity', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Glucose.Sensitivity', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Glucose.Sensitivity)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.Gluc.Sensitivity[[x]] <- result
}
Gluc.Sensitivity.lm <-
  do.call(rbind.data.frame, result.Gluc.Sensitivity)
names(Gluc.Sensitivity.lm) <-
  c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
Gluc.Sensitivity.lm$fdr <-
  p.adjust(Gluc.Sensitivity.lm$pval, method = 'fdr')
Gluc.Sensitivity.lm[Gluc.Sensitivity.lm$fdr < .1, ]

# MGS vs OGIS.2h
result.OGIS.2h <- list()
for (x in seq(from = 48, to = ncol(progression.model))) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:47])
  mod.a <-
    lm(mgs ~ OGIS.2h + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['OGIS.2h', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['OGIS.2h', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$OGIS.2h)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  
  result.OGIS.2h[[x]] <- result
}
OGIS.2h.lm <- do.call(rbind.data.frame, result.OGIS.2h)
names(OGIS.2h.lm) <-
  c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
OGIS.2h.lm$fdr <- p.adjust(OGIS.2h.lm$pval, method = 'fdr')
OGIS.2h.lm[OGIS.2h.lm$fdr < .1, ]

# MGS vs Stumvoll
result.Stumvoll <- list()
for (x in seq(from = 48, to = ncol(progression.model))) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:47])
  mod.a <-
    lm(mgs ~ Stumvoll + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Stumvoll', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Stumvoll', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Stumvoll)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  
  result.Stumvoll[[x]] <- result
}
Stumvoll.lm <- do.call(rbind.data.frame, result.Stumvoll)
names(Stumvoll.lm) <-
  c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
Stumvoll.lm$fdr <- p.adjust(Stumvoll.lm$pval, method = 'fdr')

# MGS vs Matsuda
result.Matsuda <- list()
for (x in seq(from = 48, to = ncol(progression.model))) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:47])
  mod.a <-
    lm(mgs ~ Matsuda + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Matsuda', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Matsuda', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Matsuda)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.Matsuda[[x]] <- result
}
Matsuda.lm <- do.call(rbind.data.frame, result.Matsuda)
names(Matsuda.lm) <-
  c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
Matsuda.lm$fdr <- p.adjust(Matsuda.lm$pval, method = 'fdr')
Matsuda.lm[Matsuda.lm$fdr < .1, ]

# MGS vs Clinsb
result.Clinsb <- list()
for (x in seq(from = 48, to = ncol(progression.model))) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:47])
  mod.a <-
    lm(mgs ~ Clinsb + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Clinsb', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Clinsb', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Clinsb)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.Clinsb[[x]] <- result
}
Clinsb.lm <- do.call(rbind.data.frame, result.Clinsb)
names(Clinsb.lm) <-
  c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
Clinsb.lm$fdr <- p.adjust(Clinsb.lm$pval, method = 'fdr')

# MGS vs Clins
result.Clins <- list()
for (x in seq(from = 48, to = ncol(progression.model))) {
  print(names(progression.model[x]))
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:47])
  mod.a <-
    lm(mgs ~ Clins + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Clins', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Clins', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Clins)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.Clins[[x]] <- result
}
Clins.lm <- do.call(rbind.data.frame, result.Clins)
names(Clins.lm) <-
  c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
Clins.lm$fdr <- p.adjust(Clins.lm$pval, method = 'fdr')

# MGS vs activeGLP1
result.activeGLP1 <- list()
for (x in seq(from = 48, to = ncol(progression.model))) {
  print(names(progression.model[x]))
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:47])
  mod.a <-
    lm(mgs ~ Active.GLP1.Conc.0m + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Active.GLP1.Conc.0m', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Active.GLP1.Conc.0m', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Active.GLP1.Conc.0m)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.activeGLP1[[x]] <- result
}
activeGLP1.lm <- do.call(rbind.data.frame, result.activeGLP1)
names(activeGLP1.lm) <-
  c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
activeGLP1.lm$fdr <- p.adjust(activeGLP1.lm$pval, method = 'fdr')

# MGS vs totalGLP1.0
result.totalGLP1.0 <- list()
for (x in seq(from = 48, to = ncol(progression.model))) {
  print(names(progression.model[x]))
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:47])
  mod.a <-
    lm(mgs ~ Total.GLP1.Conc.0m + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Total.GLP1.Conc.0m', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Total.GLP1.Conc.0m', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Total.GLP1.Conc.0m)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.totalGLP1.0[[x]] <- result
}
totalGLP1.0.lm <- do.call(rbind.data.frame, result.totalGLP1.0)
names(totalGLP1.0.lm) <-
  c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
totalGLP1.0.lm$fdr <- p.adjust(totalGLP1.0.lm$pval, method = 'fdr')

# MGS vs totalGLP1.0
result.totalGLP1.6 <- list()
for (x in seq(from = 48, to = ncol(progression.model))) {
  print(names(progression.model[x]))
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:47])
  mod.a <-
    lm(mgs ~ Total.GLP1.Conc.60m + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Total.GLP1.Conc.60m', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Total.GLP1.Conc.60m', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Total.GLP1.Conc.60m)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.totalGLP1.60[[x]] <- result
}
totalGLP1.60.lm <- do.call(rbind.data.frame, result.totalGLP1.60)
names(totalGLP1.60.lm) <-
  c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
totalGLP1.60.lm$fdr <- p.adjust(totalGLP1.60.lm$pval, method = 'fdr')

# MGS vs hsCRP
result.hscrp <- list()
for (x in seq(from = 48, to = ncol(progression.model))) {
  print(names(progression.model[x]))
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:47])
  mod.a <-
    lm(mgs ~ hsCRP + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['hsCRP', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['hsCRP', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$hsCRP)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.hscrp[[x]] <- result
}
hscrp.lm <- do.call(rbind.data.frame, result.hscrp)
names(hscrp.lm) <-
  c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
hscrp.lm$fdr <- p.adjust(hscrp.lm$pval, method = 'fdr')

# MGS vs Fasting.HDL
result.hdl <- list()
for (x in seq(from = 48, to = ncol(progression.model))) {
  print(names(progression.model[x]))
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:47])
  mod.a <-
    lm(mgs ~ Fasting.HDL + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Fasting.HDL', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Fasting.HDL', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Fasting.HDL)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.hdl[[x]] <- result
}
hdl.lm <- do.call(rbind.data.frame, result.hdl)
names(hdl.lm) <-
  c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
hdl.lm$fdr <- p.adjust(hdl.lm$pval, method = 'fdr')

# MGS vs Fasting.LDL
result.ldl <- list()
for (x in seq(from = 48, to = ncol(progression.model))) {
  print(names(progression.model[x]))
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:47])
  mod.a <-
    lm(mgs ~ Fasting.LDL + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Fasting.LDL', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Fasting.LDL', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Fasting.LDL)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.ldl[[x]] <- result
}
ldl.lm <- do.call(rbind.data.frame, result.ldl)
names(ldl.lm) <-
  c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
ldl.lm$fdr <- p.adjust(ldl.lm$pval, method = 'fdr')

# MGS vs Fasting.TG
result.tg <- list()
for (x in seq(from = 48, to = ncol(progression.model))) {
  print(names(progression.model[x]))
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:47])
  mod.a <-
    lm(mgs ~ Fasting.TG + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Fasting.TG', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Fasting.TG', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Fasting.TG)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.tg[[x]] <- result
}
tg.lm <- do.call(rbind.data.frame, result.tg)
names(tg.lm) <-
  c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
tg.lm$fdr <- p.adjust(tg.lm$pval, method = 'fdr')

# MGS vs Fasting.AST
result.ast <- list()
for (x in seq(from = 48, to = ncol(progression.model))) {
  print(names(progression.model[x]))
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:47])
  mod.a <-
    lm(mgs ~ Fasting.AST + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Fasting.AST', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Fasting.AST', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Fasting.AST)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.ast[[x]] <- result
}
ast.lm <- do.call(rbind.data.frame, result.ast)
names(ast.lm) <-
  c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
ast.lm$fdr <- p.adjust(ast.lm$pval, method = 'fdr')

# MGS vs Fasting.ALT
result.alt <- list()
for (x in seq(from = 48, to = ncol(progression.model))) {
  print(names(progression.model[x]))
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:47])
  mod.a <-
    lm(mgs ~ Fasting.ALT + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Fasting.ALT', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Fasting.ALT', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Fasting.ALT)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.alt[[x]] <- result
}
alt.lm <- do.call(rbind.data.frame, result.alt)
names(alt.lm) <-
  c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
asl.lm$fdr <- p.adjust(alt.lm$pval, method = 'fdr')

# MGS vs Fasting.Chol
result.chol <- list()
for (x in seq(from = 48, to = ncol(progression.model))) {
  print(names(progression.model[x]))
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:47])
  mod.a <-
    lm(mgs ~ Fasting.Chol + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Fasting.Chol', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Fasting.Chol', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Fasting.Chol)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.chol[[x]] <- result
}
chol.lm <- do.call(rbind.data.frame, result.chol)
names(chol.lm) <-
  c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
chol.lm$fdr <- p.adjust(chol.lm$pval, method = 'fdr')

# MGS vs Liver.Iron
result.liver.iron <- list()
for (x in seq(from = 48, to = ncol(progression.model))) {
  print(names(progression.model[x]))
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:47])
  mod.a <-
    lm(mgs ~ Liver.Iron + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Liver.Iron', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Liver.Iron', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Liver.Iron)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.liver.iron[[x]] <- result
}
liver.iron.lm <- do.call(rbind.data.frame, result.liver.iron)
names(liver.iron.lm) <-
  c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
liver.iron.lm$fdr <- p.adjust(liver.iron.lm$pval, method = 'fdr')

# MGS vs Panc.Iron
result.panc.iron <- list()
for (x in seq(from = 48, to = ncol(progression.model))) {
  print(names(progression.model[x]))
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:47])
  mod.a <-
    lm(mgs ~ Panc.Iron + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Panc.Iron', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Panc.Iron', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Panc.Iron)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.panc.iron[[x]] <- result
}
panc.iron.lm <- do.call(rbind.data.frame, result.panc.iron)
names(panc.iron.lm) <-
  c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
panc.iron.lm$fdr <- p.adjust(panc.iron.lm$pval, method = 'fdr')

# MGS vs IAAT
result.IAAT <- list()
for (x in seq(from = 48, to = ncol(progression.model))) {
  print(names(progression.model[x]))
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:47])
  mod.a <-
    lm(mgs ~ IAAT + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['IAAT', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['IAAT', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$IAAT)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.IAAT[[x]] <- result
}
IAAT.lm <- do.call(rbind.data.frame, result.IAAT)
names(IAAT.lm) <-
  c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
IAAT.lm$fdr <- p.adjust(IAAT.lm$pval, method = 'fdr')

# MGS vs ASAT
result.ASAT <- list()
for (x in seq(from = 48, to = ncol(progression.model))) {
  print(names(progression.model[x]))
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:47])
  mod.a <-
    lm(mgs ~ ASAT + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['ASAT', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['ASAT', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$ASAT)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.ASAT[[x]] <- result
}
ASAT.lm <- do.call(rbind.data.frame, result.ASAT)
names(ASAT.lm) <-
  c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
ASAT.lm$fdr <- p.adjust(ASAT.lm$pval, method = 'fdr')

# MGS vs TAAT
result.TAAT <- list()
for (x in seq(from = 48, to = ncol(progression.model))) {
  print(names(progression.model[x]))
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:47])
  mod.a <-
    lm(mgs ~ TAAT + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['TAAT', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['TAAT', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$TAAT)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.TAAT[[x]] <- result
}
TAAT.lm <- do.call(rbind.data.frame, result.TAAT)
names(TAAT.lm) <-
  c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
TAAT.lm$fdr <- p.adjust(TAAT.lm$pval, method = 'fdr')

# MGS vs LiverFat..
result.liverfat <- list()
for (x in seq(from = 48, to = ncol(progression.model))) {
  print(names(progression.model[x]))
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:47])
  mod.a <-
    lm(mgs ~ LiverFat.. + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['LiverFat..', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['LiverFat..', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$LiverFat..)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.liverfat[[x]] <- result
}
liverfat.lm <- do.call(rbind.data.frame, result.liverfat)
names(liverfat.lm) <-
  c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
liverfat.lm$fdr <- p.adjust(liverfat.lm$pval, method = 'fdr')

# MGS vs PancreasFat..
result.pancreasfat <- list()
for (x in seq(from = 48, to = ncol(progression.model))) {
  print(names(progression.model[x]))
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:47])
  mod.a <-
    lm(mgs ~ PancreasFat.. + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['PancreasFat..', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['PancreasFat..', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$PancreasFat..)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.pancreasfat[[x]] <- result
}
pancreasfat.lm <- do.call(rbind.data.frame, result.pancreasfat)
names(pancreasfat.lm) <-
  c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
pancreasfat.lm$fdr <- p.adjust(pancreasfat.lm$pval, method = 'fdr')

# MGS vs BP.Sys
result.BP.Sys <- list()
for (x in seq(from = 48, to = ncol(progression.model))) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:47])
  mod.a <-
    lm(mgs ~ BP.S.Mean + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['BP.S.Mean', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['BP.S.Mean', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$BP.S.Mean)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  
  result.BP.Sys[[x]] <- result
}
BP.Sys.lm <- do.call(rbind.data.frame, result.BP.Sys)
names(BP.Sys.lm) <-
  c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
BP.Sys.lm$fdr <- p.adjust(BP.Sys.lm$pval, method = 'fdr')
BP.Sys.lm[BP.Sys.lm$fdr < .1, ]

# MGS vs BP.Dias
result.BP.Dias <- list()
for (x in seq(from = 48, to = ncol(progression.model))) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:47])
  mod.a <-
    lm(mgs ~ BP.D.Mean + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['BP.D.Mean', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['BP.D.Mean', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$BP.D.Mean)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.BP.Dias[[x]] <- result
}
BP.Dias.lm <- do.call(rbind.data.frame, result.BP.Dias)
names(BP.Dias.lm) <-
  c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
BP.Dias.lm$fdr <- p.adjust(BP.Dias.lm$pval, method = 'fdr')
BP.Dias.lm[BP.Dias.lm$fdr < .1, ]

# MGS vs Value.vm.hpf.mean.BL
result.physical <- list()
for (x in seq(from = 48, to = ncol(progression.model))) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:47])
  mod.a <-
    lm(mgs ~ Value.vm.hpf.mean.BL + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Value.vm.hpf.mean.BL', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Value.vm.hpf.mean.BL', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Value.vm.hpf.mean.BL)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.physical[[x]] <- result
}
physical.lm <- do.call(rbind.data.frame, result.physical)
names(physical.lm) <-
  c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
physical.lm$fdr <- p.adjust(physical.lm$pval, method = 'fdr')

# MGS vs HDI
result.hdi<- list()
for (x in seq(from = 48, to = ncol(progression.model))) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:47])
  mod.a <-
    lm(mgs ~ HDI + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['HDI', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['HDI', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$HDI)
  print(cl.delta$estimate)
  print(cl.delta$conf.int[1]) #lower
  print(cl.delta$conf.int[2])
  delta <- cl.delta$estimate
  delta.low <- cl.delta$conf.int[1] #lower
  delta.up <- cl.delta$conf.int[2] #upper
  result <-
    c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
  result.hdi[[x]] <- result
}
hdi.lm <- do.call(rbind.data.frame, result.hdi)
names(hdi.lm) <-
  c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
hdi.lm$fdr <- p.adjust(hdi.lm$pval, method = 'fdr')

# # extract significant ones - MGS ----
# lm.outputs<-grep(".lm",names(.GlobalEnv),value=TRUE)
# results.dfs <- do.call("list",mget(lm.outputs))
# results.dfs.sig <- lapply(results.dfs, function(x){
#   x <- x[x$fdr<.1,]
#   return(data.frame('mgs' = x$mgs, 'delta' = x$delta, 'fdr' = x$fdr))
# })
# mgs.associations <-  do.call(rbind.data.frame, results.dfs.sig)
# mgs.associations$association <- gsub('.lm.[0-9][0-9]', '', row.names(mgs.associations))
# mgs.associations$association <- gsub('.lm.[0-9]', '', mgs.associations$association)
# mgs.associations$association <- gsub('.lm', '', mgs.associations$association)
# 
# spps <- sapply(mgs.associations$mgs, function(x)(
#   tax[grep(x, row.names(tax)),]$species
# ))
# mgs.associations$spps <- spps
# 
# mgs.associations$delta <- as.numeric(mgs.associations$delta)
# direction <- c()
# for (i in 1:nrow(mgs.associations)){
#   if(mgs.associations[i,]$delta > 0) {
#   direction[[i]] <- c("Positive")
# } else {
#   if(mgs.associations[i,]$delta  == 0) {
#     direction[[i]] <- c("Zero")
#   } else {
#     direction[[i]] <- c("Negative")
#   }
# }}
# mgs.associations$direction <- unlist(direction)
# mgs.associations$direction <- factor(mgs.associations$direction, levels=c('Positive', 'Negative'))
# mgs.associations$name <- paste(mgs.associations$mgs, mgs.associations$spps, sep=':')
# 
# mgs.plot <- ggplot(mgs.associations, aes(association, name, fill=delta, shape=direction)) +
#   geom_point(size=3.9)  + scale_shape_manual(values = c(24, 25)) +
#   scale_fill_gradient2(low='blue', high = 'red') +
#   theme(panel.border = element_rect(fill=NA), panel.background = element_blank(),
#         panel.grid.major = element_line(color='gray85', linetype='dashed'),
#         axis.text.x = element_text(angle=90, vjust=0.5, face='bold', color='black', size=12),
#         axis.text.y = element_text(face='bold.italic'),
#         axis.title = element_text(size = 13, face = "bold"), axis.text = element_text(face = "bold"),
#         legend.text = element_text(face = "bold"), legend.title = element_text(face = "bold"),
#         legend.background = element_blank(), legend.box.background = element_blank(), legend.key = element_blank()) +
#   labs(x = NULL, y = NULL, fill = "Cliff's delta", shape = "Direction", title = 'MGS associations')
# 
# ## extract significant ones - Species ----
# lm.outputs<-grep(".lm",names(.GlobalEnv),value=TRUE)
# results.dfs <- do.call("list",mget(lm.outputs))
# results.dfs.sig <- lapply(results.dfs, function(x){
#   x <- x[x$fdr<.1,]
#   return(data.frame('mgs' = x$mgs, 'delta' = x$delta, 'fdr' = x$fdr))
# })
# spps.associations <-  do.call(rbind.data.frame, results.dfs.sig)
# spps.associations$association <- gsub('.lm.[0-9][0-9]', '', row.names(spps.associations))
# spps.associations$association <- gsub('.lm.[0-9]', '', spps.associations$association)
# spps.associations$association <- gsub('.lm', '', spps.associations$association)
# 
# spps <- sapply(spps.associations$mgs, function(x)(
#   tax[grep(x, row.names(tax)),]$species
# ))
# spps.associations$spps <- spps
# 
# spps.associations$delta <- as.numeric(spps.associations$delta)
# direction <- c()
# for (i in 1:nrow(spps.associations)){
#   if(spps.associations[i,]$delta > 0) {
#     direction[[i]] <- c("Positive")
#   } else {
#     if(spps.associations[i,]$delta  == 0) {
#       direction[[i]] <- c("Zero")
#     } else {
#       direction[[i]] <- c("Negative")
#     }
#   }}
# spps.associations$direction <- unlist(direction)
# spps.associations$direction <- factor(spps.associations$direction, levels=c('Positive', 'Negative'))
# spps.associations$name <- paste(spps.associations$mgs, spps.associations$spps, sep=':')
# 
# spps.plot <- ggplot(spps.associations, aes(association, name, fill=delta, shape=direction)) +
#   geom_point(size=3.9)  + scale_shape_manual(values = c(24, 25)) +
#   scale_fill_gradient2(low='blue', high = 'red') +
#   theme(panel.border = element_rect(fill=NA), panel.background = element_blank(),
#         panel.grid.major = element_line(color='gray85', linetype='dashed'),
#         axis.text.x = element_text(angle=90, vjust=0.5, face='bold', color='black', size=12),
#         axis.text.y = element_text(face='bold.italic'),
#         axis.title = element_text(size = 13, face = "bold"), axis.text = element_text(face = "bold"),
#         legend.text = element_text(face = "bold"), legend.title = element_text(face = "bold"),
#         legend.background = element_blank(), legend.box.background = element_blank(), legend.key = element_blank()) +
#   labs(x = NULL, y = NULL, fill = "Cliff's delta", shape = "Direction", title = 'Species associations')
# #
# ## extract significant ones - Genus ----
# lm.outputs<-grep(".lm",names(.GlobalEnv),value=TRUE)
# results.dfs <- do.call("list",mget(lm.outputs))
# results.dfs.sig <- lapply(results.dfs, function(x){
#   x <- x[x$fdr<.1,]
#   return(data.frame('mgs' = x$mgs, 'delta' = x$delta, 'fdr' = x$fdr))
# })
# genus.associations <-  do.call(rbind.data.frame, results.dfs.sig)
# genus.associations$association <- gsub('.lm.[0-9][0-9]', '', row.names(genus.associations))
# genus.associations$association <- gsub('.lm.[0-9]', '', genus.associations$association)
# genus.associations$association <- gsub('.lm', '', genus.associations$association)
# 
# genus <- sapply(genus.associations$mgs, function(x)(
#   tax[grep(x, row.names(tax)),]$genus
# ))
# genus.associations$genus <- genus
# 
# genus.associations$delta <- as.numeric(genus.associations$delta)
# direction <- c()
# for (i in 1:nrow(genus.associations)){
#   if(genus.associations[i,]$delta > 0) {
#     direction[[i]] <- c("Positive")
#   } else {
#     if(genus.associations[i,]$delta  == 0) {
#       direction[[i]] <- c("Zero")
#     } else {
#       direction[[i]] <- c("Negative")
#     }
#   }}
# genus.associations$direction <- unlist(direction)
# genus.associations$direction <- factor(genus.associations$direction, levels=c('Positive', 'Negative'))
# genus.associations$name <- paste(genus.associations$mgs, genus.associations$genus, sep=':')
# 
# genus.plot <- ggplot(genus.associations, aes(association, name, fill=delta, shape=direction)) +
#   geom_point(size=3.9)  + scale_shape_manual(values = c(24, 25)) +
#   scale_fill_gradient2(low='blue', high = 'red') +
#   theme(panel.border = element_rect(fill=NA), panel.background = element_blank(),
#         panel.grid.major = element_line(color='gray85', linetype='dashed'),
#         axis.text.x = element_text(angle=90, vjust=0.5, face='bold', color='black', size=12),
#         axis.text.y = element_text(face='bold.italic'),
#         axis.title = element_text(size = 13, face = "bold"), axis.text = element_text(face = "bold"),
#         legend.text = element_text(face = "bold"), legend.title = element_text(face = "bold"),
#         legend.background = element_blank(), legend.box.background = element_blank(), legend.key = element_blank()) +
#   labs(x = NULL, y = NULL, fill = "Cliff's delta", shape = "Direction", title = 'Genus associations')
# 
## extract significant ones - GMM ----
lm.outputs<-grep(".lm",names(.GlobalEnv),value=TRUE)
results.dfs <- do.call("list",mget(lm.outputs))
results.dfs.sig <- lapply(results.dfs, function(x){
  x <- x[x$fdr<.1,]
  return(data.frame('mgs' = x$mgs, 'delta' = x$delta, 'fdr' = x$fdr))
})
gmm.associations <-  do.call(rbind.data.frame, results.dfs.sig)
gmm.associations$association <- gsub('.lm.[0-9][0-9]', '', row.names(gmm.associations))
gmm.associations$association <- gsub('.lm.[0-9]', '', gmm.associations$association)
gmm.associations$association <- gsub('.lm', '', gmm.associations$association)

gmm.names <- read.table('data/GMMs.v1.07.names', sep='\t')

annotation <- sapply(gmm.associations$mgs, function(x)(
  gmm.names[grep(x, gmm.names$V1),]$V2
))
gmm.associations$annotation <- annotation

gmm.associations$delta <- as.numeric(gmm.associations$delta)
direction <- c()
for (i in 1:nrow(gmm.associations)){
  if(gmm.associations[i,]$delta > 0) {
    direction[[i]] <- c("Positive")
  } else {
    if(gmm.associations[i,]$delta  == 0) {
      direction[[i]] <- c("Zero")
    } else {
      direction[[i]] <- c("Negative")
    }
  }}
gmm.associations$direction <- unlist(direction)
gmm.associations$direction <- factor(gmm.associations$direction, levels=c('Positive', 'Negative'))
gmm.associations$name <- paste(gmm.associations$mgs, gmm.associations$annotation, sep=':')

gmm.plot <- ggplot(gmm.associations, aes(association, name, fill=delta, shape=direction)) +
  geom_point(size=3.9)  + scale_shape_manual(values = c(24, 25)) +
  scale_fill_gradient2(low='blue', high = 'red') +
  theme(panel.border = element_rect(fill=NA), panel.background = element_blank(),
        panel.grid.major = element_line(color='gray85', linetype='dashed'),
        axis.text.x = element_text(angle=90, vjust=0.5, face='bold', color='black', size=12),
        axis.text.y = element_text(face='bold.italic'),
        axis.title = element_text(size = 13, face = "bold"), axis.text = element_text(face = "bold"),
        legend.text = element_text(face = "bold"), legend.title = element_text(face = "bold"),
        legend.background = element_blank(), legend.box.background = element_blank(), legend.key = element_blank()) +
  labs(x = NULL, y = NULL, fill = "Cliff's delta", shape = "Direction", title = 'GMM associations')

## export figures ----
# pdf('baselineMGSassociationstophenotype.pdf', height =18, width = 12)
#   mgs.plot
#   spps.plot
#   genus.plot
#   gmm.plot
# dev.off()

## summarize figure
mgs.associations %>% 
  group_by(association, direction) %>% 
  tally() %>% 
  mutate(n.dir = ifelse(direction=='Negative', n*(-1),n)) %>%
  ggplot(aes(direction, association, size=n, shape = direction, fill = n.dir)) + geom_point() +
  scale_shape_manual(values = c(24, 25)) +
  scale_fill_gradient2(low='blue', high = 'red') +
  # scale_x_continuous(breaks = seq(-30,30, by=5)) +
  # geom_vline(xintercept = 0, linetype='dashed') +
  theme(panel.border = element_rect(fill=NA), panel.background = element_blank(),
        panel.grid.major = element_line(color='gray85', linetype='dashed'),
        axis.text.x = element_text(face='bold', color='black', size=14),
        axis.text.y = element_text(face='bold.italic'),
        axis.title = element_text(size = 13, face = "bold"), axis.text = element_text(face = "bold"),
        legend.text = element_text(face = "bold"), legend.title = element_text(face = "bold"),
        legend.background = element_blank(), legend.box.background = element_blank(), legend.key = element_blank()) +
  labs(x = NULL, y = NULL, fill = "Number of \nassociations", 
       shape = "Direction", title = 'MGS associations') + 
  guides(size=F)

mgs <- mgs.associations %>% 
  group_by(association, direction) %>% 
  tally() %>% 
  mutate(n.dir = ifelse(direction=='Negative', n*(-1),n))# %>%
  # ggplot(aes(n.dir, association, size=n, shape = direction, fill = n.dir)) + geom_point() +
  # scale_shape_manual(values = c(24, 25)) +
  # scale_fill_gradient2(low='blue', high = 'red') +
  # scale_x_continuous(breaks = seq(-30,30, by=5)) +
  # geom_vline(xintercept = 0, linetype='dashed') +
  # theme(panel.border = element_rect(fill=NA), panel.background = element_blank(),
  #       panel.grid.major = element_line(color='gray85', linetype='dashed'),
  #       axis.text.x = element_text(face='bold', color='black'),
  #       axis.text.y = element_text(face='bold.italic'),
  #       axis.title = element_text(size = 13, face = "bold"), axis.text = element_text(face = "bold"),
  #       legend.text = element_text(face = "bold"), legend.title = element_text(face = "bold"),
  #       legend.background = element_blank(), legend.box.background = element_blank(), legend.key = element_blank()) +
  # labs(x = 'Number of associations', y = NULL, fill = "Number of \nassociations",
  #      shape = "Direction", title = 'MGS associations') +
  # guides(size=F)

spps <- spps.associations %>% 
  group_by(association, direction) %>% 
  tally() %>% 
  mutate(n.dir = ifelse(direction=='Negative', n*(-1),n)) #%>%
  # ggplot(aes(n.dir, association, size=n, shape = direction, fill = n.dir)) + geom_point() +
  # scale_shape_manual(values = c(24, 25)) +
  # scale_fill_gradient2(low='blue', high = 'red') +
  # scale_x_continuous(breaks = seq(-30,30, by=5)) +
  # geom_vline(xintercept = 0, linetype='dashed') +
  # theme(panel.border = element_rect(fill=NA), panel.background = element_blank(),
  #       panel.grid.major = element_line(color='gray85', linetype='dashed'),
  #       axis.text.x = element_text(face='bold', color='black'),
  #       axis.text.y = element_text(face='bold.italic'),
  #       axis.title = element_text(size = 13, face = "bold"), axis.text = element_text(face = "bold"),
  #       legend.text = element_text(face = "bold"), legend.title = element_text(face = "bold"),
  #       legend.background = element_blank(), legend.box.background = element_blank(), legend.key = element_blank()) +
  # labs(x = 'Number of associations', y = NULL, fill = "Number of \nassociations", 
  #      shape = "Direction", title = 'Species associations') + 
  # guides(size=F)

genus <- genus.associations %>% 
  group_by(association, direction) %>% 
  tally() %>% 
  mutate(n.dir = ifelse(direction=='Negative', n*(-1),n))# %>%
  # ggplot(aes(n.dir, association, size=n, shape = direction, fill = n.dir)) + geom_point() +
  # scale_shape_manual(values = c(24, 25)) +
  # scale_fill_gradient2(low='blue', high = 'red') +
  # scale_x_continuous(breaks = seq(-30,30, by=5)) +
  # geom_vline(xintercept = 0, linetype='dashed') +
  # theme(panel.border = element_rect(fill=NA), panel.background = element_blank(),
  #       panel.grid.major = element_line(color='gray85', linetype='dashed'),
  #       axis.text.x = element_text(face='bold', color='black'),
  #       axis.text.y = element_text(face='bold.italic'),
  #       axis.title = element_text(size = 13, face = "bold"), axis.text = element_text(face = "bold"),
  #       legend.text = element_text(face = "bold"), legend.title = element_text(face = "bold"),
  #       legend.background = element_blank(), legend.box.background = element_blank(), legend.key = element_blank()) +
  # labs(x = 'Number of associations', y = NULL, fill = "Number of \nassociations", 
  #      shape = "Direction", title = 'Genera associations') + 
  # guides(size=F)

gmm <- gmm.associations %>% 
  group_by(association, direction) %>% 
  tally() %>% 
  mutate(n.dir = ifelse(direction=='Negative', n*(-1),n))# %>%
  # ggplot(aes(n.dir, association, size=n, shape = direction, fill = n.dir)) + geom_point() +
  # scale_shape_manual(values = c(24, 25)) +
  # scale_fill_gradient2(low='blue', high = 'red') +
  # scale_x_continuous(breaks = seq(-30,30, by=5)) +
  # geom_vline(xintercept = 0, linetype='dashed') +
  # theme(panel.border = element_rect(fill=NA), panel.background = element_blank(),
  #       panel.grid.major = element_line(color='gray85', linetype='dashed'),
  #       axis.text.x = element_text(face='bold', color='black'),
  #       axis.text.y = element_text(face='bold.italic'),
  #       axis.title = element_text(size = 13, face = "bold"), axis.text = element_text(face = "bold"),
  #       legend.text = element_text(face = "bold"), legend.title = element_text(face = "bold"),
  #       legend.background = element_blank(), legend.box.background = element_blank(), legend.key = element_blank()) +
  # labs(x = 'Number of associations', y = NULL, fill = "Number of \nassociations", 
  #      shape = "Direction", title = 'GMM associations') + 
  # guides(size=F)

associations.plot <- data.frame(rbind(mgs, spps, genus, gmm), 'variable'=c(rep('MGS', 40), rep('Species', 32), 
                                                      rep('Genera', 24), rep('GMM', 9)))
associations.plot$variable <- factor(associations.plot$variable, levels=c('MGS', 'Species', 'Genera', 'GMM'))

ggplot(associations.plot, aes(n.dir, association, size=n, shape = direction, fill = n.dir)) + geom_point() +
  scale_shape_manual(values = c(24, 25)) +
  scale_fill_gradient2(low='blue', high = 'red') +
  scale_x_continuous(breaks = seq(-30,30, by=5)) +
  geom_vline(xintercept = 0, linetype='dashed') +
  facet_wrap(~variable, scales = 'free', ncol=4) +
  theme(panel.border = element_rect(fill=NA), panel.background = element_blank(),
      panel.grid.major = element_line(color='gray85', linetype='dashed'),
      axis.text.x = element_text(face='bold', color='black'),
      axis.text.y = element_text(face='bold.italic'),
      axis.title = element_text(size = 13, face = "bold"), axis.text = element_text(face = "bold"),
      legend.text = element_text(face = "bold"), legend.title = element_text(face = "bold"),
      legend.background = element_blank(), legend.box.background = element_blank(), legend.key = element_blank(),
      strip.background = element_blank(), strip.text = element_text(face='bold', size=18)) +
  labs(x = 'Number of associations', y = NULL, fill = "Number of \nassociations",
     shape = "Direction") +
  guides(size=F)
