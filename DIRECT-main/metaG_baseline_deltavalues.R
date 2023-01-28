############################################################################################################
######################## Association of baseline metaG to phenotype progression ############################
############################################################################################################

# directory
setwd('/home/scratch/margar')

# libraries
library(effsize)
library(phyloseq)

# data
qmp.bsl <-
  read.table(
    'data/WP2.1_775_IGR_MGS-QMP_igrHA_Ent_M.csv',
    header = T,
    row.names = 1
  )
qmp.bsl <- qmp.bsl[, 7:727]

delta.df <- readRDS('DeltaofLogValues.rds')
qmp.bsl <- qmp.bsl[row.names(qmp.bsl) %in% delta.df$sample, ]
delta.df$sample == row.names(qmp.bsl)
names(qmp.bsl) <- gsub('\\.', ':', names(qmp.bsl))

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

delta.df
delta.df$Gly.Cat.2 <- NULL
delta.df$ADA <- NULL

# MGS modelling -----
progression.model <- data.frame(delta.df, qmp.bsl)

# MGS vs Height
result.height <- list()
for (x in seq(from = 25, to = 745)) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
  mod.a <-
    lm(mgs ~ Height + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Height', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Height', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Height)
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
height.lm[height.lm$fdr < .1, ]
height.lm[height.lm$pval < .05, ]

# MGS vs Weight
result.weight <- list()
for (x in seq(from = 25, to = 745)) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
  mod.a <-
    lm(mgs ~ Weight + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Weight', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Weight', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Weight)
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
weight.lm[weight.lm$fdr < .1, ]
weight.lm[weight.lm$pval < .05, ]

# MGS vs Waist
result.waist <- list()
for (x in seq(from = 25, to = 745)) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
  mod.a <-
    lm(mgs ~ Waist + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Waist', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Waist', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Waist)
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
waist.lm[waist.lm$fdr < .1, ]
# waist.lm[waist.lm$pval<.05,]

# MGS vs BMI
result.bmi <- list()
for (x in seq(from = 25, to = 745)) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
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
bmi.lm[bmi.lm$fdr < .1, ]
# bmi.lm[bmi.lm$pval<.05,]

# MGS vs Fast.Glucose
result.Fast.Glucose <- list()
for (x in seq(from = 25, to = 745)) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
  mod.a <-
    lm(mgs ~ Fast.Glucose + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Fast.Glucose', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Fast.Glucose', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Fast.Glucose)
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
for (x in seq(from = 25, to = 745)) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
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
for (x in seq(from = 25, to = 745)) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
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
for (x in seq(from = 25, to = 745)) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
  mod.a <-
    lm(mgs ~ Glucose.2h + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Glucose.2h', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Glucose.2h', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Glucose.2h)
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
Glucose.2h.lm[Glucose.2h.lm$fdr < .1, ]

# MGS vs Fast.Insulin
result.Fast.Insulin <- list()
for (x in seq(from = 25, to = 745)) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
  mod.a <-
    lm(mgs ~ Fast.Insulin + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Fast.Insulin', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Fast.Insulin', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Fast.Insulin)
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
for (x in seq(from = 25, to = 745)) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
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
for (x in seq(from = 25, to = 745)) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
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
for (x in seq(from = 25, to = 745)) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
  mod.a <-
    lm(mgs ~ Basal.InsSecr + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Basal.InsSecr', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Basal.InsSecr', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Basal.InsSecr)
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
for (x in seq(from = 25, to = 745)) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
  mod.a <-
    lm(mgs ~ Total.InsSecr + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Total.InsSecr', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Total.InsSecr', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Total.InsSecr)
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
for (x in seq(from = 25, to = 745)) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
  mod.a <-
    lm(
      mgs ~ Gluc.Sensitivity + Gender + CenterID + Age,
      data = model,
      na.action = na.omit
    )
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Gluc.Sensitivity', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Gluc.Sensitivity', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Gluc.Sensitivity)
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
for (x in seq(from = 25, to = 745)) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
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

# MGS vs Matsuda
result.Matsuda <- list()
for (x in seq(from = 25, to = 745)) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
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

# MGS vs BP.Sys
result.BP.Sys <- list()
for (x in seq(from = 25, to = 745)) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
  mod.a <-
    lm(mgs ~ BP.Sys + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['BP.Sys', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['BP.Sys', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$BP.Sys)
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
for (x in seq(from = 25, to = 745)) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
  mod.a <-
    lm(mgs ~ BP.Dias + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['BP.Dias', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['BP.Dias', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$BP.Dias)
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

# save results
write.table(height.lm, 'MGS_delta_lm_height.csv', sep = ';')
write.table(weight.lm, 'MGS_delta_lm_weight.csv', sep = ';')
write.table(waist.lm, 'MGS_delta_lm_waist.csv', sep = ';')
write.table(bmi.lm, 'MGS_delta_lm_BMI.csv', sep = ';')
write.table(Fast.Glucose.lm, 'MGS_delta_lm_fastGluc.csv', sep = ';')
write.table(Basal.Glucose.lm, 'MGS_delta_lm_basGluc.csv', sep = ';')
write.table(Mean.Glucose.lm, 'MGS_delta_lm_meanGluc.csv', sep = ';')
write.table(Glucose.2h.lm, 'MGS_delta_lm_2hGluc.csv', sep = ';')
write.table(Fast.Insulin.lm, 'MGS_delta_lm_fasIns.csv', sep = ';')
write.table(Basal.Insulin.lm, 'MGS_delta_lm_basIns.csv', sep = ';')
write.table(Mean.Insulin.lm, 'MGS_delta_lm_meanIns.csv', sep = ';')
write.table(Basal.InsSecr.lm, 'MGS_delta_lm_basInsSecr.csv', sep = ';')
write.table(Total.InsSecr.lm, 'MGS_delta_lm_totalInsSecr.csv', sep = ';')
write.table(Gluc.Sensitivity.lm, 'MGS_delta_lm_GlucSens.csv', sep = ';')
write.table(OGIS.2h.lm, 'MGS_delta_lm_OGIS2h.csv', sep = ';')
write.table(Matsuda.lm, 'MGS_delta_lm_Matsuda.csv', sep = ';')
write.table(BP.Sys.lm, 'MGS_delta_lm_BPSys.csv', sep = ';')
write.table(BP.Dias.lm, 'MGS_delta_lm_BPDias.csv', sep = ';')
tax.1 <- tax[row.names(tax) %in% names(qmp.bsl), ]


# species modelling ----
physeq.spps <- tax_glom(physeq, taxrank = rank_names(physeq)[7])
spps <- as.data.frame(physeq.spps@otu_table)

progression.model <- data.frame(delta.df, spps)

# MGS vs Height
result.height <- list()
for (x in seq(from = 25, to = 276)) {
  print(names(progression.model[x]))
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
  mod.a <-
    lm(mgs ~ Height + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Height', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Height', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Height)
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
height.lm[height.lm$fdr < .1, ]
height.lm[height.lm$pval < .05, ]

# MGS vs Weight
result.weight <- list()
for (x in seq(from = 25, to = 276)) {
  print(names(progression.model[x]))
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
  mod.a <-
    lm(mgs ~ Weight + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Weight', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Weight', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Weight)
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
weight.lm[weight.lm$fdr < .1, ]
weight.lm[weight.lm$pval < .05, ]

# MGS vs Waist
result.waist <- list()
for (x in seq(from = 25, to = 276)) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
  mod.a <-
    lm(mgs ~ Waist + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Waist', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Waist', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Waist)
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
waist.lm[waist.lm$fdr < .1, ]
waist.lm[waist.lm$pval<.05,]

# MGS vs BMI
result.bmi <- list()
for (x in seq(from = 25, to = 276)) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
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
bmi.lm[bmi.lm$fdr < .1, ]
bmi.lm[bmi.lm$pval<.05,]

# MGS vs Fast.Glucose
result.Fast.Glucose <- list()
for (x in seq(from = 25, to = 276)) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
  mod.a <-
    lm(mgs ~ Fast.Glucose + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Fast.Glucose', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Fast.Glucose', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Fast.Glucose)
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
for (x in seq(from = 25, to = 276)) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
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
for (x in seq(from = 25, to = 276)) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
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
for (x in seq(from = 25, to = 276)) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
  mod.a <-
    lm(mgs ~ Glucose.2h + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Glucose.2h', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Glucose.2h', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Glucose.2h)
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
Glucose.2h.lm[Glucose.2h.lm$fdr < .1, ]

# MGS vs Fast.Insulin
result.Fast.Insulin <- list()
for (x in seq(from = 25, to = 276)) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
  mod.a <-
    lm(mgs ~ Fast.Insulin + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Fast.Insulin', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Fast.Insulin', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Fast.Insulin)
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
for (x in seq(from = 25, to = 276)) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
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
for (x in seq(from = 25, to = 276)) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
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
for (x in seq(from = 25, to = 276)) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
  mod.a <-
    lm(mgs ~ Basal.InsSecr + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Basal.InsSecr', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Basal.InsSecr', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Basal.InsSecr)
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
for (x in seq(from = 25, to = 276)) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
  mod.a <-
    lm(mgs ~ Total.InsSecr + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Total.InsSecr', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Total.InsSecr', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Total.InsSecr)
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
for (x in seq(from = 25, to = 276)) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
  mod.a <-
    lm(
      mgs ~ Gluc.Sensitivity + Gender + CenterID + Age,
      data = model,
      na.action = na.omit
    )
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['Gluc.Sensitivity', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['Gluc.Sensitivity', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$Gluc.Sensitivity)
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
for (x in seq(from = 25, to = 276)) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
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

# MGS vs Matsuda
result.Matsuda <- list()
for (x in seq(from = 25, to = 276)) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
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

# MGS vs BP.Sys
result.BP.Sys <- list()
for (x in seq(from = 25, to = 276)) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
  mod.a <-
    lm(mgs ~ BP.Sys + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['BP.Sys', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['BP.Sys', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$BP.Sys)
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
for (x in seq(from = 25, to = 276)) {
  print(names(progression.model[x]))
  # result.model$MGS[x] <- names(mgs.model.bsl[x])
  model <-
    data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
  mod.a <-
    lm(mgs ~ BP.Dias + Gender + CenterID + Age,
       data = model,
       na.action = na.omit)
  summ.mod.a <- summary(mod.a)
  print(summ.mod.a$coefficients['BP.Dias', 'Pr(>|t|)'])
  p.val <- summ.mod.a$coefficients['BP.Dias', 'Pr(>|t|)']
  cl.delta <- cliff.delta(model$mgs, model$BP.Dias)
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

# # save results
# write.table(height.lm, 'SPPS_delta_lm_height.csv', sep = ';')
# write.table(weight.lm, 'SPPS_delta_lm_weight.csv', sep = ';')
# write.table(waist.lm, 'SPPS_delta_lm_waist.csv', sep = ';')
# write.table(bmi.lm, 'SPPS_delta_lm_BMI.csv', sep = ';')
# write.table(Fast.Glucose.lm, 'SPPS_delta_lm_fastGluc.csv', sep = ';')
# write.table(Basal.Glucose.lm, 'SPPS_delta_lm_basGluc.csv', sep = ';')
# write.table(Mean.Glucose.lm, 'SPPS_delta_lm_meanGluc.csv', sep = ';')
# write.table(Glucose.2h.lm, 'SPPS_delta_lm_2hGluc.csv', sep = ';')
# write.table(Fast.Insulin.lm, 'SPPS_delta_lm_fasIns.csv', sep = ';')
# write.table(Basal.Insulin.lm, 'SPPS_delta_lm_basIns.csv', sep = ';')
# write.table(Mean.Insulin.lm, 'SPPS_delta_lm_meanIns.csv', sep = ';')
# write.table(Basal.InsSecr.lm, 'SPPS_delta_lm_basInsSecr.csv', sep = ';')
# write.table(Total.InsSecr.lm, 'SPPS_delta_lm_totalInsSecr.csv', sep = ';')
# write.table(Gluc.Sensitivity.lm, 'SPPS_delta_lm_GlucSens.csv', sep = ';')
# write.table(OGIS.2h.lm, 'SPPS_delta_lm_OGIS2h.csv', sep = ';')
# write.table(Matsuda.lm, 'SPPS_delta_lm_Matsuda.csv', sep = ';')
# write.table(BP.Sys.lm, 'SPPS_delta_lm_BPSys.csv', sep = ';')
# write.table(BP.Dias.lm, 'SPPS_delta_lm_BPDias.csv', sep = ';')
# tax.1 <- tax[row.names(tax) %in% names(spps), ]
# write.table(tax.1, 'SPPS_delta_lm_TAXONOMY.csv', sep = ';')


## extract significant ones - Species ----
lm.outputs<-grep(".lm",names(.GlobalEnv),value=TRUE)
results.dfs <- do.call("list",mget(lm.outputs))
results.dfs.sig <- lapply(results.dfs, function(x){
  x <- x[x$fdr<.1,]
  return(data.frame('mgs' = x$mgs, 'delta' = x$delta, 'fdr' = x$fdr))
})
spps.associations <-  do.call(rbind.data.frame, results.dfs.sig)
spps.associations$association <- gsub('.lm.[0-9][0-9]', '', row.names(spps.associations))
spps.associations$association <- gsub('.lm.[0-9]', '', spps.associations$association)
spps.associations$association <- gsub('.lm', '', spps.associations$association)

spps <- sapply(spps.associations$mgs, function(x)(
  tax[grep(x, row.names(tax)),]$species
))
spps.associations$spps <- spps

spps.associations$delta <- as.numeric(spps.associations$delta)
direction <- c()
for (i in 1:nrow(spps.associations)){
  if(spps.associations[i,]$delta > 0) {
    direction[[i]] <- c("Positive")
  } else {
    if(spps.associations[i,]$delta  == 0) {
      direction[[i]] <- c("Zero")
    } else {
      direction[[i]] <- c("Negative")
    }
  }}
spps.associations$direction <- unlist(direction)
spps.associations$direction <- factor(spps.associations$direction, levels=c('Positive', 'Negative'))
spps.associations$name <- paste(spps.associations$mgs, spps.associations$spps, sep=':')

spps.plot <- ggplot(spps.associations, aes(association, name, fill=delta, shape=direction)) +
  geom_point(size=3.9)  + scale_shape_manual(values = c(24, 25)) +
  scale_fill_gradient2(low='blue', high = 'red') +
  theme(panel.border = element_rect(fill=NA), panel.background = element_blank(),
        panel.grid.major = element_line(color='gray85', linetype='dashed'),
        axis.text.x = element_text(angle=90, vjust=0.5, face='bold', color='black', size=12),
        axis.text.y = element_text(face='bold.italic'),
        axis.title = element_text(size = 13, face = "bold"), axis.text = element_text(face = "bold"),
        legend.text = element_text(face = "bold"), legend.title = element_text(face = "bold"),
        legend.background = element_blank(), legend.box.background = element_blank(), legend.key = element_blank()) +
  labs(x = NULL, y = NULL, fill = "Cliff's delta", shape = "Direction", title = 'Species associations')
#



# # genus modelling ----
# physeq.gen <- tax_glom(physeq, taxrank = rank_names(physeq)[6])
# genus <- as.data.frame(physeq.gen@otu_table)
# 
# progression.model <- data.frame(delta.df, genus)
# 
# # MGS vs Height
# result.height <- list()
# for (x in seq(from = 25, to = 147)) {
#   print(names(progression.model[x]))
#   model <-
#     data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
#   mod.a <-
#     lm(mgs ~ Height + Gender + CenterID + Age,
#        data = model,
#        na.action = na.omit)
#   summ.mod.a <- summary(mod.a)
#   print(summ.mod.a$coefficients['Height', 'Pr(>|t|)'])
#   p.val <- summ.mod.a$coefficients['Height', 'Pr(>|t|)']
#   cl.delta <- cliff.delta(model$mgs, model$Height)
#   print(cl.delta$estimate)
#   print(cl.delta$conf.int[1]) #lower
#   print(cl.delta$conf.int[2])
#   delta <- cl.delta$estimate
#   delta.low <- cl.delta$conf.int[1] #lower
#   delta.up <- cl.delta$conf.int[2] #upper
#   result <-
#     c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
#   result.height[[x]] <- result
# }
# 
# height.lm <- do.call(rbind.data.frame, result.height)
# names(height.lm) <-
#   c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
# height.lm$fdr <- p.adjust(height.lm$pval, method = 'fdr')
# height.lm[height.lm$fdr < .1, ]
# height.lm[height.lm$pval < .05, ]
# 
# # MGS vs Weight
# result.weight <- list()
# for (x in seq(from = 25, to = 147)) {
#   print(names(progression.model[x]))
#   model <-
#     data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
#   mod.a <-
#     lm(mgs ~ Weight + Gender + CenterID + Age,
#        data = model,
#        na.action = na.omit)
#   summ.mod.a <- summary(mod.a)
#   print(summ.mod.a$coefficients['Weight', 'Pr(>|t|)'])
#   p.val <- summ.mod.a$coefficients['Weight', 'Pr(>|t|)']
#   cl.delta <- cliff.delta(model$mgs, model$Weight)
#   print(cl.delta$estimate)
#   print(cl.delta$conf.int[1]) #lower
#   print(cl.delta$conf.int[2])
#   delta <- cl.delta$estimate
#   delta.low <- cl.delta$conf.int[1] #lower
#   delta.up <- cl.delta$conf.int[2] #upper
#   result <-
#     c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
#   result.weight[[x]] <- result
# }
# 
# weight.lm <- do.call(rbind.data.frame, result.weight)
# names(weight.lm) <-
#   c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
# weight.lm$fdr <- p.adjust(weight.lm$pval, method = 'fdr')
# weight.lm[weight.lm$fdr < .1, ]
# weight.lm[weight.lm$pval < .05, ]
# 
# # MGS vs Waist
# result.waist <- list()
# for (x in seq(from = 25, to = 147)) {
#   print(names(progression.model[x]))
#   # result.model$MGS[x] <- names(mgs.model.bsl[x])
#   model <-
#     data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
#   mod.a <-
#     lm(mgs ~ Waist + Gender + CenterID + Age,
#        data = model,
#        na.action = na.omit)
#   summ.mod.a <- summary(mod.a)
#   print(summ.mod.a$coefficients['Waist', 'Pr(>|t|)'])
#   p.val <- summ.mod.a$coefficients['Waist', 'Pr(>|t|)']
#   cl.delta <- cliff.delta(model$mgs, model$Waist)
#   print(cl.delta$estimate)
#   print(cl.delta$conf.int[1]) #lower
#   print(cl.delta$conf.int[2])
#   delta <- cl.delta$estimate
#   delta.low <- cl.delta$conf.int[1] #lower
#   delta.up <- cl.delta$conf.int[2] #upper
#   result <-
#     c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
#   
#   result.waist[[x]] <- result
# }
# 
# waist.lm <- do.call(rbind.data.frame, result.waist)
# names(waist.lm) <-
#   c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
# waist.lm$fdr <- p.adjust(waist.lm$pval, method = 'fdr')
# waist.lm[waist.lm$fdr < .1, ]
# waist.lm[waist.lm$pval<.05,]
# 
# # MGS vs BMI
# result.bmi <- list()
# for (x in seq(from = 25, to = 147)) {
#   print(names(progression.model[x]))
#   # result.model$MGS[x] <- names(mgs.model.bsl[x])
#   model <-
#     data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
#   mod.a <-
#     lm(mgs ~ BMI + Gender + CenterID + Age,
#        data = model,
#        na.action = na.omit)
#   summ.mod.a <- summary(mod.a)
#   print(summ.mod.a$coefficients['BMI', 'Pr(>|t|)'])
#   p.val <- summ.mod.a$coefficients['BMI', 'Pr(>|t|)']
#   cl.delta <- cliff.delta(model$mgs, model$BMI)
#   print(cl.delta$estimate)
#   print(cl.delta$conf.int[1]) #lower
#   print(cl.delta$conf.int[2])
#   delta <- cl.delta$estimate
#   delta.low <- cl.delta$conf.int[1] #lower
#   delta.up <- cl.delta$conf.int[2] #upper
#   result <-
#     c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
#   result.bmi[[x]] <- result
# }
# 
# bmi.lm <- do.call(rbind.data.frame, result.bmi)
# names(bmi.lm) <- c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
# bmi.lm$fdr <- p.adjust(bmi.lm$pval, method = 'fdr')
# bmi.lm[bmi.lm$fdr < .1, ]
# bmi.lm[bmi.lm$pval<.05,]
# 
# # MGS vs Fast.Glucose
# result.Fast.Glucose <- list()
# for (x in seq(from = 25, to = 147)) {
#   print(names(progression.model[x]))
#   # result.model$MGS[x] <- names(mgs.model.bsl[x])
#   model <-
#     data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
#   mod.a <-
#     lm(mgs ~ Fast.Glucose + Gender + CenterID + Age,
#        data = model,
#        na.action = na.omit)
#   summ.mod.a <- summary(mod.a)
#   print(summ.mod.a$coefficients['Fast.Glucose', 'Pr(>|t|)'])
#   p.val <- summ.mod.a$coefficients['Fast.Glucose', 'Pr(>|t|)']
#   cl.delta <- cliff.delta(model$mgs, model$Fast.Glucose)
#   print(cl.delta$estimate)
#   print(cl.delta$conf.int[1]) #lower
#   print(cl.delta$conf.int[2])
#   delta <- cl.delta$estimate
#   delta.low <- cl.delta$conf.int[1] #lower
#   
#   delta.up <- cl.delta$conf.int[2] #upper
#   result <-
#     c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
#   result.Fast.Glucose[[x]] <- result
# }
# Fast.Glucose.lm <- do.call(rbind.data.frame, result.Fast.Glucose)
# names(Fast.Glucose.lm) <-
#   c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
# Fast.Glucose.lm$fdr <- p.adjust(Fast.Glucose.lm$pval, method = 'fdr')
# Fast.Glucose.lm[Fast.Glucose.lm$fdr < .1, ]
# 
# # MGS vs Basal.Glucose
# result.Basal.Glucose <- list()
# for (x in seq(from = 25, to = 147)) {
#   print(names(progression.model[x]))
#   # result.model$MGS[x] <- names(mgs.model.bsl[x])
#   model <-
#     data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
#   mod.a <-
#     lm(mgs ~ Basal.Glucose + Gender + CenterID + Age,
#        data = model,
#        na.action = na.omit)
#   summ.mod.a <- summary(mod.a)
#   print(summ.mod.a$coefficients['Basal.Glucose', 'Pr(>|t|)'])
#   p.val <- summ.mod.a$coefficients['Basal.Glucose', 'Pr(>|t|)']
#   cl.delta <- cliff.delta(model$mgs, model$Basal.Glucose)
#   print(cl.delta$estimate)
#   print(cl.delta$conf.int[1]) #lower
#   print(cl.delta$conf.int[2])
#   delta <- cl.delta$estimate
#   delta.low <- cl.delta$conf.int[1] #lower
#   delta.up <- cl.delta$conf.int[2] #upper
#   result <-
#     c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
#   result.Basal.Glucose[[x]] <- result
# }
# Basal.Glucose.lm <- do.call(rbind.data.frame, result.Basal.Glucose)
# names(Basal.Glucose.lm) <-
#   c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
# Basal.Glucose.lm$fdr <-
#   p.adjust(Basal.Glucose.lm$pval, method = 'fdr')
# Basal.Glucose.lm[Basal.Glucose.lm$fdr < .1, ]
# 
# # MGS vs Mean.Glucose
# result.Mean.Glucose <- list()
# for (x in seq(from = 25, to = 147)) {
#   print(names(progression.model[x]))
#   # result.model$MGS[x] <- names(mgs.model.bsl[x])
#   model <-
#     data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
#   mod.a <-
#     lm(mgs ~ Mean.Glucose + Gender + CenterID + Age,
#        data = model,
#        na.action = na.omit)
#   summ.mod.a <- summary(mod.a)
#   print(summ.mod.a$coefficients['Mean.Glucose', 'Pr(>|t|)'])
#   p.val <- summ.mod.a$coefficients['Mean.Glucose', 'Pr(>|t|)']
#   cl.delta <- cliff.delta(model$mgs, model$Mean.Glucose)
#   print(cl.delta$estimate)
#   print(cl.delta$conf.int[1]) #lower
#   print(cl.delta$conf.int[2])
#   delta <- cl.delta$estimate
#   delta.low <- cl.delta$conf.int[1] #lower
#   delta.up <- cl.delta$conf.int[2] #upper
#   result <-
#     c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
#   result.Mean.Glucose[[x]] <- result
# }
# Mean.Glucose.lm <- do.call(rbind.data.frame, result.Mean.Glucose)
# names(Mean.Glucose.lm) <-
#   c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
# Mean.Glucose.lm$fdr <- p.adjust(Mean.Glucose.lm$pval, method = 'fdr')
# Mean.Glucose.lm[Mean.Glucose.lm$fdr < .1, ]
# 
# # MGS vs Glucose.2h
# result.Glucose.2h <- list()
# for (x in seq(from = 25, to = 147)) {
#   print(names(progression.model[x]))
#   # result.model$MGS[x] <- names(mgs.model.bsl[x])
#   model <-
#     data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
#   mod.a <-
#     lm(mgs ~ Glucose.2h + Gender + CenterID + Age,
#        data = model,
#        na.action = na.omit)
#   summ.mod.a <- summary(mod.a)
#   print(summ.mod.a$coefficients['Glucose.2h', 'Pr(>|t|)'])
#   p.val <- summ.mod.a$coefficients['Glucose.2h', 'Pr(>|t|)']
#   cl.delta <- cliff.delta(model$mgs, model$Glucose.2h)
#   print(cl.delta$estimate)
#   print(cl.delta$conf.int[1]) #lower
#   print(cl.delta$conf.int[2])
#   delta <- cl.delta$estimate
#   delta.low <- cl.delta$conf.int[1] #lower
#   delta.up <- cl.delta$conf.int[2] #upper
#   result <-
#     c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
#   result.Glucose.2h[[x]] <- result
# }
# Glucose.2h.lm <- do.call(rbind.data.frame, result.Glucose.2h)
# names(Glucose.2h.lm) <-
#   c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
# Glucose.2h.lm$fdr <- p.adjust(Glucose.2h.lm$pval, method = 'fdr')
# Glucose.2h.lm[Glucose.2h.lm$fdr < .1, ]
# 
# # MGS vs Fast.Insulin
# result.Fast.Insulin <- list()
# for (x in seq(from = 25, to = 147)) {
#   print(names(progression.model[x]))
#   # result.model$MGS[x] <- names(mgs.model.bsl[x])
#   model <-
#     data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
#   mod.a <-
#     lm(mgs ~ Fast.Insulin + Gender + CenterID + Age,
#        data = model,
#        na.action = na.omit)
#   summ.mod.a <- summary(mod.a)
#   print(summ.mod.a$coefficients['Fast.Insulin', 'Pr(>|t|)'])
#   p.val <- summ.mod.a$coefficients['Fast.Insulin', 'Pr(>|t|)']
#   cl.delta <- cliff.delta(model$mgs, model$Fast.Insulin)
#   print(cl.delta$estimate)
#   print(cl.delta$conf.int[1]) #lower
#   print(cl.delta$conf.int[2])
#   delta <- cl.delta$estimate
#   delta.low <- cl.delta$conf.int[1] #lower
#   
#   delta.up <- cl.delta$conf.int[2] #upper
#   result <-
#     c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
#   result.Fast.Insulin[[x]] <- result
# }
# Fast.Insulin.lm <- do.call(rbind.data.frame, result.Fast.Insulin)
# names(Fast.Insulin.lm) <-
#   c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
# Fast.Insulin.lm$fdr <- p.adjust(Fast.Insulin.lm$pval, method = 'fdr')
# Fast.Insulin.lm[Fast.Insulin.lm$fdr < .1, ]
# 
# # MGS vs Basal.Insulin
# result.Basal.Insulin <- list()
# for (x in seq(from = 25, to = 147)) {
#   print(names(progression.model[x]))
#   # result.model$MGS[x] <- names(mgs.model.bsl[x])
#   model <-
#     data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
#   mod.a <-
#     lm(mgs ~ Basal.Insulin + Gender + CenterID + Age,
#        data = model,
#        na.action = na.omit)
#   summ.mod.a <- summary(mod.a)
#   print(summ.mod.a$coefficients['Basal.Insulin', 'Pr(>|t|)'])
#   p.val <- summ.mod.a$coefficients['Basal.Insulin', 'Pr(>|t|)']
#   cl.delta <- cliff.delta(model$mgs, model$Basal.Insulin)
#   print(cl.delta$estimate)
#   print(cl.delta$conf.int[1]) #lower
#   print(cl.delta$conf.int[2])
#   delta <- cl.delta$estimate
#   delta.low <- cl.delta$conf.int[1] #lower
#   delta.up <- cl.delta$conf.int[2] #upper
#   result <-
#     c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
#   result.Basal.Insulin[[x]] <- result
# }
# Basal.Insulin.lm <- do.call(rbind.data.frame, result.Basal.Insulin)
# names(Basal.Insulin.lm) <-
#   c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
# Basal.Insulin.lm$fdr <-
#   p.adjust(Basal.Insulin.lm$pval, method = 'fdr')
# Basal.Insulin.lm[Basal.Insulin.lm$fdr < .1, ]
# 
# # MGS vs Mean.Insulin
# result.Mean.Insulin <- list()
# for (x in seq(from = 25, to = 147)) {
#   print(names(progression.model[x]))
#   # result.model$MGS[x] <- names(mgs.model.bsl[x])
#   model <-
#     data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
#   mod.a <-
#     lm(mgs ~ Mean.Insulin + Gender + CenterID + Age,
#        data = model,
#        na.action = na.omit)
#   summ.mod.a <- summary(mod.a)
#   print(summ.mod.a$coefficients['Mean.Insulin', 'Pr(>|t|)'])
#   p.val <- summ.mod.a$coefficients['Mean.Insulin', 'Pr(>|t|)']
#   cl.delta <- cliff.delta(model$mgs, model$Mean.Insulin)
#   print(cl.delta$estimate)
#   print(cl.delta$conf.int[1]) #lower
#   print(cl.delta$conf.int[2])
#   delta <- cl.delta$estimate
#   delta.low <- cl.delta$conf.int[1] #lower
#   delta.up <- cl.delta$conf.int[2] #upper
#   result <-
#     c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
#   
#   result.Mean.Insulin[[x]] <- result
# }
# Mean.Insulin.lm <- do.call(rbind.data.frame, result.Mean.Insulin)
# names(Mean.Insulin.lm) <-
#   c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
# Mean.Insulin.lm$fdr <- p.adjust(Mean.Insulin.lm$pval, method = 'fdr')
# Mean.Insulin.lm[Mean.Insulin.lm$fdr < .1, ]
# 
# # MGS vs Basal.InsSecr
# result.Basal.InsSecr <- list()
# for (x in seq(from = 25, to = 147)) {
#   print(names(progression.model[x]))
#   # result.model$MGS[x] <- names(mgs.model.bsl[x])
#   model <-
#     data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
#   mod.a <-
#     lm(mgs ~ Basal.InsSecr + Gender + CenterID + Age,
#        data = model,
#        na.action = na.omit)
#   summ.mod.a <- summary(mod.a)
#   print(summ.mod.a$coefficients['Basal.InsSecr', 'Pr(>|t|)'])
#   p.val <- summ.mod.a$coefficients['Basal.InsSecr', 'Pr(>|t|)']
#   cl.delta <- cliff.delta(model$mgs, model$Basal.InsSecr)
#   print(cl.delta$estimate)
#   print(cl.delta$conf.int[1]) #lower
#   print(cl.delta$conf.int[2])
#   delta <- cl.delta$estimate
#   delta.low <- cl.delta$conf.int[1] #lower
#   delta.up <- cl.delta$conf.int[2] #upper
#   result <-
#     c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
#   result.Basal.InsSecr[[x]] <- result
# }
# Basal.InsSecr.lm <- do.call(rbind.data.frame, result.Basal.InsSecr)
# names(Basal.InsSecr.lm) <-
#   c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
# Basal.InsSecr.lm$fdr <-
#   p.adjust(Basal.InsSecr.lm$pval, method = 'fdr')
# Basal.InsSecr.lm[Basal.InsSecr.lm$fdr < .1, ]
# 
# # MGS vs Total.InsSecr
# result.Total.InsSecr <- list()
# for (x in seq(from = 25, to = 147)) {
#   print(names(progression.model[x]))
#   # result.model$MGS[x] <- names(mgs.model.bsl[x])
#   model <-
#     data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
#   mod.a <-
#     lm(mgs ~ Total.InsSecr + Gender + CenterID + Age,
#        data = model,
#        na.action = na.omit)
#   summ.mod.a <- summary(mod.a)
#   print(summ.mod.a$coefficients['Total.InsSecr', 'Pr(>|t|)'])
#   p.val <- summ.mod.a$coefficients['Total.InsSecr', 'Pr(>|t|)']
#   cl.delta <- cliff.delta(model$mgs, model$Total.InsSecr)
#   print(cl.delta$estimate)
#   print(cl.delta$conf.int[1]) #lower
#   print(cl.delta$conf.int[2])
#   delta <- cl.delta$estimate
#   delta.low <- cl.delta$conf.int[1] #lower
#   
#   delta.up <- cl.delta$conf.int[2] #upper
#   result <-
#     c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
#   result.Total.InsSecr[[x]] <- result
# }
# Total.InsSecr.lm <- do.call(rbind.data.frame, result.Total.InsSecr)
# names(Total.InsSecr.lm) <-
#   c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
# Total.InsSecr.lm$fdr <-
#   p.adjust(Total.InsSecr.lm$pval, method = 'fdr')
# Total.InsSecr.lm[Total.InsSecr.lm$fdr < .1, ]
# 
# # MGS vs Gluc.Sensitivity
# result.Gluc.Sensitivity <- list()
# for (x in seq(from = 25, to = 147)) {
#   print(names(progression.model[x]))
#   # result.model$MGS[x] <- names(mgs.model.bsl[x])
#   model <-
#     data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
#   mod.a <-
#     lm(
#       mgs ~ Gluc.Sensitivity + Gender + CenterID + Age,
#       data = model,
#       na.action = na.omit
#     )
#   summ.mod.a <- summary(mod.a)
#   print(summ.mod.a$coefficients['Gluc.Sensitivity', 'Pr(>|t|)'])
#   p.val <- summ.mod.a$coefficients['Gluc.Sensitivity', 'Pr(>|t|)']
#   cl.delta <- cliff.delta(model$mgs, model$Gluc.Sensitivity)
#   print(cl.delta$estimate)
#   print(cl.delta$conf.int[1]) #lower
#   print(cl.delta$conf.int[2])
#   delta <- cl.delta$estimate
#   delta.low <- cl.delta$conf.int[1] #lower
#   delta.up <- cl.delta$conf.int[2] #upper
#   result <-
#     c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
#   result.Gluc.Sensitivity[[x]] <- result
# }
# Gluc.Sensitivity.lm <-
#   do.call(rbind.data.frame, result.Gluc.Sensitivity)
# names(Gluc.Sensitivity.lm) <-
#   c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
# Gluc.Sensitivity.lm$fdr <-
#   p.adjust(Gluc.Sensitivity.lm$pval, method = 'fdr')
# Gluc.Sensitivity.lm[Gluc.Sensitivity.lm$fdr < .1, ]
# 
# # MGS vs OGIS.2h
# result.OGIS.2h <- list()
# for (x in seq(from = 25, to = 147)) {
#   print(names(progression.model[x]))
#   # result.model$MGS[x] <- names(mgs.model.bsl[x])
#   model <-
#     data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
#   mod.a <-
#     lm(mgs ~ OGIS.2h + Gender + CenterID + Age,
#        data = model,
#        na.action = na.omit)
#   summ.mod.a <- summary(mod.a)
#   print(summ.mod.a$coefficients['OGIS.2h', 'Pr(>|t|)'])
#   p.val <- summ.mod.a$coefficients['OGIS.2h', 'Pr(>|t|)']
#   cl.delta <- cliff.delta(model$mgs, model$OGIS.2h)
#   print(cl.delta$estimate)
#   print(cl.delta$conf.int[1]) #lower
#   print(cl.delta$conf.int[2])
#   delta <- cl.delta$estimate
#   delta.low <- cl.delta$conf.int[1] #lower
#   delta.up <- cl.delta$conf.int[2] #upper
#   result <-
#     c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
#   
#   result.OGIS.2h[[x]] <- result
# }
# OGIS.2h.lm <- do.call(rbind.data.frame, result.OGIS.2h)
# names(OGIS.2h.lm) <-
#   c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
# OGIS.2h.lm$fdr <- p.adjust(OGIS.2h.lm$pval, method = 'fdr')
# OGIS.2h.lm[OGIS.2h.lm$fdr < .1, ]
# 
# # MGS vs Matsuda
# result.Matsuda <- list()
# for (x in seq(from = 25, to = 147)) {
#   print(names(progression.model[x]))
#   # result.model$MGS[x] <- names(mgs.model.bsl[x])
#   model <-
#     data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
#   mod.a <-
#     lm(mgs ~ Matsuda + Gender + CenterID + Age,
#        data = model,
#        na.action = na.omit)
#   summ.mod.a <- summary(mod.a)
#   print(summ.mod.a$coefficients['Matsuda', 'Pr(>|t|)'])
#   p.val <- summ.mod.a$coefficients['Matsuda', 'Pr(>|t|)']
#   cl.delta <- cliff.delta(model$mgs, model$Matsuda)
#   print(cl.delta$estimate)
#   print(cl.delta$conf.int[1]) #lower
#   print(cl.delta$conf.int[2])
#   delta <- cl.delta$estimate
#   delta.low <- cl.delta$conf.int[1] #lower
#   delta.up <- cl.delta$conf.int[2] #upper
#   result <-
#     c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
#   result.Matsuda[[x]] <- result
# }
# Matsuda.lm <- do.call(rbind.data.frame, result.Matsuda)
# names(Matsuda.lm) <-
#   c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
# Matsuda.lm$fdr <- p.adjust(Matsuda.lm$pval, method = 'fdr')
# Matsuda.lm[Matsuda.lm$fdr < .1, ]
# 
# # MGS vs BP.Sys
# result.BP.Sys <- list()
# for (x in seq(from = 25, to = 147)) {
#   print(names(progression.model[x]))
#   # result.model$MGS[x] <- names(mgs.model.bsl[x])
#   model <-
#     data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
#   mod.a <-
#     lm(mgs ~ BP.Sys + Gender + CenterID + Age,
#        data = model,
#        na.action = na.omit)
#   summ.mod.a <- summary(mod.a)
#   print(summ.mod.a$coefficients['BP.Sys', 'Pr(>|t|)'])
#   p.val <- summ.mod.a$coefficients['BP.Sys', 'Pr(>|t|)']
#   cl.delta <- cliff.delta(model$mgs, model$BP.Sys)
#   print(cl.delta$estimate)
#   print(cl.delta$conf.int[1]) #lower
#   print(cl.delta$conf.int[2])
#   delta <- cl.delta$estimate
#   delta.low <- cl.delta$conf.int[1] #lower
#   delta.up <- cl.delta$conf.int[2] #upper
#   result <-
#     c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
#   
#   result.BP.Sys[[x]] <- result
# }
# BP.Sys.lm <- do.call(rbind.data.frame, result.BP.Sys)
# names(BP.Sys.lm) <-
#   c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
# BP.Sys.lm$fdr <- p.adjust(BP.Sys.lm$pval, method = 'fdr')
# BP.Sys.lm[BP.Sys.lm$fdr < .1, ]
# 
# # MGS vs BP.Dias
# result.BP.Dias <- list()
# for (x in seq(from = 25, to = 147)) {
#   print(names(progression.model[x]))
#   # result.model$MGS[x] <- names(mgs.model.bsl[x])
#   model <-
#     data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
#   mod.a <-
#     lm(mgs ~ BP.Dias + Gender + CenterID + Age,
#        data = model,
#        na.action = na.omit)
#   summ.mod.a <- summary(mod.a)
#   print(summ.mod.a$coefficients['BP.Dias', 'Pr(>|t|)'])
#   p.val <- summ.mod.a$coefficients['BP.Dias', 'Pr(>|t|)']
#   cl.delta <- cliff.delta(model$mgs, model$BP.Dias)
#   print(cl.delta$estimate)
#   print(cl.delta$conf.int[1]) #lower
#   print(cl.delta$conf.int[2])
#   delta <- cl.delta$estimate
#   delta.low <- cl.delta$conf.int[1] #lower
#   delta.up <- cl.delta$conf.int[2] #upper
#   result <-
#     c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
#   result.BP.Dias[[x]] <- result
# }
# BP.Dias.lm <- do.call(rbind.data.frame, result.BP.Dias)
# names(BP.Dias.lm) <-
#   c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
# BP.Dias.lm$fdr <- p.adjust(BP.Dias.lm$pval, method = 'fdr')
# BP.Dias.lm[BP.Dias.lm$fdr < .1, ]
# 
# # save results
# write.table(height.lm, 'GENUS_delta_lm_height.csv', sep = ';')
# write.table(weight.lm, 'GENUS_delta_lm_weight.csv', sep = ';')
# write.table(waist.lm, 'GENUS_delta_lm_waist.csv', sep = ';')
# write.table(bmi.lm, 'GENUS_delta_lm_BMI.csv', sep = ';')
# write.table(Fast.Glucose.lm, 'GENUS_delta_lm_fastGluc.csv', sep = ';')
# write.table(Basal.Glucose.lm, 'GENUS_delta_lm_basGluc.csv', sep = ';')
# write.table(Mean.Glucose.lm, 'GENUS_delta_lm_meanGluc.csv', sep = ';')
# write.table(Glucose.2h.lm, 'GENUS_delta_lm_2hGluc.csv', sep = ';')
# write.table(Fast.Insulin.lm, 'GENUS_delta_lm_fasIns.csv', sep = ';')
# write.table(Basal.Insulin.lm, 'GENUS_delta_lm_basIns.csv', sep = ';')
# write.table(Mean.Insulin.lm, 'GENUS_delta_lm_meanIns.csv', sep = ';')
# write.table(Basal.InsSecr.lm, 'GENUS_delta_lm_basInsSecr.csv', sep = ';')
# write.table(Total.InsSecr.lm, 'GENUS_delta_lm_totalInsSecr.csv', sep = ';')
# write.table(Gluc.Sensitivity.lm, 'GENUS_delta_lm_GlucSens.csv', sep = ';')
# write.table(OGIS.2h.lm, 'GENUS_delta_lm_OGIS2h.csv', sep = ';')
# write.table(Matsuda.lm, 'GENUS_delta_lm_Matsuda.csv', sep = ';')
# write.table(BP.Sys.lm, 'GENUS_delta_lm_BPSys.csv', sep = ';')
# write.table(BP.Dias.lm, 'GENUS_delta_lm_BPDias.csv', sep = ';')
# tax.1 <- tax[row.names(tax) %in% names(genus), ]
# write.table(tax.1, 'GENUS_delta_lm_TAXONOMY.csv', sep = ';')
# 
# # GMM modelling ----
# GMMs <- read.table('data/WP2.1_775_IGR_GMM-QMP_igrHA_Ent_M.csv', header=T, sep='\t', row.names=1)
# GMMs <- GMMs[,7:103]
# GMMs <- GMMs[row.names(GMMs)%in%delta.df$sample,]
# 
# progression.model <- data.frame(delta.df, GMMs)
# 
# # MGS vs Height
# result.height <- list()
# for (x in seq(from = 25, to = 121)) {
#   print(names(progression.model[x]))
#   model <-
#     data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
#   mod.a <-
#     lm(mgs ~ Height + Gender + CenterID + Age,
#        data = model,
#        na.action = na.omit)
#   summ.mod.a <- summary(mod.a)
#   print(summ.mod.a$coefficients['Height', 'Pr(>|t|)'])
#   p.val <- summ.mod.a$coefficients['Height', 'Pr(>|t|)']
#   cl.delta <- cliff.delta(model$mgs, model$Height)
#   print(cl.delta$estimate)
#   print(cl.delta$conf.int[1]) #lower
#   print(cl.delta$conf.int[2])
#   delta <- cl.delta$estimate
#   delta.low <- cl.delta$conf.int[1] #lower
#   delta.up <- cl.delta$conf.int[2] #upper
#   result <-
#     c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
#   result.height[[x]] <- result
# }
# 
# height.lm <- do.call(rbind.data.frame, result.height)
# names(height.lm) <-
#   c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
# height.lm$fdr <- p.adjust(height.lm$pval, method = 'fdr')
# height.lm[height.lm$fdr < .1, ]
# height.lm[height.lm$pval < .05, ]
# 
# # MGS vs Weight
# result.weight <- list()
# for (x in seq(from = 25, to = 121)) {
#   print(names(progression.model[x]))
#   model <-
#     data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
#   mod.a <-
#     lm(mgs ~ Weight + Gender + CenterID + Age,
#        data = model,
#        na.action = na.omit)
#   summ.mod.a <- summary(mod.a)
#   print(summ.mod.a$coefficients['Weight', 'Pr(>|t|)'])
#   p.val <- summ.mod.a$coefficients['Weight', 'Pr(>|t|)']
#   cl.delta <- cliff.delta(model$mgs, model$Weight)
#   print(cl.delta$estimate)
#   print(cl.delta$conf.int[1]) #lower
#   print(cl.delta$conf.int[2])
#   delta <- cl.delta$estimate
#   delta.low <- cl.delta$conf.int[1] #lower
#   delta.up <- cl.delta$conf.int[2] #upper
#   result <-
#     c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
#   result.weight[[x]] <- result
# }
# 
# weight.lm <- do.call(rbind.data.frame, result.weight)
# names(weight.lm) <-
#   c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
# weight.lm$fdr <- p.adjust(weight.lm$pval, method = 'fdr')
# weight.lm[weight.lm$fdr < .1, ]
# weight.lm[weight.lm$pval < .05, ]
# 
# # MGS vs Waist
# result.waist <- list()
# for (x in seq(from = 25, to = 121)) {
#   print(names(progression.model[x]))
#   # result.model$MGS[x] <- names(mgs.model.bsl[x])
#   model <-
#     data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
#   mod.a <-
#     lm(mgs ~ Waist + Gender + CenterID + Age,
#        data = model,
#        na.action = na.omit)
#   summ.mod.a <- summary(mod.a)
#   print(summ.mod.a$coefficients['Waist', 'Pr(>|t|)'])
#   p.val <- summ.mod.a$coefficients['Waist', 'Pr(>|t|)']
#   cl.delta <- cliff.delta(model$mgs, model$Waist)
#   print(cl.delta$estimate)
#   print(cl.delta$conf.int[1]) #lower
#   print(cl.delta$conf.int[2])
#   delta <- cl.delta$estimate
#   delta.low <- cl.delta$conf.int[1] #lower
#   delta.up <- cl.delta$conf.int[2] #upper
#   result <-
#     c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
#   
#   result.waist[[x]] <- result
# }
# 
# waist.lm <- do.call(rbind.data.frame, result.waist)
# names(waist.lm) <-
#   c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
# waist.lm$fdr <- p.adjust(waist.lm$pval, method = 'fdr')
# waist.lm[waist.lm$fdr < .1, ]
# waist.lm[waist.lm$pval<.05,]
# 
# # MGS vs BMI
# result.bmi <- list()
# for (x in seq(from = 25, to = 121)) {
#   print(names(progression.model[x]))
#   # result.model$MGS[x] <- names(mgs.model.bsl[x])
#   model <-
#     data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
#   mod.a <-
#     lm(mgs ~ BMI + Gender + CenterID + Age,
#        data = model,
#        na.action = na.omit)
#   summ.mod.a <- summary(mod.a)
#   print(summ.mod.a$coefficients['BMI', 'Pr(>|t|)'])
#   p.val <- summ.mod.a$coefficients['BMI', 'Pr(>|t|)']
#   cl.delta <- cliff.delta(model$mgs, model$BMI)
#   print(cl.delta$estimate)
#   print(cl.delta$conf.int[1]) #lower
#   print(cl.delta$conf.int[2])
#   delta <- cl.delta$estimate
#   delta.low <- cl.delta$conf.int[1] #lower
#   delta.up <- cl.delta$conf.int[2] #upper
#   result <-
#     c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
#   result.bmi[[x]] <- result
# }
# 
# bmi.lm <- do.call(rbind.data.frame, result.bmi)
# names(bmi.lm) <- c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
# bmi.lm$fdr <- p.adjust(bmi.lm$pval, method = 'fdr')
# bmi.lm[bmi.lm$fdr < .1, ]
# bmi.lm[bmi.lm$pval<.05,]
# 
# # MGS vs Fast.Glucose
# result.Fast.Glucose <- list()
# for (x in seq(from = 25, to = 121)) {
#   print(names(progression.model[x]))
#   # result.model$MGS[x] <- names(mgs.model.bsl[x])
#   model <-
#     data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
#   mod.a <-
#     lm(mgs ~ Fast.Glucose + Gender + CenterID + Age,
#        data = model,
#        na.action = na.omit)
#   summ.mod.a <- summary(mod.a)
#   print(summ.mod.a$coefficients['Fast.Glucose', 'Pr(>|t|)'])
#   p.val <- summ.mod.a$coefficients['Fast.Glucose', 'Pr(>|t|)']
#   cl.delta <- cliff.delta(model$mgs, model$Fast.Glucose)
#   print(cl.delta$estimate)
#   print(cl.delta$conf.int[1]) #lower
#   print(cl.delta$conf.int[2])
#   delta <- cl.delta$estimate
#   delta.low <- cl.delta$conf.int[1] #lower
#   
#   delta.up <- cl.delta$conf.int[2] #upper
#   result <-
#     c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
#   result.Fast.Glucose[[x]] <- result
# }
# Fast.Glucose.lm <- do.call(rbind.data.frame, result.Fast.Glucose)
# names(Fast.Glucose.lm) <-
#   c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
# Fast.Glucose.lm$fdr <- p.adjust(Fast.Glucose.lm$pval, method = 'fdr')
# Fast.Glucose.lm[Fast.Glucose.lm$fdr < .1, ]
# 
# # MGS vs Basal.Glucose
# result.Basal.Glucose <- list()
# for (x in seq(from = 25, to = 121)) {
#   print(names(progression.model[x]))
#   # result.model$MGS[x] <- names(mgs.model.bsl[x])
#   model <-
#     data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
#   mod.a <-
#     lm(mgs ~ Basal.Glucose + Gender + CenterID + Age,
#        data = model,
#        na.action = na.omit)
#   summ.mod.a <- summary(mod.a)
#   print(summ.mod.a$coefficients['Basal.Glucose', 'Pr(>|t|)'])
#   p.val <- summ.mod.a$coefficients['Basal.Glucose', 'Pr(>|t|)']
#   cl.delta <- cliff.delta(model$mgs, model$Basal.Glucose)
#   print(cl.delta$estimate)
#   print(cl.delta$conf.int[1]) #lower
#   print(cl.delta$conf.int[2])
#   delta <- cl.delta$estimate
#   delta.low <- cl.delta$conf.int[1] #lower
#   delta.up <- cl.delta$conf.int[2] #upper
#   result <-
#     c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
#   result.Basal.Glucose[[x]] <- result
# }
# Basal.Glucose.lm <- do.call(rbind.data.frame, result.Basal.Glucose)
# names(Basal.Glucose.lm) <-
#   c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
# Basal.Glucose.lm$fdr <-
#   p.adjust(Basal.Glucose.lm$pval, method = 'fdr')
# Basal.Glucose.lm[Basal.Glucose.lm$fdr < .1, ]
# 
# # MGS vs Mean.Glucose
# result.Mean.Glucose <- list()
# for (x in seq(from = 25, to = 121)) {
#   print(names(progression.model[x]))
#   # result.model$MGS[x] <- names(mgs.model.bsl[x])
#   model <-
#     data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
#   mod.a <-
#     lm(mgs ~ Mean.Glucose + Gender + CenterID + Age,
#        data = model,
#        na.action = na.omit)
#   summ.mod.a <- summary(mod.a)
#   print(summ.mod.a$coefficients['Mean.Glucose', 'Pr(>|t|)'])
#   p.val <- summ.mod.a$coefficients['Mean.Glucose', 'Pr(>|t|)']
#   cl.delta <- cliff.delta(model$mgs, model$Mean.Glucose)
#   print(cl.delta$estimate)
#   print(cl.delta$conf.int[1]) #lower
#   print(cl.delta$conf.int[2])
#   delta <- cl.delta$estimate
#   delta.low <- cl.delta$conf.int[1] #lower
#   delta.up <- cl.delta$conf.int[2] #upper
#   result <-
#     c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
#   result.Mean.Glucose[[x]] <- result
# }
# Mean.Glucose.lm <- do.call(rbind.data.frame, result.Mean.Glucose)
# names(Mean.Glucose.lm) <-
#   c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
# Mean.Glucose.lm$fdr <- p.adjust(Mean.Glucose.lm$pval, method = 'fdr')
# Mean.Glucose.lm[Mean.Glucose.lm$fdr < .1, ]
# 
# # MGS vs Glucose.2h
# result.Glucose.2h <- list()
# for (x in seq(from = 25, to = 121)) {
#   print(names(progression.model[x]))
#   # result.model$MGS[x] <- names(mgs.model.bsl[x])
#   model <-
#     data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
#   mod.a <-
#     lm(mgs ~ Glucose.2h + Gender + CenterID + Age,
#        data = model,
#        na.action = na.omit)
#   summ.mod.a <- summary(mod.a)
#   print(summ.mod.a$coefficients['Glucose.2h', 'Pr(>|t|)'])
#   p.val <- summ.mod.a$coefficients['Glucose.2h', 'Pr(>|t|)']
#   cl.delta <- cliff.delta(model$mgs, model$Glucose.2h)
#   print(cl.delta$estimate)
#   print(cl.delta$conf.int[1]) #lower
#   print(cl.delta$conf.int[2])
#   delta <- cl.delta$estimate
#   delta.low <- cl.delta$conf.int[1] #lower
#   delta.up <- cl.delta$conf.int[2] #upper
#   result <-
#     c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
#   result.Glucose.2h[[x]] <- result
# }
# Glucose.2h.lm <- do.call(rbind.data.frame, result.Glucose.2h)
# names(Glucose.2h.lm) <-
#   c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
# Glucose.2h.lm$fdr <- p.adjust(Glucose.2h.lm$pval, method = 'fdr')
# Glucose.2h.lm[Glucose.2h.lm$fdr < .1, ]
# 
# # MGS vs Fast.Insulin
# result.Fast.Insulin <- list()
# for (x in seq(from = 25, to = 121)) {
#   print(names(progression.model[x]))
#   # result.model$MGS[x] <- names(mgs.model.bsl[x])
#   model <-
#     data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
#   mod.a <-
#     lm(mgs ~ Fast.Insulin + Gender + CenterID + Age,
#        data = model,
#        na.action = na.omit)
#   summ.mod.a <- summary(mod.a)
#   print(summ.mod.a$coefficients['Fast.Insulin', 'Pr(>|t|)'])
#   p.val <- summ.mod.a$coefficients['Fast.Insulin', 'Pr(>|t|)']
#   cl.delta <- cliff.delta(model$mgs, model$Fast.Insulin)
#   print(cl.delta$estimate)
#   print(cl.delta$conf.int[1]) #lower
#   print(cl.delta$conf.int[2])
#   delta <- cl.delta$estimate
#   delta.low <- cl.delta$conf.int[1] #lower
#   
#   delta.up <- cl.delta$conf.int[2] #upper
#   result <-
#     c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
#   result.Fast.Insulin[[x]] <- result
# }
# Fast.Insulin.lm <- do.call(rbind.data.frame, result.Fast.Insulin)
# names(Fast.Insulin.lm) <-
#   c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
# Fast.Insulin.lm$fdr <- p.adjust(Fast.Insulin.lm$pval, method = 'fdr')
# Fast.Insulin.lm[Fast.Insulin.lm$fdr < .1, ]
# 
# # MGS vs Basal.Insulin
# result.Basal.Insulin <- list()
# for (x in seq(from = 25, to = 121)) {
#   print(names(progression.model[x]))
#   # result.model$MGS[x] <- names(mgs.model.bsl[x])
#   model <-
#     data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
#   mod.a <-
#     lm(mgs ~ Basal.Insulin + Gender + CenterID + Age,
#        data = model,
#        na.action = na.omit)
#   summ.mod.a <- summary(mod.a)
#   print(summ.mod.a$coefficients['Basal.Insulin', 'Pr(>|t|)'])
#   p.val <- summ.mod.a$coefficients['Basal.Insulin', 'Pr(>|t|)']
#   cl.delta <- cliff.delta(model$mgs, model$Basal.Insulin)
#   print(cl.delta$estimate)
#   print(cl.delta$conf.int[1]) #lower
#   print(cl.delta$conf.int[2])
#   delta <- cl.delta$estimate
#   delta.low <- cl.delta$conf.int[1] #lower
#   delta.up <- cl.delta$conf.int[2] #upper
#   result <-
#     c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
#   result.Basal.Insulin[[x]] <- result
# }
# Basal.Insulin.lm <- do.call(rbind.data.frame, result.Basal.Insulin)
# names(Basal.Insulin.lm) <-
#   c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
# Basal.Insulin.lm$fdr <-
#   p.adjust(Basal.Insulin.lm$pval, method = 'fdr')
# Basal.Insulin.lm[Basal.Insulin.lm$fdr < .1, ]
# 
# # MGS vs Mean.Insulin
# result.Mean.Insulin <- list()
# for (x in seq(from = 25, to = 121)) {
#   print(names(progression.model[x]))
#   # result.model$MGS[x] <- names(mgs.model.bsl[x])
#   model <-
#     data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
#   mod.a <-
#     lm(mgs ~ Mean.Insulin + Gender + CenterID + Age,
#        data = model,
#        na.action = na.omit)
#   summ.mod.a <- summary(mod.a)
#   print(summ.mod.a$coefficients['Mean.Insulin', 'Pr(>|t|)'])
#   p.val <- summ.mod.a$coefficients['Mean.Insulin', 'Pr(>|t|)']
#   cl.delta <- cliff.delta(model$mgs, model$Mean.Insulin)
#   print(cl.delta$estimate)
#   print(cl.delta$conf.int[1]) #lower
#   print(cl.delta$conf.int[2])
#   delta <- cl.delta$estimate
#   delta.low <- cl.delta$conf.int[1] #lower
#   delta.up <- cl.delta$conf.int[2] #upper
#   result <-
#     c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
#   
#   result.Mean.Insulin[[x]] <- result
# }
# Mean.Insulin.lm <- do.call(rbind.data.frame, result.Mean.Insulin)
# names(Mean.Insulin.lm) <-
#   c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
# Mean.Insulin.lm$fdr <- p.adjust(Mean.Insulin.lm$pval, method = 'fdr')
# Mean.Insulin.lm[Mean.Insulin.lm$fdr < .1, ]
# 
# # MGS vs Basal.InsSecr
# result.Basal.InsSecr <- list()
# for (x in seq(from = 25, to = 121)) {
#   print(names(progression.model[x]))
#   # result.model$MGS[x] <- names(mgs.model.bsl[x])
#   model <-
#     data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
#   mod.a <-
#     lm(mgs ~ Basal.InsSecr + Gender + CenterID + Age,
#        data = model,
#        na.action = na.omit)
#   summ.mod.a <- summary(mod.a)
#   print(summ.mod.a$coefficients['Basal.InsSecr', 'Pr(>|t|)'])
#   p.val <- summ.mod.a$coefficients['Basal.InsSecr', 'Pr(>|t|)']
#   cl.delta <- cliff.delta(model$mgs, model$Basal.InsSecr)
#   print(cl.delta$estimate)
#   print(cl.delta$conf.int[1]) #lower
#   print(cl.delta$conf.int[2])
#   delta <- cl.delta$estimate
#   delta.low <- cl.delta$conf.int[1] #lower
#   delta.up <- cl.delta$conf.int[2] #upper
#   result <-
#     c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
#   result.Basal.InsSecr[[x]] <- result
# }
# Basal.InsSecr.lm <- do.call(rbind.data.frame, result.Basal.InsSecr)
# names(Basal.InsSecr.lm) <-
#   c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
# Basal.InsSecr.lm$fdr <-
#   p.adjust(Basal.InsSecr.lm$pval, method = 'fdr')
# Basal.InsSecr.lm[Basal.InsSecr.lm$fdr < .1, ]
# 
# # MGS vs Total.InsSecr
# result.Total.InsSecr <- list()
# for (x in seq(from = 25, to = 121)) {
#   print(names(progression.model[x]))
#   # result.model$MGS[x] <- names(mgs.model.bsl[x])
#   model <-
#     data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
#   mod.a <-
#     lm(mgs ~ Total.InsSecr + Gender + CenterID + Age,
#        data = model,
#        na.action = na.omit)
#   summ.mod.a <- summary(mod.a)
#   print(summ.mod.a$coefficients['Total.InsSecr', 'Pr(>|t|)'])
#   p.val <- summ.mod.a$coefficients['Total.InsSecr', 'Pr(>|t|)']
#   cl.delta <- cliff.delta(model$mgs, model$Total.InsSecr)
#   print(cl.delta$estimate)
#   print(cl.delta$conf.int[1]) #lower
#   print(cl.delta$conf.int[2])
#   delta <- cl.delta$estimate
#   delta.low <- cl.delta$conf.int[1] #lower
#   
#   delta.up <- cl.delta$conf.int[2] #upper
#   result <-
#     c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
#   result.Total.InsSecr[[x]] <- result
# }
# Total.InsSecr.lm <- do.call(rbind.data.frame, result.Total.InsSecr)
# names(Total.InsSecr.lm) <-
#   c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
# Total.InsSecr.lm$fdr <-
#   p.adjust(Total.InsSecr.lm$pval, method = 'fdr')
# Total.InsSecr.lm[Total.InsSecr.lm$fdr < .1, ]
# 
# # MGS vs Gluc.Sensitivity
# result.Gluc.Sensitivity <- list()
# for (x in seq(from = 25, to = 121)) {
#   print(names(progression.model[x]))
#   # result.model$MGS[x] <- names(mgs.model.bsl[x])
#   model <-
#     data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
#   mod.a <-
#     lm(
#       mgs ~ Gluc.Sensitivity + Gender + CenterID + Age,
#       data = model,
#       na.action = na.omit
#     )
#   summ.mod.a <- summary(mod.a)
#   print(summ.mod.a$coefficients['Gluc.Sensitivity', 'Pr(>|t|)'])
#   p.val <- summ.mod.a$coefficients['Gluc.Sensitivity', 'Pr(>|t|)']
#   cl.delta <- cliff.delta(model$mgs, model$Gluc.Sensitivity)
#   print(cl.delta$estimate)
#   print(cl.delta$conf.int[1]) #lower
#   print(cl.delta$conf.int[2])
#   delta <- cl.delta$estimate
#   delta.low <- cl.delta$conf.int[1] #lower
#   delta.up <- cl.delta$conf.int[2] #upper
#   result <-
#     c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
#   result.Gluc.Sensitivity[[x]] <- result
# }
# Gluc.Sensitivity.lm <-
#   do.call(rbind.data.frame, result.Gluc.Sensitivity)
# names(Gluc.Sensitivity.lm) <-
#   c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
# Gluc.Sensitivity.lm$fdr <-
#   p.adjust(Gluc.Sensitivity.lm$pval, method = 'fdr')
# Gluc.Sensitivity.lm[Gluc.Sensitivity.lm$fdr < .1, ]
# 
# # MGS vs OGIS.2h
# result.OGIS.2h <- list()
# for (x in seq(from = 25, to = 121)) {
#   print(names(progression.model[x]))
#   # result.model$MGS[x] <- names(mgs.model.bsl[x])
#   model <-
#     data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
#   mod.a <-
#     lm(mgs ~ OGIS.2h + Gender + CenterID + Age,
#        data = model,
#        na.action = na.omit)
#   summ.mod.a <- summary(mod.a)
#   print(summ.mod.a$coefficients['OGIS.2h', 'Pr(>|t|)'])
#   p.val <- summ.mod.a$coefficients['OGIS.2h', 'Pr(>|t|)']
#   cl.delta <- cliff.delta(model$mgs, model$OGIS.2h)
#   print(cl.delta$estimate)
#   print(cl.delta$conf.int[1]) #lower
#   print(cl.delta$conf.int[2])
#   delta <- cl.delta$estimate
#   delta.low <- cl.delta$conf.int[1] #lower
#   delta.up <- cl.delta$conf.int[2] #upper
#   result <-
#     c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
#   
#   result.OGIS.2h[[x]] <- result
# }
# OGIS.2h.lm <- do.call(rbind.data.frame, result.OGIS.2h)
# names(OGIS.2h.lm) <-
#   c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
# OGIS.2h.lm$fdr <- p.adjust(OGIS.2h.lm$pval, method = 'fdr')
# OGIS.2h.lm[OGIS.2h.lm$fdr < .1, ]
# 
# # MGS vs Matsuda
# result.Matsuda <- list()
# for (x in seq(from = 25, to = 121)) {
#   print(names(progression.model[x]))
#   # result.model$MGS[x] <- names(mgs.model.bsl[x])
#   model <-
#     data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
#   mod.a <-
#     lm(mgs ~ Matsuda + Gender + CenterID + Age,
#        data = model,
#        na.action = na.omit)
#   summ.mod.a <- summary(mod.a)
#   print(summ.mod.a$coefficients['Matsuda', 'Pr(>|t|)'])
#   p.val <- summ.mod.a$coefficients['Matsuda', 'Pr(>|t|)']
#   cl.delta <- cliff.delta(model$mgs, model$Matsuda)
#   print(cl.delta$estimate)
#   print(cl.delta$conf.int[1]) #lower
#   print(cl.delta$conf.int[2])
#   delta <- cl.delta$estimate
#   delta.low <- cl.delta$conf.int[1] #lower
#   delta.up <- cl.delta$conf.int[2] #upper
#   result <-
#     c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
#   result.Matsuda[[x]] <- result
# }
# Matsuda.lm <- do.call(rbind.data.frame, result.Matsuda)
# names(Matsuda.lm) <-
#   c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
# Matsuda.lm$fdr <- p.adjust(Matsuda.lm$pval, method = 'fdr')
# Matsuda.lm[Matsuda.lm$fdr < .1, ]
# 
# # MGS vs BP.Sys
# result.BP.Sys <- list()
# for (x in seq(from = 25, to = 121)) {
#   print(names(progression.model[x]))
#   # result.model$MGS[x] <- names(mgs.model.bsl[x])
#   model <-
#     data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
#   mod.a <-
#     lm(mgs ~ BP.Sys + Gender + CenterID + Age,
#        data = model,
#        na.action = na.omit)
#   summ.mod.a <- summary(mod.a)
#   print(summ.mod.a$coefficients['BP.Sys', 'Pr(>|t|)'])
#   p.val <- summ.mod.a$coefficients['BP.Sys', 'Pr(>|t|)']
#   cl.delta <- cliff.delta(model$mgs, model$BP.Sys)
#   print(cl.delta$estimate)
#   print(cl.delta$conf.int[1]) #lower
#   print(cl.delta$conf.int[2])
#   delta <- cl.delta$estimate
#   delta.low <- cl.delta$conf.int[1] #lower
#   delta.up <- cl.delta$conf.int[2] #upper
#   result <-
#     c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
#   
#   result.BP.Sys[[x]] <- result
# }
# BP.Sys.lm <- do.call(rbind.data.frame, result.BP.Sys)
# names(BP.Sys.lm) <-
#   c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
# BP.Sys.lm$fdr <- p.adjust(BP.Sys.lm$pval, method = 'fdr')
# BP.Sys.lm[BP.Sys.lm$fdr < .1, ]
# 
# # MGS vs BP.Dias
# result.BP.Dias <- list()
# for (x in seq(from = 25, to = 121)) {
#   print(names(progression.model[x]))
#   # result.model$MGS[x] <- names(mgs.model.bsl[x])
#   model <-
#     data.frame('mgs' = progression.model[, x], progression.model[, 1:24])
#   mod.a <-
#     lm(mgs ~ BP.Dias + Gender + CenterID + Age,
#        data = model,
#        na.action = na.omit)
#   summ.mod.a <- summary(mod.a)
#   print(summ.mod.a$coefficients['BP.Dias', 'Pr(>|t|)'])
#   p.val <- summ.mod.a$coefficients['BP.Dias', 'Pr(>|t|)']
#   cl.delta <- cliff.delta(model$mgs, model$BP.Dias)
#   print(cl.delta$estimate)
#   print(cl.delta$conf.int[1]) #lower
#   print(cl.delta$conf.int[2])
#   delta <- cl.delta$estimate
#   delta.low <- cl.delta$conf.int[1] #lower
#   delta.up <- cl.delta$conf.int[2] #upper
#   result <-
#     c(names(progression.model)[x], p.val, delta, delta.low, delta.up)
#   result.BP.Dias[[x]] <- result
# }
# BP.Dias.lm <- do.call(rbind.data.frame, result.BP.Dias)
# names(BP.Dias.lm) <-
#   c('mgs', 'pval', 'delta', 'delta.low', 'delta.up')
# BP.Dias.lm$fdr <- p.adjust(BP.Dias.lm$pval, method = 'fdr')
# BP.Dias.lm[BP.Dias.lm$fdr < .1, ]
# 
# # save results
# write.table(height.lm, 'GMM_delta_lm_height.csv', sep = ';')
# write.table(weight.lm, 'GMM_delta_lm_weight.csv', sep = ';')
# write.table(waist.lm, 'GMM_delta_lm_waist.csv', sep = ';')
# write.table(bmi.lm, 'GMM_delta_lm_BMI.csv', sep = ';')
# write.table(Fast.Glucose.lm, 'GMM_delta_lm_fastGluc.csv', sep = ';')
# write.table(Basal.Glucose.lm, 'GMM_delta_lm_basGluc.csv', sep = ';')
# write.table(Mean.Glucose.lm, 'GMM_delta_lm_meanGluc.csv', sep = ';')
# write.table(Glucose.2h.lm, 'GMM_delta_lm_2hGluc.csv', sep = ';')
# write.table(Fast.Insulin.lm, 'GMM_delta_lm_fasIns.csv', sep = ';')
# write.table(Basal.Insulin.lm, 'GMM_delta_lm_basIns.csv', sep = ';')
# write.table(Mean.Insulin.lm, 'GMM_delta_lm_meanIns.csv', sep = ';')
# write.table(Basal.InsSecr.lm, 'GMM_delta_lm_basInsSecr.csv', sep = ';')
# write.table(Total.InsSecr.lm, 'GMM_delta_lm_totalInsSecr.csv', sep = ';')
# write.table(Gluc.Sensitivity.lm, 'GMM_delta_lm_GlucSens.csv', sep = ';')
# write.table(OGIS.2h.lm, 'GMM_delta_lm_OGIS2h.csv', sep = ';')
# write.table(Matsuda.lm, 'GMM_delta_lm_Matsuda.csv', sep = ';')
# write.table(BP.Sys.lm, 'GMM_delta_lm_BPSys.csv', sep = ';')
# write.table(BP.Dias.lm, 'GMM_delta_lm_BPDias.csv', sep = ';')