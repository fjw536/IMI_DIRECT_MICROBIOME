if (!requireNamespace("BiocManager")) {
  install.packages("BiocManager")
}

packages <- c("ggplot2", "ggthemes", "tidyverse", "survival", "magrittr", "rms", "scales", "stringr")

is.installed <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) {
    BiocManager::install(new.pkg, ask=F)
  }
  sapply(pkg, require, character.only = TRUE)
}
is.installed(packages)

wideScreen <- function(howWide=Sys.getenv("COLUMNS")) {
  options(width=as.integer(howWide))
}
wideScreen()

theme_set(theme_tufte(base_family = "sans", base_size = 18) + theme(panel.border = element_rect(colour = "black", fill = NA), axis.text = element_text(colour = "black", size = 18), plot.background = element_rect(colour = "white")))

load("/Pedersen_2022/2015_60_Salomaa_Jain_dataFR02_FU17_2020-11-15.RData")
#Subset to data which includes the fecal samples
FR02 <- FR02[!is.na(FR02$Barcode),]

richness_data <- read.csv("/Pedersen_2022/richness_means.txt", header=F, sep=" ", col.names=c("Barcode", "Richness"))
needed_features <- c("PREVAL_DIAB", "BL_USE_RX_J01", "GRAVID", "BL_AGE", "BMI", "MEN", "SYSTM", "KOL", "HDL", "CURR_SMOKE", "TRIG", "FR02_GLUK_NOLLA", "FR02_GLUK_120", "HBA1C", "INCIDENT_DIAB_T2", "DIAB_T2_AGEDIFF", "Richness")
main_data <- merge(FR02, richness_data, by = "Barcode")
main_data <- main_data[,colnames(main_data) %in% needed_features]
main_data <- main_data[-which(main_data$BL_USE_RX_J01 == 1 | main_data$GRAVID==2),]
main_data <- main_data[which(main_data$PREVAL_DIAB == 0),]
#Remove participants with diabetes indicator values over set guidelines: FR02_GLUK_NOLLA >= 7, FR02_GLUK_120 >= 11.1 & HBA1C >= 48 (ignore NA values)
main_data <- main_data[which(main_data$FR02_GLUK_NOLLA < 7 | is.na(main_data$FR02_GLUK_NOLLA)),]
main_data <- main_data[which(main_data$FR02_GLUK_120 < 11.1 | is.na(main_data$FR02_GLUK_120)),]
main_data <- main_data[which(main_data$HBA1C < 48 | is.na(main_data$HBA1C)),]
main_data$NON_HDL <- main_data$KOL - main_data$HDL  


richness_tertiles <- quantile(main_data$Richness, c(0:3/3))
main_data$richness_tert = with(main_data, cut(Richness, richness_tertiles, include.lowest = T, labels = c("LGR", "MGR", "HGR")))
main_data$richness_tert_num = as.numeric(as.factor(main_data$richness_tert))
main_data$richness_tert_num = ifelse(main_data$richness_tert_num == 1,1,0)

cox_wrapper <- function(data,
                        predictors,
                        covariates,
                        status,
                        time_to_event,
                        alpha_level,
                        normalize,
                        test_ph_assumption) {
  if(normalize) {  
    if(class(data[, predictors]) == "numeric") {
      x <- data[, predictors]
      data[, predictors] <- (x - mean(x, na.rm = T))/sd(x, na.rm = T) 
    } else {
      data[, predictors] <- apply(data[, predictors], 2, FUN = function(x) {(x - mean(x, na.rm = T))/sd(x, na.rm = T) })
    }
  }
  ## Formulas ***************************
  linear_formulas <- lapply(predictors, function(x) {
    formula_data <- deparse(substitute(data))
    formula <- paste0("Surv(", formula_data, "$", time_to_event, ", ",formula_data,"$", status, ") ~ ",paste(covariates, collapse = "+"), " + ",x)
    return(formula)
  }) %>% 
    set_names(predictors)
  ## Cox regression *********************
  print("Cox")
  linear_cox_fit <- lapply(linear_formulas, function(x) {
    coxph(as.formula(x), data=data, x=TRUE)
  })
  ## Check PH assumptions ****************
  if(test_ph_assumption) {
    print("PH assumptions")
    ph_assumption <-  lapply(predictors, function(m) {
      test <- cox.zph(linear_cox_fit[[m]])
      p_values <- test$table[, "p"]
      # significant cases
      x <- which(p_values < 1)
      if(length(x) == 0) {
        return(NULL)
      }
      df <- data.frame(feature = m, variable_not_ph = names(x), p_value = p_values[x])
    }) %>% 
      do.call(rbind, .) %>%
      mutate(p_adj = p.adjust(p_value, "BH")) %>% 
      filter(p_value < alpha_level)
  }
  ## Results *****************************
  print("Results")
  results <- lapply(predictors, function(x) {
    df <- summary(linear_cox_fit[[x]])$coefficients %>% as.data.frame()
    df <- df[nrow(df), ] %>% 
      select(coef, "se(coef)", "z", "Pr(>|z|)") %>% 
      set_colnames(c("coef", "se_coef", "test_stat_value", "p")) %>% 
      mutate(test_stat = "Wald")
    df <- df %>% 
      mutate(predictor = x)
  }) %>% 
    do.call(rbind, .) 
  # Multiple testing correction
  results <- results %>%
    mutate(P_adjusted = p.adjust(p, "BH")) %>% 
    ungroup() %>% 
    group_by(predictor) 
  # Results in neat form for presentation
  neat_results <- results %>% 
    # filter(p == min(p)) %>% 
    ungroup() %>% 
    mutate(HR = round(exp(coef),3)) %>% 
    mutate(HR_lower_95 = round(exp(coef - 1.96*se_coef), 3),
           HR_upper_95 = round(exp(coef + 1.96*se_coef), 3),
           P = round(p, 5),
           Coefficient = round(coef, 3),
           "Coefficient SE" = round(se_coef, 3)) %>% 
    mutate(HR = paste0(HR, " (95% CI, ", HR_lower_95, "-", HR_upper_95, ")")) %>% 
    select(Predictor = predictor, Coefficient, "Coefficient SE", HR, "p","P_adjusted", "test_stat_value", "test_stat") %>%
    mutate(HR = ifelse(is.na(Coefficient), NA, HR))  %>% 
    filter(P_adjusted < alpha_level) %>% 
    arrange(P_adjusted) %>% 
    set_colnames(c("Predictor", "Coefficient", "Coefficient SE", "HR","P-value" ,"P (adjusted)", "Test Statistic Value", "Test Statistic"))
  # Results in a form more convenient for further manipulations
  results <- results %>% 
    ungroup %>% 
    mutate(PH = exp(coef)) %>% 
    mutate(p_adj = P_adjusted) %>% 
    mutate(direction = ifelse(coef < 0, "negative", "positive"))
  if(nrow(neat_results) == 0) {
    return(list(results = results))
  }
  if(test_ph_assumption) {
    if(nrow(neat_results) == 0) {
      return(list(results = results, ph_assumption = ph_assumption))
    }
    return(list(neat_results = neat_results, 
                results = results,
                ph_assumption = ph_assumption))
  }
  return(list(neat_results = neat_results, results = results, model = linear_cox_fit))
}

#Cox regression for incident type 2 diabetes, accounting for age and sex

#set variables
alpha_level <- 0.05 #to filter
status <- "INCIDENT_DIAB_T2"
time_to_event <- "DIAB_T2_AGEDIFF"
predictors <- c("Richness")
covariates <- c("BL_AGE", "MEN")
splines <- TRUE
normalize <- TRUE
test_ph_assumption <- FALSE
#Cox regression with previously defined function
set.seed(42)
age_sex_cox <- cox_wrapper(data = main_data,
                        predictors = predictors,
                        covariates = covariates,
                        alpha_level = alpha_level,
                        status = status,
                        time_to_event = time_to_event,
                        normalize = normalize,
                        test_ph_assumption = test_ph_assumption)

neat_cox <- data.frame(age_sex_cox$neat_results)

#Plot hazard ratio
HRdf <- data.frame(Richness = neat_cox$Predictor[which(neat_cox$P.value < 0.05)],
                      HR = str_split(neat_cox$HR," ")[[1]][1],
                      HR1 = str_split(str_split(neat_cox$HR, " ")[[1]][4], "-")[[1]][1],
                      HR2 = gsub(")", "", str_split(str_split(neat_cox$HR, " ")[[1]][4], "-")[[1]][2]))

p <- ggplot(data = HRdf, aes(y = Richness, x = as.numeric(as.character(HR)))) + 
      geom_pointrange(aes(xmin=as.numeric(as.character(HR1)), xmax=as.numeric(as.character(HR2))), lwd = 1) +
	  scale_x_continuous(limits = c(0.5, 1.5)) +
      xlab("Hazard ratio") +
      geom_vline(xintercept=c(1.0), linetype="dotted") +
	  theme(axis.title.y = element_blank())
ggsave("/Pedersen_2022/HR_plot.svg", plot=p, units="in", width=6, height=2.5)

#Plot Kaplan-Meier curves
kp_predictors <- neat_cox$Predictor[which(neat_cox$P.value < 0.05)]
kp_covariates <- covariates
kp_time_to_event <- time_to_event
kp_status <- status
kp_data <- main_data[,which(colnames(main_data) %in% c(kp_status, kp_time_to_event, kp_predictors, kp_covariates, "richness_tert"))]
kp_time <- seq(0, max(kp_data$DIAB_T2_AGEDIFF), by = .01)
kp_list <- list(NULL)
for (time in 1:length(kp_time)) {
kp_table <- lapply(kp_predictors, function(x) {
    return_table <- data.frame(groupkm(kp_data[x], Surv(kp_data$DIAB_T2_AGEDIFF, kp_data$INCIDENT_DIAB_T2), g=3, u=kp_time[time], pl=FALSE))
	return_table$Predictor <- x
	return_table$tertile <- c(1:3)
    return(return_table)
  })
kp_table <- do.call(rbind, kp_table)
kp_table$time <- kp_time[time]
kp_list[[time]] <- kp_table
}

kp_list <- do.call(rbind, kp_list)

p <-  ggplot(data = kp_list, aes(y = KM, x = time, group = factor(tertile))) + 
	  geom_line(aes(color = factor(tertile)), linewidth = 2) +
	  scale_color_manual(values = c("#eeee00", "#bebebe", "#0000ff"), labels = c("LGR", "MGR", "HGR")) +
	  scale_y_continuous(breaks = pretty_breaks()) +
	  guides(color = guide_legend(override.aes = list(linewidth = 6)), fill = "none") +
      xlab("Time (years)") + ylab("Probability of diabetes-free survival") +
	  labs(color = "Gene\nRichness")
ggsave("/Pedersen_2022/KP_plot.svg", plot=p, units="cm", width=30, height=20)

#Plot distributions of the tertiles (for inlays in the KP-plot)
p <- ggplot(kp_data, aes(x = Richness, group = richness_tert, fill = richness_tert)) + 
	  geom_density(alpha = 0.5) + scale_x_continuous(limits = c(min(kp_data$Richness)*0.95, max(kp_data$Richness)*1.05), labels = label_scientific()) +
	  scale_y_continuous(labels = label_scientific()) +
	  scale_fill_manual(values = c("#eeee00", "#bebebe", "#0000ff")) +
	  guides(fill = "none") +
	  xlab("Gene Richness") + ylab("Density")
ggsave("/Pedersen_2022/KP_plot_tertiles.svg", plot=p, units="cm", width=15, height=10)

#calculate P-values for KP for the groups
KP_logrank_test <- survdiff(Surv(kp_data$DIAB_T2_AGEDIFF, kp_data$INCIDENT_DIAB_T2) ~ richness_tert, data = kp_data)