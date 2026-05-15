# Figure 4a
# Variance of metabolites explained by various factors

metab_0m_775 <- readRDS('Inputs/metab_0m_775.rdata')

met_panel <- metab_0m_775[, sig_traj$feature]

# Ensure samples match metadata
common_ids <- Reduce(intersect, list(
  rownames(met_panel),
  rownames(bac_0m_775),
  rownames(meta_0m_775)
))

met_panel <- met_panel[common_ids, ]
bac_0m_775 <- bac_0m_775[common_ids, ]
meta_0m_775 <- meta_0m_775[common_ids, ]


prev <- colMeans(bac_0m_775 > 0, na.rm = TRUE)
bac_0m_775_filt <- bac_0m_775[, prev >= 0.10]

# kept_species <- colnames(bac_0m_775_filt)


library(vegan)

bac_dist <- vegdist(bac_0m_775_filt, method = "bray")

bac_pcoa <- cmdscale(bac_dist, k = 5, eig = TRUE)

bac_pcs <- as.data.frame(bac_pcoa$points)
colnames(bac_pcs) <- paste0("Microbiome_PC", 1:ncol(bac_pcs))


meta_cov <- meta_0m_775 %>%
  select(
    Age,
    Sex,
    `Body mass index (kg/m2)`,
    Smoking.Status.BL,
    Value.vm.hpf.mean.BL,
    HDI
  ) %>%
  mutate(
    SEX = factor(Sex),
    BMI = as.numeric(`Body mass index (kg/m2)`),
    smoking = factor(Smoking.Status.BL),
    exercise = as.numeric(Value.vm.hpf.mean.BL),
    diet = as.numeric(HDI)
  )


# Multivariate metabolite distance
met_dist <- vegdist(met_panel, method = "euclidean")

library(vegan)
# Step 1 - Identify those significantly contributing intestinal taxa
res <- lapply(colnames(bac_0m_775_filt), function(taxon) {

  f <- as.formula(
    paste0("met_dist ~ `", taxon, "`")
  )
  set.seed(10)
  ad <- adonis2(
    f,
    data = bac_0m_775_filt,
    permutations = 999,
    by = "margin"
  )

  data.frame(
    taxon  = taxon,
    R2     = ad$R2[1],
    F      = ad$F[1],
    pvalue = ad$`Pr(>F)`[1]
  )
})
res <- do.call(rbind, res)
res$qvalue <- p.adjust(res$pvalue, method = "BH")
sig_taxa <- subset(res, pvalue < 0.05)


# Step 2
taxa_step2 <- sig_taxa$taxon
rhs <- paste0("`", taxa_step2, "`", collapse = " + ")

form_step2 <- as.formula(
  paste("met_dist ~", rhs)
)

set.seed(10)
div_species_step2 <- adonis2(
  form_step2,
  data = bac_0m_775_filt,
  distance = 'bray',
  permutations = 999,
  by = 'margin'
)

div_species_step2_sig <- div_species_step2 %>% filter(`Pr(>F)` < 0.05)

rhs_step3 <- paste0(rownames(div_species_step2_sig), collapse = " + ")

form_step3 <- as.formula(
  paste("met_dist ~", rhs_step3)
)

div_species_step3 <- adonis2(
  form_step3,
  data = bac_0m_775_filt,
  distance = 'bray',
  permutations = 999,
  by = 'margin'
)

div_species_final <- adonis2(
  form_step3,
  data = bac_0m_775_filt,
  distance = 'bray',
  permutations = 999
)
div_species_final <- div_species_step3
saveRDS(div_species_step3, file = 'Outputs/Microbiota_expalined_variance_of_met_traj.rdata')

# Intrinsic host factors
set.seed(10)
div_demo <- adonis2(met_dist ~ Age + Sex + BMI, data = meta_cov,  permutations = 1000, by= NULL)

# Lifestyle factors
set.seed(10)
div_life <- adonis2(met_dist ~ smoking + exercise, data = meta_cov,  permutations = 999, by= NULL)

# Diets
diet <- readRDS('Inputs/Diet_record_0m_739.rdata')
set.seed(10)
diet_res_step1 <- lapply(colnames(diet), function(taxon) {

  f <- as.formula(
    paste0("met_panel[rownames(diet), ] ~ ", taxon)
  )
  # set.seeds(10)
  ad <- adonis2(
    f,
    data = diet,
    permutations = 999,
    distance = 'bray',
    by = "margin"
  )

  data.frame(
    taxon  = taxon,
    R2     = ad$R2[1],
    F      = ad$F[1],
    pvalue = ad$`Pr(>F)`[1]
  )
})

diet_res_step1 <- do.call(rbind, diet_res_step1)
diet_res_step1$qvalue <- p.adjust(diet_res_step1$pvalue, method = "BH")
sig_diet <- subset(diet_res_step1, pvalue < 0.05)

# Step 2
diet_step2 <- sig_diet$taxon
rhs <- paste0(diet_step2, collapse = " + ")

form_step2 <- as.formula(
  paste("met_panel[rownames(diet), ] ~", rhs)
)

set.seed(10)
div_diet_step2 <- adonis2(
  form_step2,
  data = diet,
  distance = 'bray',
  permutations = 999,
  by = 'margin'
) # Takes long time to finish

div_diet_step2_sig <- div_diet_step2 %>% filter(`Pr(>F)` < 0.05)

rhs_step3 <- paste0(rownames(div_diet_step2_sig), collapse = " + ")

form_step3 <- as.formula(
  paste("met_panel[rownames(diet), ] ~", rhs_step3)
)

div_diet_step3 <- adonis2(
  form_step3,
  data = diet,
  distance = 'bray',
  permutations = 999,
  by = 'margin'
)

div_diet_step3_sig <- div_diet_step3 %>% filter(`Pr(>F)` < 0.05)

rhs_step4 <- paste0(rownames(div_diet_step3_sig), collapse = " + ")

form_step4 <- as.formula(
  paste("met_panel[rownames(diet), ] ~", rhs_step4)
)

div_diet_step4 <- adonis2(
  form_step4,
  data = diet,
  distance = 'bray',
  permutations = 999,
  by = 'margin'
) 

div_diet_final <- adonis2(
  form_step4,
  data = diet,
  distance = 'bray',
  permutations = 999
)

Explained_variance <- list(intrinsic_factor = div_demo,
                           lifestyle = div_life,
                           diet = div_diet_final,
                           gut_microbiota = div_species_final)

saveRDS(Explained_variance, file = 'All_factors_explained_variance_in_met_traj.rdata')

# Visualize the explained variance by host and microbiota factors

extract_model_adonis <- function(ad, label) {
  data.frame(
    Factor = label,
    Df     = ad$Df[1],
    R2     = ad$R2[1],
    F      = ad$F[1],
    Pvalue = ad$`Pr(>F)`[1],
    row.names = NULL
  )
}

library(dplyr)

explained_df <- bind_rows(
  extract_model_adonis(Explained_variance$intrinsic_factor, "Age, sex, and BMI"),
  extract_model_adonis(Explained_variance$lifestyle,        "Smoking and physical activity"),
  extract_model_adonis(Explained_variance$diet,             "Diet"),
  extract_model_adonis(Explained_variance$gut_microbiota,   "Gut microbiota")
)

explained_df <- explained_df %>%
  mutate(
    Factor = factor(
      Factor,
      levels = c(
        "Gut microbiota",
        "Age, sex, and BMI",
        "Diet",
        "Smoking and physical activity"
      )
    )
  )

library(ggplot2)

library(ggsci)

p_r2 <- ggplot(explained_df, aes(x = Factor, y = R2, fill = Factor)) +
  geom_col(
    width = 0.65
  ) +
  geom_text(
    aes(label = sprintf("%.1f%%", R2 * 100)),
    vjust = 2,
    size = 4,
    color = 'white'
  ) +
  scale_fill_manual(values = c("Gut microbiota" = '#2a6ebb',
                    "Diet" = '#f0ab00',
                    "Age, sex, and BMI" = '#c50084',
                    "Smoking and physical activity" = '#7d5cc6')) +   # 🔑 Science-style colors
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.15))
  ) +
  labs(
    x = NULL,
    y = "Explained variance (R²)"
  ) +
  # coord_flip() +
  theme_article(base_size = 13) +
  theme(
  axis.text.x = element_text(
    angle = 60,
    # vjust = 0.5,
    hjust = 1,
    color = "black"
  ),
  axis.text.y = element_text(color = "black"),
  legend.position = "none"
)

p_r2

pdf('Outputs/plot_overall_explained_var_met_traj.pdf', width = 4, height = 5)
print(p_r2)
dev.off()


# Figure 4b
##########################################################################################
##  Explained variance in each metabolite by different host and gut microbiota factors  ##
##########################################################################################
# Step 1️⃣ Define sample intersections explicitly
# sample IDs
ids_metab <- rownames(metab_0m_775)
ids_meta  <- rownames(meta_0m_775)
ids_bac   <- rownames(bac_0m_775_filt)
ids_diet  <- rownames(diet)   # 739 × 50 matrix

# Step 2️⃣ Prepare aligned datasets
# Intrinsic + lifestyle (non-diet)

ids_intrinsic <- Reduce(intersect, list(ids_metab, ids_meta))

Y_intrinsic <- metab_0m_775[ids_intrinsic, sig_traj$feature]

meta_intrinsic <- meta_0m_775[ids_intrinsic, ] %>%
  dplyr::select(Age, Sex, `Body mass index (kg/m2)`, Smoking.Status.BL, Value.vm.hpf.mean.BL)

meta_intrinsic$Sex <- factor(meta_intrinsic$Sex)

# Diet (separate!)
ids_diet_use <- Reduce(intersect, list(ids_metab, ids_diet))

Y_diet <- metab_0m_775[ids_diet_use, sig_traj$feature]
X_diet <- scale(diet[ids_diet_use, ])

# Gut microbiota
ids_bac_use <- Reduce(intersect, list(ids_metab, ids_bac))

Y_bac <- metab_0m_775[ids_bac_use, sig_traj$feature]
X_bac <- scale(bac_0m_775_filt[ids_bac_use, ])

# Step 3️⃣ LASSO + adjusted R² (factor-specific)
library(glmnet)

fit_lasso_adjR2 <- function(y, X) {

  cvfit <- cv.glmnet(
    x = X,
    y = y,
    alpha = 1,
    nfolds = 10,
    standardize = FALSE
  )

  beta <- coef(cvfit, s = "lambda.min")
  selected <- rownames(beta)[beta[,1] != 0]
  selected <- setdiff(selected, "(Intercept)")

  if (length(selected) == 0) return(0)

  df <- data.frame(
    y = y,
    X[, selected, drop = FALSE]
  )

  summary(lm(y ~ ., data = df))$adj.r.squared
}

# Step 4️⃣ Per-metabolite variance explained by each factor group
library(purrr)
library(dplyr)

var_explained_metab <- map_dfr(
  sig_traj$feature,
  function(met) {

    tibble(
      metabolite = met,

      intrinsic = fit_lasso_adjR2(
        Y_intrinsic[, met],
        model.matrix(~ Age + Sex + `Body mass index (kg/m2)`, meta_intrinsic)[, -1]
      ),

      lifestyle = fit_lasso_adjR2(
        Y_intrinsic[, met],
        model.matrix(~ Smoking.Status.BL + Value.vm.hpf.mean.BL, meta_intrinsic)[, -1]
      ),

      diet = fit_lasso_adjR2(
        Y_diet[, met],
        X_diet
      ),

      gut_microbiota = fit_lasso_adjR2(
        Y_bac[, met],
        X_bac
      )
    )
  }
)

# Step 5️⃣ Long format for visualization
var_explained_long <- var_explained_metab %>%
  tidyr::pivot_longer(
    cols = -metabolite,
    names_to = "factor",
    values_to = "adj_R2"
  )

var_explained_long$Chem <- metabolome_annotation[var_explained_long$metabolite, 'SHORT_NAME']

# Step 6️⃣ Visualization (clean & interpretable)

library(dplyr)
library(forcats)

plot_var_stack <- var_explained_long %>%
  mutate(
    adj_R2 = pmax(adj_R2, 0),  # avoid negative bars
    factor = factor(
      factor,
      levels = rev(c(
        "intrinsic",
        "lifestyle",
        "diet",
        "gut_microbiota"
      )
    ))
  ) %>%
  group_by(Chem) %>%
  mutate(total_R2 = sum(adj_R2)) %>%
  ungroup() %>%
  mutate(
    Chem = fct_reorder(Chem, total_R2)
  )

library(ggplot2)
library(ggsci)
# desired metabolite order from sig_traj
chem_order <- sig_traj %>%
  arrange(estimate) %>%            # or desc(estimate) if needed
  pull(feature_label) %>%
  unique()

plot_var_stack <- plot_var_stack %>%
  mutate(
    Chem = factor(Chem, levels = chem_order)
  )


p_stack <- ggplot(
  plot_var_stack,
  aes(x = Chem, y = adj_R2, fill = factor)
) +
  geom_col(width = 0.75) +
  coord_flip() +
  scale_fill_bmj(
    name = "Factors",
    labels = c(
      intrinsic      = "Age, sex, and BMI",
      lifestyle      = "Smonking and physical activity",
      diet           = "Diet",
      gut_microbiota = "Gut microbiota"
    )
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = NULL,
    y = "Explained variance (R²)"
  ) +
  theme_article(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 9),
    legend.position = "right",
    axis.text = element_text(color = 'black'))

p_stack

metabolite_traj_explained_var_left <- plot_summary_metab_traj_gly_48m
metabolite_traj_explained_var_right <- p_stack + theme(axis.text.y = element_blank(), legend.position = 'none')


metabolite_traj_explained_var <- egg::ggarrange(metabolite_traj_explained_var_left, metabolite_traj_explained_var_right, nrow = 1)
grid.arrange(
  p_r2,
  metabolite_traj_explained_var,
  nrow = 2,
  heights = c(1, 3)   # optional: control relative height
)
pdf('Outputs/Metab_traj_summary_and_explained_variance.pdf', width = 9, height = 7)
print(metabolite_traj_explained_var)
dev.off()