# Figure 5a
# Load libraries
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(broom)
library(multcomp)

# Data preparation for differential analysis
library(data.table)
library(stringr)
library(dplyr)
library(tidyr)

meta_fl <- readRDS('../Inputs/Followup_metadata.Rdata')
meta_fl <- meta_fl %>% filter(Gly.Cat.48m != 'nk')
rownames(meta_fl) <- meta_fl$rownames
long_dt <- fread("/Users/fjw536/Downloads/DIRECT_baseline.combined_profile(1).tsv", 
                 sep = "\t", header = T) %>%
  mutate(SGB = gsub(".*s__","s__", lineage ), 
         genus = gsub(".*g__","g__", lineage ), 
         genus = gsub("\\|s__.*", "", genus),
         kingdom = gsub("\\|.*", "", lineage))

wide_dt <- long_dt %>%
  subset( select = c("SGB", "rel_abundance", "sampleID")) %>%
  pivot_wider(names_from = sampleID, 
              values_from = rel_abundance, values_fill = 0)

wide_dt <- column_to_rownames(wide_dt, 'SGB')
wide_dt_clr <- clr(wide_dt) %>% as.data.frame()

# bac_dt: rows containing "t__SGB"
bac_dt <- wide_dt_clr %>%
  filter(str_detect(rownames(wide_dt_clr), "t__SGB"))

# vir_dt: all remaining rows
vir_dt <- wide_dt_clr %>%
  filter(!str_detect(rownames(wide_dt_clr), "t__SGB"))

bac_dt <- as.data.frame(bac_dt)
rownames(bac_dt) <- bac_dt$SGB
bac_dt <- bac_dt[, -1]
bac_dt <- t(bac_dt) %>% data.frame()

vir_dt <- as.data.frame(vir_dt)
rownames(vir_dt) <- vir_dt$SGB
vir_dt <- vir_dt[, -1]
vir_dt <- t(vir_dt) %>% data.frame()


bac_counts <- rowSums(bac_dt != 0)
vir_counts <- rowSums(vir_dt != 0)

# Combine into one dataframe
count_df <- data.frame(
  sample_id = rownames(bac_dt),
  bac_nonzero = bac_counts,
  vir_nonzero = vir_counts,
  row.names = NULL
)

library(ggplot2)

ggplot(count_df, aes(x = bac_nonzero, y = vir_nonzero)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  theme_minimal(base_size = 14) +
  labs(
    x = "Number of detected bacterial species",
    y = "Number of detected viral species"
  )



library(compositions)  # for clr
library(dplyr)

# --- Step 1: CLR transform ---
# Add small pseudocount to avoid log(0)
bac_clr <- clr(bac_dt)
vir_clr <- clr(vir_dt)

# --- Step 2: Spearman correlations ---
# Convert to matrices for efficiency
# Convert to matrices
bac_mat <- as.matrix(t(bac_clr))
vir_mat <- as.matrix(t(vir_clr))

# Filter by prevalence >= 20%
bac_prevalence <- colMeans(bac_mat != 0)
vir_prevalence <- colMeans(vir_mat != 0)

bac_keep <- names(bac_prevalence[bac_prevalence >= 0.1])
vir_keep <- names(vir_prevalence[vir_prevalence >= 0.1])

bac_mat_filt <- bac_mat[, bac_keep, drop = FALSE] %>% as.data.frame()
vir_mat_filt <- vir_mat[, vir_keep, drop = FALSE] %>% as.data.frame()
# Example input:
# df: a dataframe with columns
#   sampleID, feature1, feature2, ..., featureN, Group, Age, Sex
# Group should be factor with levels NGR, IGR, T2D

# Ensure Group is a factor
dif_df <- cbind(bac_mat_filt[rownames(meta_fl), ], meta_fl$Sex, meta_fl$bsl_age, meta_fl$bsl_centerID, meta_fl$Gly.Cat.48m)
colnames(dif_df)[457:460] <- c('Sex', 'Age', 'bsl_centerID', 'Gly.Cat.48m')

dif_df$Gly.Cat.48m <- factor(dif_df$Gly.Cat.48m, levels = c("NGR", "IGR", "T2D"))

# Function to run ANOVA + Tukey HSD for one feature
run_anova_tukey <- function(feature, data) {
  formula <- as.formula(paste(feature, "~ Gly.Cat.48m + Age + Sex + bsl_centerID"))
  fit <- aov(formula, data = data)
  
  # ANOVA table
  anova_res <- tidy(anova(fit))
  
  # Tukey HSD only if Group is significant
  tukey_res <- TukeyHSD(fit, "Gly.Cat.48m") %>%
    broom::tidy() %>%
    mutate(feature = feature)
  
  return(list(anova = anova_res, tukey = tukey_res))
}

# Apply across all microbiome features
features <- colnames(dif_df)[!(colnames(dif_df) %in% c("Gly.Cat.48m","Age","Sex", 'bsl_centerID'))]

results <- lapply(features, run_anova_tukey, data = dif_df)

# Combine into dataframes
anova_results <- bind_rows(lapply(results, function(x) x$anova), .id = "feature_id")
tukey_results <- bind_rows(lapply(results, function(x) x$tukey))

# Filter Tukey results for significant pairwise comparisons
tukey_sig <- tukey_results %>% filter(adj.p.value < 0.1)

# View results
head(anova_results)
head(tukey_sig)





library(dplyr)
library(tidyr)
library(pheatmap)
library(tibble)

# -------------------------------
# Step 1. Compute group means
# -------------------------------
group_means <- dif_df %>%
  dplyr::select(-Age, -Sex, -bsl_centerID, -Gly.Cat.48m) %>%  # drop metadata + group var
  bind_cols(Gly.Cat.48m = dif_df$Gly.Cat.48m) %>%             # add group var back
  group_by(Gly.Cat.48m) %>%
  summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop") %>%
  pivot_longer(-Gly.Cat.48m, names_to = "feature", values_to = "mean") %>%
  pivot_wider(names_from = Gly.Cat.48m, values_from = mean)

# -------------------------------
# Step 1b. Keep only features in tukey_sig
# -------------------------------
sig_features <- unique(tukey_sig$feature)

group_means <- group_means %>%
  filter(feature %in% sig_features)

# Convert to matrix
mat <- as.matrix(group_means[, -1])
rownames(mat) <- group_means$feature

# Scale to z-scores by row
mat_scaled <- t(scale(t(mat)))

# -------------------------------
# Step 2. Create annotation table
# -------------------------------
sig_contrasts <- tukey_sig %>%
  filter(adj.p.value < 0.1) %>%
  dplyr::select(feature, contrast)

anno <- sig_contrasts %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = contrast, values_from = value, values_fill = 0) %>%
  right_join(data.frame(feature = rownames(mat_scaled)), by = "feature") %>%
  replace(is.na(.), 0)

row_anno <- anno %>% column_to_rownames("feature")

# -------------------------------
# Step 4. Plot heatmap
# -------------------------------



library(ComplexHeatmap)
library(circlize)

# Define heatmap colors
col_fun <- colorRamp2(c(-1.5, 0, 1.5), c("#2e756e", "white", "#945e25"))

# Convert row annotation to data frame
row_anno_df <- as.data.frame(row_anno)

# Define annotation colors
ann_colors <- list(
  `IGR-NGR` = c("0" = "white", "1" = "#3cb4e9"),
  `T2D-NGR` = c("0" = "white", "1" = "#f486ad"),
  `T2D-IGR` = c("0" = "white", "1" = "#ecce88")
)

# Row annotation object with fixed width
ha <- rowAnnotation(
  df = row_anno_df,
  col = ann_colors,
  border = T,
  annotation_name_gp = gpar(fontsize = 10)    # make this equal to one heatmap column width
)


# Remove "s__" at the start
rownames(mat_scaled) <- gsub("^s__", "", rownames(mat_scaled))

# Remove ".t__XXXX" at the end
rownames(mat_scaled) <- gsub("\\.t__.*$", "", rownames(mat_scaled))


# 3. Replace "_" with space
rownames(mat_scaled) <- gsub("_", " ", rownames(mat_scaled))


# ha_col_feat = ComplexHeatmap::HeatmapAnnotation(Features=annot_col_feat$Features, annotation_name_side = "left", show_annotation_name = F, show_legend=c(T))
# Build heatmap
set.seed(100)
ht <- Heatmap(
  mat_scaled,
  name = "zscore",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = T,
  show_column_names = TRUE,
  right_annotation = ha,
  border = T,
  row_split = 2,
  row_title = NULL,
  column_names_gp = gpar(fontsize = 10),
  row_names_gp = gpar(fontsize = 10, fontface = "italic"),
  width = unit(ncol(mat_scaled), "cm"))
ht


pdf('Contrasted_bac_baseline.pdf', height = 10, width = 15)
print(ht)
dev.off()








# Ensure Group is a factor
vir_dif_df <- cbind(vir_mat_filt[rownames(meta_fl), ], meta_fl$Sex, meta_fl$bsl_age, meta_fl$bsl_centerID, meta_fl$Gly.Cat.48m)
colnames(vir_dif_df)[598:601] <- c('Sex', 'Age', 'bsl_centerID', 'Gly.Cat.48m')

vir_dif_df$Gly.Cat.48m <- factor(vir_dif_df$Gly.Cat.48m, levels = c("NGR", "IGR", "T2D"))

# Function to run ANOVA + Tukey HSD for one feature
run_anova_tukey <- function(feature, data) {
  formula <- as.formula(paste(feature, "~ Gly.Cat.48m + Age + Sex + bsl_centerID"))
  fit <- aov(formula, data = data)
  
  # ANOVA table
  anova_res <- tidy(anova(fit))
  
  # Tukey HSD only if Group is significant
  tukey_res <- TukeyHSD(fit, "Gly.Cat.48m") %>%
    broom::tidy() %>%
    mutate(feature = feature)
  
  return(list(anova = anova_res, tukey = tukey_res))
}

# Apply across all microbiome features
features <- colnames(vir_dif_df)[!(colnames(vir_dif_df) %in% c("Gly.Cat.48m","Age","Sex", 'bsl_centerID'))]

results <- lapply(features, run_anova_tukey, data = vir_dif_df)

# Combine into dataframes
vir_anova_results <- bind_rows(lapply(results, function(x) x$anova), .id = "feature_id")
vir_tukey_results <- bind_rows(lapply(results, function(x) x$tukey))

# Filter Tukey results for significant pairwise comparisons
vir_tukey_sig <- vir_tukey_results %>% filter(adj.p.value < 0.1)

# View results
head(anova_results)
head(tukey_sig)





library(dplyr)
library(tidyr)
library(pheatmap)
library(tibble)

# -------------------------------
# Step 1. Compute group means
# -------------------------------
group_means <- vir_dif_df %>%
  dplyr::select(-Age, -Sex, -bsl_centerID, -Gly.Cat.48m) %>%  # drop metadata + group var
  bind_cols(Gly.Cat.48m = dif_df$Gly.Cat.48m) %>%             # add group var back
  group_by(Gly.Cat.48m) %>%
  summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop") %>%
  pivot_longer(-Gly.Cat.48m, names_to = "feature", values_to = "mean") %>%
  pivot_wider(names_from = Gly.Cat.48m, values_from = mean)

# -------------------------------
# Step 1b. Keep only features in tukey_sig
# -------------------------------
sig_features <- unique(vir_tukey_sig$feature)

group_means <- group_means %>%
  filter(feature %in% sig_features)

# Convert to matrix
mat <- as.matrix(group_means[, -1])
rownames(mat) <- group_means$feature

# Scale to z-scores by row
mat_scaled <- t(scale(t(mat)))

# -------------------------------
# Step 2. Create annotation table
# -------------------------------
sig_contrasts <- vir_tukey_sig %>%
  filter(adj.p.value < 0.1) %>%
  dplyr::select(feature, contrast)

anno <- sig_contrasts %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = contrast, values_from = value, values_fill = 0) %>%
  right_join(data.frame(feature = rownames(mat_scaled)), by = "feature") %>%
  replace(is.na(.), 0)

row_anno <- anno %>% column_to_rownames("feature")

library(bannerCommenter)
##################################################################
##           Add host genus for contrasted viral SGBs           ##
##################################################################
vir_meta <- read.delim('/Users/fjw536/Downloads/v1.1/Marker-MAGu_virus_DB_v1.1_metadata.tsv', header = T)

library(dplyr)
library(stringr)

vir_meta <- vir_meta %>%
  # extract the vSGB_xxx from lineage
  mutate(SGB_id = str_extract(lineage, "vSGB_[0-9]+"),
         
         # extract the genus name after g__ in iphop_host_genus
         genus_taxa = str_extract(iphop_host_genus, "g__[^;]+"),
         genus_taxa = str_remove(genus_taxa, "^g__"))

rownames(vir_meta) <- vir_meta$SGB_id # Rename rows

rownames(row_anno) <- gsub("^s__", "", rownames(row_anno)) # Remove 's__' from the row names of data frame row_anno

row_anno$'Host genus' <- vir_meta[rownames(row_anno), 'genus_taxa']

row_anno$`Host genus`[is.na(row_anno$`Host genus`)] <- "Unknown" # Replace 'NA' by 'Unknown'

row_anno$'Virulence' <- vir_meta[rownames(row_anno), 'Virulence_score'] # Add virulence score for each vSGB

# Define virulence
row_anno <- row_anno %>%
  mutate(Virulence = case_when(
    is.na(Virulence) ~ "Unknown",
    Virulence > 0.5  ~ "Virulent",
    Virulence <= 0.5 ~ "Non-virulent"
  ))


# Define heatmap colors
col_fun <- colorRamp2(c(-1.5, 0, 1.5), c("#2e756e", "white", "#945e25"))

# Convert row annotation to data frame
row_anno_df <- as.data.frame(row_anno)

# Define annotation colors
ann_colors <- list(
  `IGR-NGR` = c("0" = "white", "1" = "#3cb4e9"),
  `T2D-NGR` = c("0" = "white", "1" = "#f486ad"),
  `T2D-IGR` = c("0" = "white", "1" = "#ecce88"),
  `Virulence` = c("Virulent" = "red", "Non-virulent" = "blue", "Unknown" = "grey")
)

# Row annotation object with fixed width
ha <- rowAnnotation(
  df = row_anno_df,
  col = ann_colors,
  border = T,
  annotation_name_gp = gpar(fontsize = 10)    # make this equal to one heatmap column width
)


# Remove "s__" at the start
rownames(mat_scaled) <- gsub("^s__", "", rownames(mat_scaled))


# ha_col_feat = ComplexHeatmap::HeatmapAnnotation(Features=annot_col_feat$Features, annotation_name_side = "left", show_annotation_name = F, show_legend=c(T))
# Build heatmap
set.seed(100)
ht_vir <- Heatmap(
  mat_scaled,
  name = "zscore",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = T,
  show_column_names = TRUE,
  right_annotation = ha,
  border = T,
  row_split = 3,
  row_title = NULL,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  width = unit(ncol(mat_scaled), "cm"))
ht_vir


pdf('Contrasted_vir_baseline.pdf', height = 10, width = 15)
print(ht_vir)
dev.off()


##################################################################
##                    Contrasted metabolites                    ##
##################################################################
metabolome <- read.csv('../DIRECT_baseline_metabolome.txt', sep = ' ', check.names = F)
metabolome_anno <- read.csv('../DIRECT_baseline_metabolome_annotation.txt', sep = ' ')
rownames(metabolome_anno) <- metabolome_anno$CHEM_ID
meta_lifestyle <- read.delim('../Meta_data/Lifestyle_score.txt', header = T, row.names = 1)
metabolome$'rowname' <- meta_lifestyle[rownames(meta_lifestyle), 'SampleID']
rownames(metabolome) <- metabolome$rowname
metabolome <- metabolome[, -ncol(metabolome)]
colnames(metabolome) <- gsub("^X", "", colnames(metabolome))
met_dt_fl <- metabolome[rownames(meta_fl), ]
colnames(met_dt_fl) <- metabolome_anno[colnames(met_dt_fl), 'CHEMICAL_NAME']

# Ensure Group is a factor
metabolite_dif_df <- cbind(met_dt_fl[rownames(meta_fl), ], meta_fl$Sex, meta_fl$bsl_age, meta_fl$bsl_centerID, meta_fl$Gly.Cat.48m)
colnames(metabolite_dif_df)[974:977] <- c('Sex', 'Age', 'bsl_centerID', 'Gly.Cat.48m')

metabolite_dif_df$Gly.Cat.48m <- factor(metabolite_dif_df$Gly.Cat.48m, levels = c("NGR", "IGR", "T2D"))

# Function to run ANOVA + Tukey HSD for one feature
run_anova_tukey <- function(feature, data) {
  formula <- as.formula(paste0("`", feature, "` ~ Gly.Cat.48m + Age + Sex + bsl_centerID"))
  fit <- aov(formula, data = data)
  
  # ANOVA table
  anova_res <- tidy(anova(fit))
  
  # Tukey HSD only if Group is significant
  tukey_res <- TukeyHSD(fit, "Gly.Cat.48m") %>%
    broom::tidy() %>%
    mutate(feature = feature)
  
  return(list(anova = anova_res, tukey = tukey_res))
}

# Apply across all microbiome features
metabolite_features <- colnames(metabolite_dif_df)[!(colnames(metabolite_dif_df) %in% c("Gly.Cat.48m","Age","Sex", 'bsl_centerID'))]

metabolite_results <- lapply(metabolite_features, run_anova_tukey, data = metabolite_dif_df)

# Combine into dataframes
metabolite_anova_results <- bind_rows(lapply(metabolite_results, function(x) x$anova), .id = "feature_id")
metabolite_tukey_results <- bind_rows(lapply(metabolite_results, function(x) x$tukey))

# Filter Tukey results for significant pairwise comparisons
metabolite_tukey_sig <- metabolite_tukey_results %>% filter(adj.p.value < 0.05)


library(dplyr)
library(tidyr)
library(pheatmap)
library(tibble)

# -------------------------------
# Step 1. Compute group means
# -------------------------------
metabolite_group_means <- metabolite_dif_df %>%
  dplyr::select(-Age, -Sex, -bsl_centerID, -Gly.Cat.48m) %>%  # drop metadata + group var
  bind_cols(Gly.Cat.48m = metabolite_dif_df$Gly.Cat.48m) %>%             # add group var back
  group_by(Gly.Cat.48m) %>%
  summarise(across(everything(), mean, na.rm = TRUE), .groups = "drop") %>%
  pivot_longer(-Gly.Cat.48m, names_to = "feature", values_to = "mean") %>%
  pivot_wider(names_from = Gly.Cat.48m, values_from = mean)

# -------------------------------
# Step 1b. Keep only features in tukey_sig
# -------------------------------
metabolite_sig_features <- unique(metabolite_tukey_sig$feature)

metabolite_group_means <- metabolite_group_means %>%
  filter(feature %in% metabolite_sig_features)

# Convert to matrix
metabolite_mat <- as.matrix(metabolite_group_means[, -1])
rownames(metabolite_mat) <- metabolite_group_means$feature

# Scale to z-scores by row
metabolite_mat_scaled <- t(scale(t(metabolite_mat)))

# -------------------------------
# Step 2. Create annotation table
# -------------------------------
metabolite_sig_contrasts <- metabolite_tukey_sig %>%
  filter(adj.p.value < 0.1) %>%
  dplyr::select(feature, contrast)

metabolite_anno <- metabolite_sig_contrasts %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = contrast, values_from = value, values_fill = 0) %>%
  right_join(data.frame(feature = rownames(metabolite_mat_scaled)), by = "feature") %>%
  replace(is.na(.), 0)

metabolite_row_anno <- metabolite_anno %>% column_to_rownames("feature")

rownames(metabolome_anno) <- metabolome_anno$CHEMICAL_NAME
metabolite_row_anno$'Functional category' <- metabolome_anno[rownames(metabolite_row_anno), 'SUPER_PATHWAY']


# Define heatmap colors
col_fun <- colorRamp2(c(-1.5, 0, 1.5), c("#2e756e", "white", "#945e25"))

# Convert row annotation to data frame
metabolite_row_anno_df <- as.data.frame(metabolite_row_anno)

# Define annotation colors
ann_colors <- list(
  `IGR-NGR` = c("0" = "white", "1" = "#3cb4e9"),
  `T2D-NGR` = c("0" = "white", "1" = "#f486ad"),
  `T2D-IGR` = c("0" = "white", "1" = "#ecce88")
)

# Row annotation object with fixed width
metabolite_ha <- rowAnnotation(
  df = metabolite_row_anno_df,
  col = ann_colors,
  border = T,
  annotation_name_gp = gpar(fontsize = 10)    # make this equal to one heatmap column width
)

# Build heatmap
set.seed(100)
ht_metabolite <- Heatmap(
  metabolite_mat_scaled,
  name = "z score",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = TRUE,
  right_annotation = metabolite_ha,
  border = T,
  # row_split = 3,
  row_title = NULL,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  width = unit(ncol(metabolite_mat_scaled), "cm"))
ht_metabolite



#################################################################
##        Make circular plot for contrasted metabolites        ##
#################################################################
# Load libraries
library(circlize)

# -------------------------------
# Step 1: Cluster rows of hp_mat
# -------------------------------
dist_mat <- dist(hp_mat[, 1:3])              # distance on metabolite z-scores
hc <- hclust(dist_mat, method = "ward.D2")   # hierarchical clustering
row_order <- hc$order                        # extract row order

# Step 2: Reorder hp_mat by clustering
hp_mat_ordered <- hp_mat[row_order, ]

# -------------------------------
# Step 3: Define color mappings
# -------------------------------
col_fun  <- colorRamp2(c(-1.5, 0, 1.5), c("#2e756e", "white", "#945e25"))
col_fun1 <- colorRamp2(c(0, 1), c("white", "#3cb4e9"))
col_fun2 <- colorRamp2(c(0, 1), c("white", "#f486ad"))
col_fun3 <- colorRamp2(c(0, 1), c("white", "#ecce88"))

# -------------------------------
# Step 4: Initialize circos
# -------------------------------
pdf('Contrasted_metabolites_baseline.pdf', height = 6, width = 6)
circos.clear()
circos.par(
  start.degree = 90,
  gap.degree = 90,
  track.margin = c(0.01, 0.01),
  clock.wise = FALSE
)

# -------------------------------
# Step 5: Add annotation tracks (contrasts)
# -------------------------------
circos.heatmap(
  hp_mat_ordered[, 4, drop = FALSE],
  col = col_fun1,
  # dend.side = "outside",
  # rownames.side = "outside",
  cluster = FALSE,
  track.height = 0.05
)
circos.track(track.index = get.current.track.index(), bg.border = "black")

circos.heatmap(
  hp_mat_ordered[, 5, drop = FALSE],
  col = col_fun2,
  # dend.side = "outside",
  cluster = FALSE,
  track.height = 0.05
)
circos.track(track.index = get.current.track.index(), bg.border = "black")

circos.heatmap(
  hp_mat_ordered[, 6, drop = FALSE],
  col = col_fun3,
  # dend.side = "outside",
  cluster = FALSE,
  track.height = 0.05
)
circos.track(track.index = get.current.track.index(), bg.border = "black")

# -------------------------------
# Step 6: Add metabolite z-score layer (inner)
# -------------------------------
circos.heatmap(
  hp_mat_ordered[, 1:3, drop = FALSE],
  col = col_fun,
  cluster = TRUE,
  track.height = 0.3
)
circos.track(track.index = get.current.track.index(), bg.border = "black")



# Continuous legend for metabolite z-scores
lgd_zscore <- Legend(
  col_fun = col_fun,
  title = "z score",
  title_position = "topcenter",
  legend_height = unit(4, "cm"),
  direction = "horizontal"
)

# Binary legends for contrasts
lgd_igr <- Legend(
  labels = c( "IGR versus NGR"),
  legend_gp = gpar(fill = c("#3cb4e9")),
  direction = "horizontal"
)

lgd_t2dngr <- Legend(
  labels = c("T2D versus NGR"),
  legend_gp = gpar(fill = c( "#f486ad")),
  direction = "horizontal"
)

lgd_t2digr <- Legend(
  labels = c( "T2D versus IGR"),
  legend_gp = gpar(fill = c( "#ecce88")),
  direction = "horizontal"
)

# Pack and draw all legends on the right
# lgd_list_category <- packLegend(lgd_igr, lgd_t2dngr, lgd_t2digr, direction = "vertical")
lgd_list <- packLegend(lgd_zscore, lgd_igr, lgd_t2dngr, lgd_t2digr, direction = "vertical")
# Draw at the bottom center
draw(lgd_list, x = unit(0.7, "npc"), y = unit(0.95, "npc"), just = c("center", "right"))

dev.off()

#################################################################
##  Identify features associated with T2D compared to NGR/IGR  ##
#################################################################

dif_df$'bsl_cell_count' <- meta_fl[rownames(dif_df), 'bsl_cell_count']
vir_dif_df$'bsl_cell_count' <- meta_fl[rownames(vir_dif_df), 'bsl_cell_count']

# Bacterial species
dif_df <- dif_df %>%
  mutate(T2D_status = if_else(Gly.Cat.48m == "T2D", TRUE, FALSE))

dif_df <- dif_df %>%
  mutate(IGR_T2D_status = if_else(Gly.Cat.48m == "NGR", FALSE, TRUE))

vir_dif_df <- vir_dif_df %>%
  mutate(T2D_status = if_else(Gly.Cat.48m == "T2D", TRUE, FALSE))

vir_dif_df <- vir_dif_df %>%
  mutate(IGR_T2D_status = if_else(Gly.Cat.48m == "NGR", FALSE, TRUE))

metabolite_dif_df <- metabolite_dif_df %>%
  mutate(T2D_status = if_else(Gly.Cat.48m == "T2D", TRUE, FALSE))

metabolite_dif_df <- metabolite_dif_df %>%
  mutate(IGR_T2D_status = if_else(Gly.Cat.48m == "NGR", FALSE, TRUE))


bac_48m_res <- Maaslin2(input_data     = bac_mat_filt[rownames(dif_df),], 
                        input_metadata = dif_df, 
                        analysis_method = 'LM', 
                        min_prevalence = 0.1,
                        min_abundance = 0,
                        normalization  = "NONE",
                        transform = 'NONE',
                        output         = "met_48m_res", 
                        fixed_effects  =  c("IGR_T2D_status", 'Age', 'Sex','bsl_centerID', 'bsl_cell_count'),
                        reference = "bsl_centerID,Center_13",
                        max_significance = 0.1,
                        standardize = T,
                        plot_scatter = F)
sig_bac <- bac_48m_res$results %>% filter(value == 'TRUE' & pval < 0.05)

sig_bac <- sig_bac %>%
  mutate(Significance = case_when(
    pval < 0.05 & qval < 0.1 ~ "qvalue < 0.1",
    pval < 0.05 & qval >= 0.1 ~ "pvalue < 0.05",
    TRUE ~ "Not significant"
  ))

# Visualization

# Calculate the confidence intervals
sig_bac$conf_low <- sig_bac$coef - 1.96 * sig_bac$stderr
sig_bac$conf_high <- sig_bac$coef + 1.96 * sig_bac$stderr

# sig_bac$'bac_taxa' <- bac_meta[sig_bac$feature, 2]

# sig_bac <- sig_bac[1:50, ]
sig_bac <- sig_bac[order(-sig_bac$coef), ]

sig_bac$feature <- sig_bac$feature %>%
  gsub("^s__", "", .) %>%            # remove s__ at the start
  gsub("\\.t__.*$", "", .) %>%       # remove .t__XXXX at the end
  gsub("_", " ", .)                  # replace _ with space


p_bac_coeff <- ggplot(sig_bac, aes(y = reorder(feature, coef), x = coef, xmin = conf_low, xmax = conf_high)) +
  geom_errorbar(aes(color = Significance), width=0, linewidth=1) +
  geom_point(aes(color = Significance), shape = 23, size = 2, fill = "white", stroke = 1.2) +
  scale_color_manual(values = c("qvalue < 0.1" = "orange", "pvalue < 0.05" = "dodgerblue")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.4) +
  # geom_rect(aes(xmin = -Inf, xmax = 1, ymin = -Inf, ymax = Inf), fill = 'NA', color = 'black') +
  # theme_minimal() +
  # scale_x_continuous(breaks = pretty(sig_res$coef, n = 10)) + # X axis ticks
  # scale_y_discrete(limits = rev(levels(reorder(sig_res$feature, sig_res$coef)))) +
  labs(x = "Effect Size (IGR+T2D versus NGR)", y = NA) +
  theme_bw() +
  theme(panel.border = element_rect(size=1), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_line(colour = "black", size = 1),
        axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10, face = 'italic'),
        axis.title.x = element_text(color = "black", size = 10),
        axis.title.y=element_blank(),
        legend.position = c(0.75, 0.2),
        legend.text = element_text(colour = "black", size = 8),
        legend.title = element_text(colour = "black", size = 10))
pdf('Bac_species_IGR_T2D_vs_NGR.pdf', height = 3.5, width = 5)
print(p_bac_coeff)
dev.off()




bac_48m_res <- Maaslin2(input_data     = bac_mat_filt[rownames(dif_df),], 
                        input_metadata = dif_df, 
                        analysis_method = 'LM', 
                        min_prevalence = 0.1,
                        min_abundance = 0,
                        normalization  = "NONE",
                        transform = 'NONE',
                        output         = "met_48m_res", 
                        fixed_effects  =  c("T2D_status", 'Age', 'Sex','bsl_centerID', 'bsl_cell_count'),
                        reference = "bsl_centerID,Center_13",
                        max_significance = 0.1,
                        standardize = T,
                        plot_scatter = F)
sig_bac <- bac_48m_res$results %>% filter(value == 'TRUE' & pval < 0.05)

sig_bac <- sig_bac %>%
  mutate(Significance = case_when(
    pval < 0.05 & qval < 0.1 ~ "qvalue < 0.1",
    pval < 0.05 & qval >= 0.1 ~ "pvalue < 0.05",
    TRUE ~ "Not significant"
  ))

# Visualization

# Calculate the confidence intervals
sig_bac$conf_low <- sig_bac$coef - 1.96 * sig_bac$stderr
sig_bac$conf_high <- sig_bac$coef + 1.96 * sig_bac$stderr

# sig_bac$'bac_taxa' <- bac_meta[sig_bac$feature, 2]

# sig_bac <- sig_bac[1:50, ]
sig_bac <- sig_bac[order(-sig_bac$coef), ]

sig_bac$feature <- sig_bac$feature %>%
  gsub("^s__", "", .) %>%            # remove s__ at the start
  gsub("\\.t__.*$", "", .) %>%       # remove .t__XXXX at the end
  gsub("_", " ", .)                  # replace _ with space


p_bac_T2D_vs_IGR_NGR <- ggplot(sig_bac, aes(y = reorder(feature, coef), x = coef, xmin = conf_low, xmax = conf_high)) +
  geom_errorbar(aes(color = Significance), width=0, linewidth=1) +
  geom_point(aes(color = Significance), shape = 23, size = 2, fill = "white", stroke = 1.2) +
  scale_color_manual(values = c("qvalue < 0.1" = "orange", "pvalue < 0.05" = "dodgerblue")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.4) +
  # geom_rect(aes(xmin = -Inf, xmax = 1, ymin = -Inf, ymax = Inf), fill = 'NA', color = 'black') +
  # theme_minimal() +
  # scale_x_continuous(breaks = pretty(sig_res$coef, n = 10)) + # X axis ticks
  # scale_y_discrete(limits = rev(levels(reorder(sig_res$feature, sig_res$coef)))) +
  labs(x = "Effect Size (T2D versus IGR+NGR)", y = NA) +
  theme_bw() +
  theme(panel.border = element_rect(size=1), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_line(colour = "black", size = 1),
        axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10, face = 'italic'),
        axis.title.x = element_text(color = "black", size = 10),
        axis.title.y=element_blank(),
        legend.position = c(0.75, 0.2),
        legend.text = element_text(colour = "black", size = 8),
        legend.title = element_text(colour = "black", size = 10))
pdf('Bac_species_T2D_vs_IGR_NGR.pdf', height = 3.5, width = 5)
print(p_bac_T2D_vs_IGR_NGR)
dev.off()



vir_48m_res <- Maaslin2(input_data     = vir_mat_filt[rownames(dif_df),], 
                        input_metadata = dif_df, 
                        analysis_method = 'LM', 
                        min_prevalence = 0.1,
                        min_abundance = 0,
                        normalization  = "NONE",
                        transform = 'NONE',
                        output         = "met_48m_res", 
                        fixed_effects  =  c("IGR_T2D_status", 'Age', 'Sex','bsl_centerID', 'bsl_cell_count'),
                        reference = "bsl_centerID,Center_13",
                        max_significance = 0.1,
                        standardize = T,
                        plot_scatter = F)
sig_vir <- vir_48m_res$results %>% filter(value == 'TRUE' & pval < 0.05)

sig_vir <- sig_vir %>%
  mutate(Significance = case_when(
    pval < 0.05 & qval < 0.1 ~ "qvalue < 0.1",
    pval < 0.05 & qval >= 0.1 ~ "pvalue < 0.05",
    TRUE ~ "Not significant"
  ))

# Visualization

# Calculate the confidence intervals
sig_vir$conf_low <- sig_vir$coef - 1.96 * sig_vir$stderr
sig_vir$conf_high <- sig_vir$coef + 1.96 * sig_vir$stderr

# sig_bac$'bac_taxa' <- bac_meta[sig_bac$feature, 2]

# sig_bac <- sig_bac[1:50, ]
sig_vir <- sig_vir[order(-sig_vir$coef), ]

sig_vir$feature <- sig_vir$feature %>%
  gsub("^s__", "", .)          # remove s__ at the start

p_vir_T2D_IGR_vs_NGR <- ggplot(sig_vir, aes(y = reorder(feature, coef), x = coef, xmin = conf_low, xmax = conf_high)) +
  geom_errorbar(aes(color = Significance), width=0, linewidth=1) +
  geom_point(aes(color = Significance), shape = 23, size = 2, fill = "white", stroke = 1.2) +
  scale_color_manual(values = c("qvalue < 0.1" = "orange", "pvalue < 0.05" = "dodgerblue")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.4) +
  # geom_rect(aes(xmin = -Inf, xmax = 1, ymin = -Inf, ymax = Inf), fill = 'NA', color = 'black') +
  # theme_minimal() +
  # scale_x_continuous(breaks = pretty(sig_res$coef, n = 10)) + # X axis ticks
  # scale_y_discrete(limits = rev(levels(reorder(sig_res$feature, sig_res$coef)))) +
  labs(x = "Effect Size (IGR+T2D versus NGR)", y = NA) +
  theme_bw() +
  theme(panel.border = element_rect(size=1), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_line(colour = "black", size = 1),
        axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10, face = 'italic'),
        axis.title.x = element_text(color = "black", size = 10),
        axis.title.y=element_blank(),
        legend.position = c(0.75, 0.2),
        legend.text = element_text(colour = "black", size = 8),
        legend.title = element_text(colour = "black", size = 10))
pdf('Vir_species_IGR_T2D_vs_NGR.pdf', height = 3.5, width = 5)
print(p_vir_T2D_IGR_vs_NGR)
dev.off()


vir_48m_res <- Maaslin2(input_data     = vir_mat_filt[rownames(dif_df),], 
                        input_metadata = dif_df, 
                        analysis_method = 'LM', 
                        min_prevalence = 0.1,
                        min_abundance = 0,
                        normalization  = "NONE",
                        transform = 'NONE',
                        output         = "met_48m_res", 
                        fixed_effects  =  c("T2D_status", 'Age', 'Sex','bsl_centerID', 'bsl_cell_count'),
                        reference = "bsl_centerID,Center_13",
                        max_significance = 0.1,
                        standardize = T,
                        plot_scatter = F)
sig_vir <- vir_48m_res$results %>% filter(value == 'TRUE' & pval < 0.05)

sig_vir <- sig_vir %>%
  mutate(Significance = case_when(
    pval < 0.05 & qval < 0.1 ~ "qvalue < 0.1",
    pval < 0.05 & qval >= 0.1 ~ "pvalue < 0.05",
    TRUE ~ "Not significant"
  ))

# Visualization

# Calculate the confidence intervals
sig_vir$conf_low <- sig_vir$coef - 1.96 * sig_vir$stderr
sig_vir$conf_high <- sig_vir$coef + 1.96 * sig_vir$stderr

# sig_bac$'bac_taxa' <- bac_meta[sig_bac$feature, 2]

# sig_bac <- sig_bac[1:50, ]
sig_vir <- sig_vir[order(-sig_vir$coef), ]

sig_vir$feature <- sig_vir$feature %>%
  gsub("^s__", "", .)          # remove s__ at the start


p_vir_T2D_vs_IGR_NGR <- ggplot(sig_vir, aes(y = reorder(feature, coef), x = coef, xmin = conf_low, xmax = conf_high)) +
  geom_errorbar(aes(color = Significance), width=0, linewidth=1) +
  geom_point(aes(color = Significance), shape = 23, size = 2, fill = "white", stroke = 1.2) +
  scale_color_manual(values = c("qvalue < 0.1" = "orange", "pvalue < 0.05" = "dodgerblue")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.4) +
  # geom_rect(aes(xmin = -Inf, xmax = 1, ymin = -Inf, ymax = Inf), fill = 'NA', color = 'black') +
  # theme_minimal() +
  # scale_x_continuous(breaks = pretty(sig_res$coef, n = 10)) + # X axis ticks
  # scale_y_discrete(limits = rev(levels(reorder(sig_res$feature, sig_res$coef)))) +
  labs(x = "Effect Size (T2D versus IGR+NGR)", y = NA) +
  theme_bw() +
  theme(panel.border = element_rect(size=1), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_line(colour = "black", size = 1),
        axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10, face = 'italic'),
        axis.title.x = element_text(color = "black", size = 10),
        axis.title.y=element_blank(),
        legend.position = c(0.75, 0.2),
        legend.text = element_text(colour = "black", size = 8),
        legend.title = element_text(colour = "black", size = 10))
pdf('Bac_vir_T2D_vs_IGR_NGR.pdf', height = 3.5, width = 5)
print(p_vir_T2D_vs_IGR_NGR)
dev.off()



metabolome <- read.csv('../DIRECT_baseline_metabolome.txt', sep = ' ', check.names = F)
metabolome_anno <- read.csv('../DIRECT_baseline_metabolome_annotation.txt', sep = ' ')
rownames(metabolome_anno) <- metabolome_anno$CHEM_ID
meta_lifestyle <- read.delim('../Meta_data/Lifestyle_score.txt', header = T, row.names = 1)
metabolome$'rowname' <- meta_lifestyle[rownames(meta_lifestyle), 'SampleID']
rownames(metabolome) <- metabolome$rowname
metabolome <- metabolome[, -ncol(metabolome)]
colnames(metabolome) <- gsub("^X", "", colnames(metabolome))
met_dt_fl <- metabolome[rownames(meta_fl), ]
# colnames(met_dt_fl) <- metabolome_anno[colnames(met_dt_fl), 'CHEMICAL_NAME']

met_48m_res <- Maaslin2(input_data     = met_dt_fl[rownames(dif_df),], 
                        input_metadata = dif_df, 
                        analysis_method = 'LM', 
                        min_prevalence = 0.1,
                        min_abundance = 0,
                        normalization  = "NONE",
                        transform = 'LOG',
                        output         = "met_48m_res", 
                        fixed_effects  =  c("IGR_T2D_status", 'Age', 'Sex','bsl_centerID', 'bsl_cell_count'),
                        reference = "bsl_centerID,Center_13",
                        max_significance = 0.1,
                        standardize = T,
                        plot_scatter = F)
sig_met <- met_48m_res$results %>% filter(value == 'TRUE' & qval < 0.05)

# Calculate the confidence intervals
sig_met$conf_low <- sig_met$coef - 1.96 * sig_met$stderr
sig_met$conf_high <- sig_met$coef + 1.96 * sig_met$stderr

sig_met$feature <- sig_met$feature %>% gsub("^X", "", .)

sig_met$'Metabolite' <- metabolome_anno[sig_met$feature, 'CHEMICAL_NAME']

# sig_bac$'bac_taxa' <- bac_meta[sig_bac$feature, 2]

# sig_bac <- sig_bac[1:50, ]
sig_met <- sig_met[order(-sig_met$coef), ]

# sig_vir$feature <- sig_vir$feature %>%
#   gsub("^s__", "", .)          # remove s__ at the start

p_metabolite_T2D_IGR_vs_NGR <- ggplot(sig_met, aes(y = reorder(Metabolite, coef), x = coef, xmin = conf_low, xmax = conf_high)) +
  geom_errorbar(aes(color = -log10(qval)), width=0, linewidth=0.5) +
  geom_point(aes(color = -log10(qval)), shape = 23, size = 1, fill = "white", stroke = 1) +
  scale_color_gradient(low = "dodgerblue", high = "orange") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.4) +
  # geom_rect(aes(xmin = -Inf, xmax = 1, ymin = -Inf, ymax = Inf), fill = 'NA', color = 'black') +
  # theme_minimal() +
  # scale_x_continuous(breaks = pretty(sig_res$coef, n = 10)) + # X axis ticks
  # scale_y_discrete(limits = rev(levels(reorder(sig_res$feature, sig_res$coef)))) +
  labs(x = "Effect Size (IGR+T2D versus NGR)", y = NA) +
  theme_bw() +
  theme(panel.border = element_rect(size=1), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_line(colour = "black", size = 1),
        axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        axis.title.x = element_text(color = "black", size = 10),
        axis.title.y=element_blank(),
        legend.position = c(0.8,0.1),
        legend.text = element_text(colour = "black", size = 8),
        legend.title = element_text(colour = "black", size = 10))


p_metabolite_T2D_IGR_vs_NGR

pdf('Metabolite_T2D_IGR_vs_NGR.pdf', height = 12, width = 7)
print(p_metabolite_T2D_IGR_vs_NGR)
dev.off()







met_48m_res <- Maaslin2(input_data     = met_dt_fl[rownames(dif_df),], 
                        input_metadata = dif_df, 
                        analysis_method = 'LM', 
                        min_prevalence = 0.1,
                        min_abundance = 0,
                        normalization  = "NONE",
                        transform = 'LOG',
                        output         = "met_48m_res", 
                        fixed_effects  =  c("T2D_status", 'Age', 'Sex','bsl_centerID', 'bsl_cell_count'),
                        reference = "bsl_centerID,Center_13",
                        max_significance = 0.1,
                        standardize = T,
                        plot_scatter = F)
sig_met <- met_48m_res$results %>% filter(value == 'TRUE' & qval < 0.05)

# Calculate the confidence intervals
sig_met$conf_low <- sig_met$coef - 1.96 * sig_met$stderr
sig_met$conf_high <- sig_met$coef + 1.96 * sig_met$stderr

sig_met$feature <- sig_met$feature %>% gsub("^X", "", .)

sig_met$'Metabolite' <- metabolome_anno[sig_met$feature, 'CHEMICAL_NAME']

# sig_bac$'bac_taxa' <- bac_meta[sig_bac$feature, 2]

# sig_bac <- sig_bac[1:50, ]
sig_met <- sig_met[order(-sig_met$coef), ]

# sig_vir$feature <- sig_vir$feature %>%
#   gsub("^s__", "", .)          # remove s__ at the start

p_metabolite_T2D_vs_IGR_NGR <- ggplot(sig_met, aes(y = reorder(Metabolite, coef), x = coef, xmin = conf_low, xmax = conf_high)) +
  geom_errorbar(aes(color = -log10(qval)), width=0, linewidth=0.5) +
  geom_point(aes(color = -log10(qval)), shape = 23, size = 1, fill = "white", stroke = 1) +
  scale_color_gradient(low = "dodgerblue", high = "orange") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.4) +
  # geom_rect(aes(xmin = -Inf, xmax = 1, ymin = -Inf, ymax = Inf), fill = 'NA', color = 'black') +
  # theme_minimal() +
  # scale_x_continuous(breaks = pretty(sig_res$coef, n = 10)) + # X axis ticks
  # scale_y_discrete(limits = rev(levels(reorder(sig_res$feature, sig_res$coef)))) +
  labs(x = "Effect Size (T2D versus IGR+NGR)", y = NA) +
  theme_bw() +
  theme(panel.border = element_rect(size=1), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_line(colour = "black", size = 1),
        axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        axis.title.x = element_text(color = "black", size = 10),
        axis.title.y=element_blank(),
        legend.position = c(0.8,0.1),
        legend.text = element_text(colour = "black", size = 8),
        legend.title = element_text(colour = "black", size = 10))


p_metabolite_T2D_vs_IGR_NGR

pdf('Metabolite_T2D_vs_IGR_NGR.pdf', height = 14, width = 7)
print(p_metabolite_T2D_vs_IGR_NGR)
dev.off()

# Figure 5b
# Visualization of the prediction performance of different data types

IGR_T2D <- read.delim('Inputs/Iterative_AUC_IGR_vs_T2D.txt', header = T, row.names = 1, check.names = F)
NGR_IGR <- read.delim('Inputs/Iterative_AUC_NGR_vs_T2D.txt', header = T, row.names = 1, check.names = F)
NGR_T2D <- read.delim('Inputs/Iterative_AUC_NGR_vs_T2D.txt', header = T, row.names = 1, check.names = F)

prediction_column_names <- c('Risk factors',
                             'Plasma metabolites',
                             'Bacterial species',
                             'Bacterial subspecies',
                             'Viral species',
                             'Risk factors + Plasma metabolites',
                             'Risk factors + Bacterial species',
                             'Risk factors + Bacterial subspecies',
                             'Risk factors + Viral species')

colnames(IGR_T2D) <- prediction_column_names
colnames(NGR_IGR) <- prediction_column_names
colnames(NGR_T2D) <- prediction_column_names

prediction_order <- c(
  "Risk factors",
  "Plasma metabolites",
  "Bacterial species",
  "Bacterial subspecies",
  "Viral species",
  "Risk factors + Plasma metabolites",
  "Risk factors + Bacterial species",
  "Risk factors + Bacterial subspecies",
  "Risk factors + Viral species"
)

plot_iterative_auc_box_jitter <- function(
  auc_df,
  comparison_name = "",
  x_label = "AUC",
  jitter_width = 0.15,
  model_order
) {
  
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(egg)
  
  ## ---- reshape to long format ----
  auc_long <- auc_df %>%
    rownames_to_column("Iteration") %>%
    pivot_longer(
      cols = -Iteration,
      names_to = "Model",
      values_to = "AUC"
    ) %>%
    mutate(
      Model = factor(Model, levels = rev(model_order))
    )
  
  ## ---- plot ----
  ggplot(
    auc_long,
    aes(x = Model, y = AUC)
  ) +
    geom_boxplot(
      outlier.shape = NA,
      width = 0.6,
      fill = "grey85",
      color = "black"
    ) +
    geom_jitter(
      width = jitter_width,
      size = 1.8,
      alpha = 0.7
    ) +
    coord_flip() +
    labs(
      title = comparison_name,
      x = NULL,
      y = x_label
    ) +
    theme_article(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.y = element_text(size = 11)
    )
}


p_IGR_T2D <- plot_iterative_auc_box_jitter(
  IGR_T2D,
  comparison_name = "IGR vs T2D",
  model_order = prediction_order
)

p_NGR_IGR <- plot_iterative_auc_box_jitter(
  NGR_IGR,
  comparison_name = "NGR vs IGR",
  model_order = prediction_order
)

p_NGR_T2D <- plot_iterative_auc_box_jitter(
  NGR_T2D,
  comparison_name = "NGR vs T2D",
  model_order = prediction_order
)

library(egg)

plot_iteration_AUCs <- ggarrange(
  p_NGR_IGR,
  p_IGR_T2D,
  p_NGR_T2D,
  ncol = 3,
  labels = c("a", "b", "c")
)


pdf('Outputs/ML_iteration_AUCs.pdf', width = 15, height = 6)
print(plot_iteration_AUCs)
dev.off()
