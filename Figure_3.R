# Figure 3a
# OR for gene richness (binary)
library(dplyr)
library(broom)

df <- read.delim('/Users/fjw536/Library/CloudStorage/OneDrive-UniversityofCopenhagen/Desktop/DIRECT/DIRECT_reanalysis/gene_richness_df.txt', header = T, row.names = 1, check.names = F) # Load the dataset and reclassifiy the LGC and HGC

# OR for gene richness (continuous, per SD)

df <- df %>%
  mutate(
    gene_richness_z = scale(bsl_gene_richness)
  )



model_gr_cont_un <- glm(
  T2D_48m ~ gene_richness_z, 
  data = df,
  family = binomial
)

model_gr_cont_un <- tidy(model_gr_cont_un) %>%
  filter(term == "gene_richness_z") %>%
  mutate(
    OR = exp(estimate),
    CI_low = exp(estimate - 1.96 * std.error),
    CI_high = exp(estimate + 1.96 * std.error)
  ) %>%
  dplyr::select(OR, CI_low, CI_high, p.value)

model_gr_cont_un

model_gr_cont_sex <- glm(
  T2D_48m ~ gene_richness_z + Sex, 
  data = df,
  family = binomial
)

or_gr_cont_sex <- tidy(model_gr_cont_sex) %>%
  filter(term == "gene_richness_z") %>%
  mutate(
    OR = exp(estimate),
    CI_low = exp(estimate - 1.96 * std.error),
    CI_high = exp(estimate + 1.96 * std.error)
  ) %>%
  dplyr::select(OR, CI_low, CI_high, p.value)

or_gr_cont_sex

or_gr_cont_1 <- tidy(model_gr_cont_1) %>%
  filter(term == "gene_richness_z") %>%
  mutate(
    OR = exp(estimate),
    CI_low = exp(estimate - 1.96 * std.error),
    CI_high = exp(estimate + 1.96 * std.error)
  ) %>%
  dplyr::select(OR, CI_low, CI_high, p.value)


model_gr_cont_1 <- glm(
  T2D_48m ~ gene_richness_z + bsl_age + Sex + bsl_centerID, 
  data = df,
  family = binomial
)
or_gr_cont_1 <- tidy(model_gr_cont_1) %>%
  filter(term == "gene_richness_z") %>%
  mutate(
    OR = exp(estimate),
    CI_low = exp(estimate - 1.96 * std.error),
    CI_high = exp(estimate + 1.96 * std.error)
  ) %>%
  dplyr::select(OR, CI_low, CI_high, p.value)

or_gr_cont_1

model_gr_cont_2 <- glm(
  T2D_48m ~ gene_richness_z + bsl_age + Sex + bsl_centerID + bsl_cell_count + BMI + Waist.Hip,
  data = df,
  family = binomial
)
or_gr_cont_2 <- tidy(model_gr_cont_2) %>%
  filter(term == "gene_richness_z") %>%
  mutate(
    OR = exp(estimate),
    CI_low = exp(estimate - 1.96 * std.error),
    CI_high = exp(estimate + 1.96 * std.error)
  ) %>%
  dplyr::select(OR, CI_low, CI_high, p.value)

or_gr_cont_2

model_gr_cont_3 <- glm(
  T2D_48m ~ gene_richness_z + bsl_age + Sex + bsl_centerID + bsl_cell_count + BMI + Waist.Hip + lifestyle_score,
  data = df,
  family = binomial
)
or_gr_cont_3 <- tidy(model_gr_cont_3) %>%
  filter(term == "gene_richness_z") %>%
  mutate(
    OR = exp(estimate),
    CI_low = exp(estimate - 1.96 * std.error),
    CI_high = exp(estimate + 1.96 * std.error)
  ) %>%
  dplyr::select(OR, CI_low, CI_high, p.value)

or_gr_cont_3


library(broom)
library(dplyr)

extract_or <- function(model, label) {
  tidy(model) %>%
    filter(term == "gene_richness_z") %>%
    mutate(
      OR = exp(estimate),
      CI_low = exp(estimate - 1.96 * std.error),
      CI_high = exp(estimate + 1.96 * std.error),
      model = label
    ) %>%
    dplyr::select(model, OR, CI_low, CI_high, p.value)
}

or_all_models <- bind_rows(
  extract_or(model_gr_cont_un, "Unadjusted"),
  extract_or(model_gr_cont_sex, "Model 1"),
  extract_or(model_gr_cont_1, "Model 2"),
  extract_or(model_gr_cont_2, "Model 3"),
  extract_or(model_gr_cont_3, "Model 4")
)

plot_df <- or_all_models %>%
  mutate(
    model = factor(
      model,
      levels = rev(c(
        "Unadjusted",
        "Model 1",
        "Model 2",
        "Model 3"
      )
    )),
    label_or = sprintf("%.2f (%.2f–%.2f)", OR, CI_low, CI_high),
    label_p = ifelse(
      p.value < 0.001,
      "p < 0.001",
      paste0("p = ", formatC(p.value, format = "f", digits = 3))
    )
  )

library(ggplot2)

plot_OR <- ggplot(plot_df, aes(y = model, x = OR)) +
  
  # Reference line
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  
  # CI bars
  geom_errorbarh(
    aes(xmin = CI_low, xmax = CI_high),
    height = 0.22,
    size = 0.8,
    color = "black"
  ) +
  
  # Point estimates
  geom_point(
    size = 3.6,
    shape = 21,
    fill = "#2C7FB8",
    color = "black"
  ) +
  
  # OR (CI) text
  geom_text(
    aes(label = label_or),
    x = max(plot_df$CI_high) * 1.15,
    hjust = 0,
    size = 3.8
  ) +
  # Log-scaled x-axis (ORs remain ORs)
  scale_x_log10(
    breaks = c(0.5, 0.75, 1, 1.25),
    limits = c(
      min(plot_df$CI_low) * 0.8,
      max(plot_df$CI_high) * 1.5
    )
  ) +
  
  labs(
    x = "Odds ratio\n(per 1 SD increase in microbial gene family richness)",
    y = NULL
    # title = "Baseline microbial gene family richness and\n risk of incident type 2 diabetes",
    # subtitle = "Unadjusted and sequentially adjusted logistic regression models"
  ) +
  
  theme_linedraw(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(color = "black"),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    axis.text.y = element_text(size = 13, color = 'black'),
    plot.margin = margin(10, 90, 10, 10)
  )

pdf('Outputs/Plot_OR_baseline_microbial_generichness.pdf', height = 3, width = 5)
print(plot_OR)
dev.off()


# Figure 3c
# Load necessary libraries
library(NbClust)
library(fpc)
library(cluster)
library(effsize)
setwd('/Users/fjw536/OneDrive - University of Copenhagen/Desktop/DIRECT/DIRECT_reanalysis/')
lifestyle <- read.delim('/Users/fjw536/Library/CloudStorage/OneDrive-UniversityofCopenhagen/Desktop/DIRECT/DIRECT_reanalysis/Meta_data/Lifestyle_score.txt', header = T, row.names = 2)
# Load the data
data <- read.table("Meta_data/Meta_for_clustering.txt", header = TRUE, sep = "\t", row.names = 1)
data$'lifescore' <- lifestyle[rownames(data), 'LifestyleScore']

saveRDS(data, 'Meta_data/Meta_for_clustering.Rdata')

data <- meta_0m_775[, c(2, 10:30, 34:49)]

# Split data by gender
men_data <- data[data$Sex == "Male", ]
women_data <- data[data$Sex == "Female", ]

# Function to pre-process the data
preprocess_data <- function(data) {
  # Winsorization
  for (col in colnames(data)) {
    if (is.numeric(data[[col]])) {
      data[[col]][is.na(data[[col]])] <- median(data[[col]], na.rm = TRUE)  # Replace NAs with the median
      data[[col]] <- pmin(pmax(data[[col]], quantile(data[[col]], 0.01, na.rm = TRUE)), quantile(data[[col]], 0.99, na.rm = TRUE))
    }
  }
  
  # Log transformation for non-normally distributed variables
  for (col in colnames(data)) {
    if (is.numeric(data[[col]])) {
      shapiro_test <- shapiro.test(data[[col]])
      if (shapiro_test$p.value < 0.05) {
        min_val <- min(data[[col]], na.rm = TRUE)
        if (min_val <= 0) {
          data[[col]] <- data[[col]] + abs(min_val) + 1  # Shift values to be positive
        }
        data[[col]] <- log10(data[[col]] + 1)  # log10(x+1) to handle zero values
      }
    }
  }
  
  # Scaling
  for (col in colnames(data)) {
    if (is.numeric(data[[col]])) {
      data[[col]] <- scale(data[[col]])
    }
  }
  
  return(data)
}

# Pre-process the data, excluding the first column which is the gender
men_data <- preprocess_data(men_data[, 2:ncol(men_data)])
women_data <- preprocess_data(women_data[, 2:ncol(women_data)])


# Combine pre-processed data
combined_data <- rbind(men_data, women_data)

# Remove non-numeric columns for clustering
numeric_data <- combined_data[sapply(combined_data, is.numeric)]

# Function to determine the optimal number of clusters using a specific index with a limited number of clusters
determine_optimal_clusters <- function(data, index, min_clusters = 2, max_clusters = 5) {
  nbclust_result <- NbClust(data, distance = "euclidean", min.nc = min_clusters, max.nc = max_clusters, method = "kmeans", index = index)
  print(str(nbclust_result))  # Print the structure of the result
  # Check if Best.nc is a matrix or a list and handle accordingly
  if (is.matrix(nbclust_result$Best.nc)) {
    return(nbclust_result$Best.nc[1, "Number_clusters"])
  } else if (is.list(nbclust_result$Best.nc)) {
    return(nbclust_result$Best.nc$Number_clusters)
  } else if (is.numeric(nbclust_result$Best.nc)) {
    return(nbclust_result$Best.nc)
  } else {
    stop("Unexpected structure of NbClust result.")
  }
}


# Determine optimal number of clusters using DB index
set.seed(123)
optimal_clusters_db <- determine_optimal_clusters(numeric_data, "db", max_clusters = 5)

# Determine optimal number of clusters using KL index
set.seed(123)
optimal_clusters_kl <- determine_optimal_clusters(numeric_data, "kl", max_clusters = 5)


# Choose the optimal number of clusters based on the indices
best_num_clusters <- max(optimal_clusters_kl)
best_num_clusters <- 4

# Perform k-means clustering
set.seed(123)
kmeans_result <- kmeans(numeric_data, centers = best_num_clusters, nstart = 25)

# Evaluate cluster stability using the Jaccard similarity coefficient
cluster_stability <- clusterboot(numeric_data, clustermethod = kmeansCBI, k = best_num_clusters, B = 100)

# Check stability scores
stability_scores <- cluster_stability$bootmean
stable_clusters <- which(stability_scores > 0.8)

# Print results
print(paste("Number of clusters:", best_num_clusters))
print("Stability scores:")
print(stability_scores)
print("Stable clusters (stability score > 0.85):")
print(stable_clusters)

# Identify feature variables based on enrichment
for (i in 1:best_num_clusters) {
  cluster_data <- numeric_data[kmeans_result$cluster == i, ]
  for (col in colnames(cluster_data)) {
    if (is.numeric(cluster_data[[col]])) {
      wilcox_test <- wilcox.test(cluster_data[[col]], numeric_data[[col]])
      cliff_delta_result <- cliff.delta(cluster_data[[col]], numeric_data[[col]])
      effect_size <- cliff_delta_result$estimate
      if (wilcox_test$p.value < 0.05 && abs(effect_size) > 0) {
        print(paste("Cluster", i, "enriched variable:", col, "p-value:", wilcox_test$p.value, "Cliff's Delta effect size:", effect_size))
      }
    }
  }
}


numeric_data$'Cluster_data' <- kmeans_result$cluster[rownames(numeric_data)]
data$'Cluster_data' <- kmeans_result$cluster[rownames(data)]

numeric_data <- numeric_data %>%
  mutate(Cluster_data = case_when(
    Cluster_data == 1 ~ 2,
    Cluster_data == 2 ~ 4,
    Cluster_data == 3 ~ 1,
    Cluster_data == 4 ~ 3,
    TRUE ~ Cluster_data
  ))


data <- data %>%
  mutate(Cluster_data = case_when(
    Cluster_data == 1 ~ 2,
    Cluster_data == 2 ~ 4,
    Cluster_data == 3 ~ 1,
    Cluster_data == 4 ~ 3,
    TRUE ~ Cluster_data
  ))


# meta_fl$'bsl_gene_richness' <- meta_bsl[rownames(meta_fl), 'Gene_richness']





numeric_data <- numeric_data %>% mutate(Cluster_data = paste0('MC_', Cluster_data))
data <- data %>% mutate(Cluster_data = paste0('MC_', Cluster_data))

saveRDS(data, '../DIRECT_reanalysis/IMI_DIRECT_reanalysis_Jun_25_2025/Outputs/Meta_w_clustering.Rdata')


bsl_48m <- read.csv('Meta_data/WP2.1_775-48m_ClinVar_M.csv', sep = '\t', check.names = F)
rownames(bsl_48m) <- bsl_48m$SampleID
data$'Gly.Cat.48m' <- bsl_48m[rownames(data), 'Gly.Cat.48']

data$'rownames' <- rownames(data)

saveRDS(data, file = 'Meta_data/Unsupervise_clustering_baseline_775.Rdata')
# data <- readRDS('Meta_data/Unsupervise_clustering_baseline_775.Rdata')


library(ggplot2)
library(reshape2)
library(circlize)
library(ComplexHeatmap)

# Assuming the dataframe is saved as 'df' in R
# Convert the data frame to a matrix, excluding the 'Cluster_data' column for the heatmap
numeric_data <- numeric_data %>% arrange(Cluster_data)
df_matrix <- as.matrix(t(numeric_data[1:37]))
colnames(df_matrix) <- rownames(data)
df_matrix<- apply(df_matrix, 2, as.numeric)
rownames(df_matrix) <- colnames(numeric_data)[1:37]

# Create row annotations for Cluster_data
annotation <- data.frame(Cluster = factor(numeric_data$Cluster_data))
rownames(annotation) <- colnames(df_matrix)

# Create a color mapping for the clusters
# cluster_colors <- structure(c("#99c48a", "#87adce", "#f4b583", "#e3867f"), names = levels(annotation$Cluster))
cluster_colors <- structure(c("#99c48a", '#769bc5', "#f4b583", '#d8746a'), names = levels(annotation$Cluster))
# Create the top annotation object without showing the annotation name
top_annotation <- HeatmapAnnotation(
  df = annotation,
  col = list(Cluster = cluster_colors),
  show_annotation_name = FALSE,
  show_legend = FALSE
)



annotation_colors <- list(
  Cluster = c("MC_1" = "#99c48a", "MC_2" = "#769bc5", "MC_3" = "#f4b583", "MC_4" = "#d8746a")
)

# Generate the heatmap
col_fun <- colorRamp2(c(-2, 0, 2), c("#08699c", "white", "#c30a2f"))

heatmap_legend_param <- list(
  title_gp = gpar(fontsize = 12),  # Increase the font size of the title
  labels_gp = gpar(fontsize = 12),
  legend_height = unit(4, "cm"),
  title_position = 'topcenter'# Increase the font size of the labels
)

# df_matrix_scaled <- t(apply(df_matrix, 1, scale))
heatmap <- Heatmap(df_matrix,
                   name = "Scaled\n value",
                   col = col_fun,
                   show_column_names = FALSE,
                   # clustering_distance_rows = "kendall",
                   cluster_columns = FALSE,  # Disable automatic column clustering
                   cluster_rows = T,
                   show_row_dend = F,
                   top_annotation = top_annotation,
                   column_split = annotation$Cluster,
                   heatmap_legend_param = heatmap_legend_param)  # Apply the legend parameters
draw(heatmap, heatmap_legend_side = "right",  merge_legends = F)

pdf('/Users/fjw536/Library/CloudStorage/OneDrive-UniversityofCopenhagen/Desktop/DIRECT/DIRECT_reanalysis/IMI_DIRECT_reanalysis_Jun_25_2025/IMI_DIRECT_Figures/Outputs/Metabolic_cluster_heatmap.pdf', width = 7, height = 7)
print(heatmap)
dev.off()

# Figure 3d
#####################################################################################
##  High-risk metabolic clusters (MC_3 + MC_4) vs low-risk clusters (MC_1 + MC_2)  ##
##              overall and stratified by gene richness (low vs high).             ##
#####################################################################################

library(dplyr)
library(broom)
library(ggplot2)

df <- df %>%
  mutate(
    Cluster_highrisk = factor(
      Cluster_highrisk,
      levels = c("MC_1_2", "MC_3_4")  # MC_1_2 is reference
    ),
    Gene_richness_binary = factor(
      Gene_richness_binary,
      levels = c("Low gene richness", "High gene richness")
    ),
    Sex = factor(Sex)
  )


extract_or_cluster <- function(model, group_label) {
  tidy(model) %>%
    filter(term == "Cluster_highriskMC_3_4") %>%
    mutate(
      OR = exp(estimate),
      CI_low = exp(estimate - 1.96 * std.error),
      CI_high = exp(estimate + 1.96 * std.error),
      group = group_label
    ) %>%
    dplyr::select(group, OR, CI_low, CI_high, p.value)
}


model_all <- glm(
  T2D_48m ~ Cluster_highrisk + bsl_age + Sex + bsl_centerID,
  data = df,
  family = binomial
)

or_all <- extract_or_cluster(
  model_all,
  group_label = "Overall cohort"
)


df_low <- df %>% filter(Gene_richness_binary == "Low gene richness")

model_low <- glm(
  T2D_48m ~ Cluster_highrisk + bsl_age + Sex + bsl_centerID,
  data = df_low,
  family = binomial
)

or_low <- extract_or_cluster(
  model_low,
  group_label = "Low gene richness"
)


df_high <- df %>% filter(Gene_richness_binary == "High gene richness")

model_high <- glm(
  T2D_48m ~ Cluster_highrisk + bsl_age + Sex + bsl_centerID,
  data = df_high,
  family = binomial
)

or_high <- extract_or_cluster(
  model_high,
  group_label = "High gene richness"
)


or_summary <- bind_rows(or_all, or_low, or_high)

or_summary


plot_df <- or_summary %>%
  mutate(
    group = factor(
      group,
      levels = rev(c("Overall cohort", "Low gene richness", "High gene richness"))
    ),
    label_or = sprintf("%.2f (%.2f–%.2f)", OR, CI_low, CI_high),
    label_p = ifelse(
      p.value < 0.001,
      "p < 0.001",
      paste0("p = ", formatC(p.value, format = "f", digits = 3))
    )
  )


plot_OR_strat_gene_richness <- ggplot(plot_df, aes(y = group, x = OR)) +
  
  # Reference line
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  
  # CI bars
  geom_errorbarh(
    aes(xmin = CI_low, xmax = CI_high),
    height = 0.22,
    size = 0.9,
    color = "black"
  ) +
  
  # Point estimates
  geom_point(
    size = 4,
    shape = 21,
    fill = "#D73027",
    color = "black"
  ) +
  
  # OR (CI) text
  geom_text(
    aes(label = label_or),
    x = max(plot_df$CI_high) * 1.15,
    hjust = 0,
    size = 4
  ) +
  
  # p-value text
  geom_text(
    aes(label = label_p),
    x = max(plot_df$CI_high) * 1.7,
    hjust = 0,
    size = 4
  ) +
  
  # Log scale for ORs (ORs themselves unchanged)
  scale_x_log10(
    breaks = c(0.5, 1, 2, 5, 10),
    limits = c(
      min(plot_df$CI_low) * 0.8,
      max(plot_df$CI_high) * 2
    )
  ) +
  
  coord_cartesian(clip = "off") +
  
  labs(
    x = "Odds ratio for incident type 2 diabetes (MC_3–4 vs MC_1–2)",
    y = NULL
    # title = "High-risk metabolic clusters and incident type 2 diabetes",
    # subtitle = "Overall and stratified by microbial gene richness"
  ) +
  
  theme_linedraw(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 13, color = "black"),
    plot.title = element_text(face = "bold"),
    plot.margin = margin(10, 120, 10, 10)
  )

pdf('Outputs/Plot_OR_strat_gene_richness_binary.pdf', height = 3, width = 5)
print(plot_OR_strat_gene_richness)
dev.off()