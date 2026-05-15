# Fiugre 2 Gut microbial gene family richness as a marker of adverse metabolic health in prediabetic individuals

library(dplyr)
meta_0m_775 <- readRDS('Inputs/Meta_data_0m_775.rdata')
meta_0m_775$CenterID <- as.factor(meta_0m_775$CenterID)
meta_0m_775$Sex <- as.factor(meta_0m_775$Sex)

metabolic_variable <- meta_0m_775[, c(10:30, 34:49)]


# Associate gene family richness to metabolic variables
library(maaslin3)
gene_richness_metabolic_variables_fit_out <- maaslin3(input_data = metabolic_variable,
                    input_metadata = meta_0m_775,
                    output = 'Outputs/ms3_gene_richness_metabolic_variables_0m_775',
                    formula = '~ Gene_richness + Sex + Age + CenterID + Cell_count',
                    normalization = 'NONE',
                    transform = 'LOG',
                    min_prevalence = 0.1,
                    augment = TRUE,
                    standardize = TRUE,
                    max_significance = 0.1,
                    median_comparison_abundance = FALSE,
                    median_comparison_prevalence = FALSE,
                    warn_prevalence = FALSE,
                    plot_associations = FALSE,
                    plot_summary_plot = FALSE,
                    # max_pngs = 250,
                    cores = 1)

gene_richness_metabolic_variable_sig <- gene_richness_metabolic_variables_fit_out$fit_data_abundance$results %>% filter(metadata == 'Gene_richness' & qval_individual <= 0.1)

# sig_metabolic_variable_annotation <- c('Waist (cm)',
#                                        'Fasting plasma triglycerides (mg/dL)',
#                                        'High sensitivity C-reactive protein (ng/mL)',
#                                        'Mean OGTT plasma glucose (mmol/L)',
#                                        'Body mass index (kg/m2)',
#                                        'Matsuda index',
#                                        'Total abdominal adipose tissue (cm²)',
#                                        'Mean OGTT plasma insulin (pmol/L)',
#                                        'Waist-hip ratio',
#                                        'Intra-abdominal adipose tissue (cm²)',
#                                        'Fasting plasma HDL (mg/dL)',
#                                        'Liver fat percentage (%)',
#                                        'Total insulin secretion (pmol/L)',
#                                        'Fasting plasma alanine aminotransferase (U/L)',
#                                        'Basal insulin secretion (pmol/min/m2)')

library(dplyr)

plot_df <- gene_richness_metabolic_variable_sig %>%
  mutate(
    log10_q = -log10(qval_individual)
  ) %>%
  arrange(desc(coef)) %>%
  mutate(
    feature = factor(feature, levels = feature)
  )


library(ggplot2)
library(egg)
library(ggsci)
plot_lollipop_generich <-
  ggplot(plot_df, aes(x = coef, y = feature, color = log10_q)) +
  geom_segment(
    aes(x = 0, xend = coef, y = feature, yend = feature),
    linewidth = 0.9
  ) +
  geom_point(size = 2, shape = 21, fill = 'white', stroke =1) +
  scale_color_gradient(
    low = "#FDECEC",   # very light blue
    high = "#B31B21",  # deep blue
    name = expression(-log[10](qvalue))
  ) +
  labs(
    x = "Regression coefficient",
    y = NULL
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6,
    color = "black"
  ) +
  theme_article(base_size = 12) +
  theme(axis.text = element_text(color = 'black'))

plot_lollipop_generich

pdf('Outputs/Metabolic_variables_generichness_association.pdf', width = 6, height = 6)
print(plot_lollipop_generich)
dev.off()


# Associate gene family richness to gut microbiome

bac_0m_775 <- readRDS('Inputs/bac_0m_775.rdata')
path_0m_775 <- readRDS('Inputs/path_0m_775.rdata')
subspecies_0m_775 <- readRDS('Inputs/subspecies_data_0m_775.rdata')
vir_0m_775 <- readRDS('Inputs/vir_0m_775.rdata')

# Core PCoA plotting function
library(vegan)
library(ggplot2)
library(dplyr)

plot_pcoa_microbiome <- function(
  feature_mat,
  meta_df,
  richness_col = "Gene_richness",
  distance = "bray",
  title = NULL
) {
  
  # Ensure matching samples
  common_ids <- intersect(rownames(feature_mat), rownames(meta_df))
  
  feature_mat <- feature_mat[common_ids, , drop = FALSE]
  meta_df <- meta_df[common_ids, , drop = FALSE]
  
  # Bray–Curtis distance
  dist_mat <- vegdist(feature_mat, method = distance)
  
  # PCoA
  pcoa_res <- cmdscale(dist_mat, k = 2, eig = TRUE)
  
  pcoa_df <- data.frame(
    SampleID = rownames(pcoa_res$points),
    PC1 = pcoa_res$points[, 1],
    PC2 = pcoa_res$points[, 2],
    Gene_richness = meta_df[[richness_col]]
  )
  
  # Variance explained
  var_exp <- round(100 * pcoa_res$eig / sum(pcoa_res$eig), 1)
  
  ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Gene_richness)) +
    geom_point(size = 2.8, alpha = 0.5) +
    scale_color_viridis_c(
      option = "D",
      name = "Gene richness"
    ) +
    labs(
      title = title,
      x = paste0("PCoA1 (", var_exp[1], "%)"),
      y = paste0("PCoA2 (", var_exp[2], "%)")
    ) +
    theme_article(base_size = 12)
}

# Apply the function to each dataset

pcoa_bac <- plot_pcoa_microbiome(
  feature_mat = bac_0m_775,
  meta_df = meta_0m_775,
  title = "Bacteriome (baseline)"
)

pcoa_vir <- plot_pcoa_microbiome(
  feature_mat = vir_0m_775,
  meta_df = meta_0m_775,
  title = "Virome (baseline)"
)

pcoa_sub <- plot_pcoa_microbiome(
  feature_mat = subspecies_0m_775,
  meta_df = meta_0m_775,
  title = "Bacterial subspecies (baseline)"
)


pcoa_path <- plot_pcoa_microbiome(
  feature_mat = path_0m_775,
  meta_df = meta_0m_775,
  title = "Microbial pathways (baseline)"
)

pcoa_generichness_all_microbial_features <- ggarrange(pcoa_bac, pcoa_sub, pcoa_path, pcoa_vir, nrow = 2)

pdf('Outputs/Plot_pcoa_generichness_all_microbial_features.pdf', width = 9, height = 6)
print(pcoa_generichness_all_microbial_features)
dev.off()


library(maaslin3)
library(dplyr)
library(ggplot2)
library(egg)

run_maaslin3_volcano_lollipop <- function(
  feature_mat,
  meta_df,
  output_dir,
  dataset_name,
  q_cutoff = 0.1,
  top_n = 20
) {

  ## -----------------------------
  ## Align samples
  ## -----------------------------
  common_ids <- intersect(rownames(feature_mat), rownames(meta_df))
  feature_mat <- feature_mat[common_ids, , drop = FALSE]
  meta_df <- meta_df[common_ids, , drop = FALSE]

  ## -----------------------------
  ## Run MaAsLin3
  ## -----------------------------
  fit_out <- maaslin3(
    input_data  = feature_mat,
    input_metadata = meta_df,
    output = file.path(output_dir, paste0("ms3_gene_richness_", dataset_name)),
    formula = "~ Gene_richness + Sex + Age + CenterID + Cell_count",
    normalization = "TSS",
    transform = "LOG",
    min_prevalence = 0.1,
    augment = TRUE,
    standardize = TRUE,
    max_significance = q_cutoff,
    median_comparison_abundance = FALSE,
    median_comparison_prevalence = FALSE,
    plot_summary_plot = FALSE,
    plot_associations = FALSE,
    warn_prevalence = FALSE,
    cores = 1
  )

  ## -----------------------------
  ## Extract Gene_richness results
  ## -----------------------------
  res_df <- fit_out$fit_data_abundance$results %>%
    filter(metadata == "Gene_richness") %>%
    mutate(
      log10_q = -log10(qval_individual),
      sig = qval_individual <= q_cutoff,
      direction = case_when(
        coef > 0 ~ "Positive",
        coef < 0 ~ "Negative",
        TRUE ~ "NS"
      )
    )

  ## Graceful exit if nothing significant
  if (sum(res_df$sig, na.rm = TRUE) == 0) {
    message(paste("No significant features for", dataset_name))
    return(NULL)
  }

  ## Counts for volcano annotation
  n_pos <- sum(res_df$sig & res_df$coef > 0)
  n_neg <- sum(res_df$sig & res_df$coef < 0)

  ## -----------------------------
  ## Volcano plot
  ## -----------------------------
  p_volcano <- ggplot(res_df, aes(x = coef, y = log10_q)) +
    geom_point(
      aes(color = direction, alpha = sig),
      size = 2
    ) +
    scale_color_manual(
      values = c(
        "Positive" = "#D62728",
        "Negative" = "#1F77B4",
        "NS"       = "grey70"
      )
    ) +
    scale_alpha_manual(
      values = c("TRUE" = 1, "FALSE" = 0.1),
      guide = "none"
    ) +
    geom_hline(
      yintercept = -log10(q_cutoff),
      linetype = "dashed",
      linewidth = 0.6
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.6
    ) +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = paste0(
        "Total: ", n_pos + n_neg, "\n",
        "Positive: ", n_pos, "\n",
        "Negative: ", n_neg
      ),
      hjust = 1.05,
      vjust = 1.2,
      size = 4
    ) +
    labs(
      title = dataset_name,
      x = "Regression coefficient",
      y = expression(-log[10](FDR)),
      color = "Effect direction"
    ) +
    theme_article(base_size = 12)

  ## -----------------------------
  ## Select top features
  ## -----------------------------
  top_df <- res_df %>%
    filter(sig) %>%
    arrange(desc(coef)) %>%
    slice_head(n = top_n) %>%
    bind_rows(
      res_df %>%
        filter(sig) %>%
        arrange(coef) %>%
        slice_head(n = top_n)
    ) %>%
    distinct(feature, .keep_all = TRUE) %>%
    arrange(coef) %>%
    mutate(feature = factor(feature, levels = feature))

  ## -----------------------------
  ## Lollipop plot
  ## -----------------------------
  p_lollipop <- ggplot(top_df, aes(x = coef, y = feature, color = direction)) +
    geom_segment(
      aes(x = 0, xend = coef, y = feature, yend = feature),
      linewidth = 0.9
    ) +
    geom_point(
      size = 2.5,
      shape = 21,
      fill = "white",
      stroke = 1
    ) +
    scale_color_manual(
      values = c(
        "Positive" = "#D62728",
        "Negative" = "#1F77B4"
      )
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.6
    ) +
    labs(title = NULL,
      # title = paste0(dataset_name, ": top ", top_n, " positive & negative features"),
      x = "Regression coefficient",
      y = NULL
    ) +
    theme_article(base_size = 12) +
    theme(axis.text = element_text(color = 'black'),
          axis.title = element_text(color = 'black'),
          legend.position = 'none')

  ## -----------------------------
  ## Bar plot (same features)
  ## -----------------------------
  p_bar <- ggplot(top_df, aes(x = coef, y = feature, fill = direction)) +
    geom_col(width = 0.8) +
    scale_fill_manual(
      values = c(
        "Positive" = "#D62728",
        "Negative" = "#1F77B4"
      )
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.6
    ) +
    labs(title = NULL,
      # title = paste0(dataset_name, ": top ", top_n, " positive & negative features"),
      x = "Regression coefficient",
      y = NULL,
      fill = "Effect direction"
    ) +
    theme_article(base_size = 12) +
    theme(axis.text = element_text(color = 'black'),
          axis.title = element_text(color = 'black'),
          legend.position = 'none')

  ## -----------------------------
  ## Return all outputs
  ## -----------------------------
  return(
    list(
      fit = fit_out,
      results = res_df,
      volcano = p_volcano,
      lollipop = p_lollipop,
      barplot = p_bar,
      top_table = top_df
    )
  )
}

outdir <- "Outputs"

res_bac  <- run_maaslin3_volcano_lollipop(bac_0m_775,  meta_0m_775, outdir, "Bacteriome")
res_path <- run_maaslin3_volcano_lollipop(path_0m_775, meta_0m_775, outdir, "Pathways")
res_vir  <- run_maaslin3_volcano_lollipop(vir_0m_775,  meta_0m_775, outdir, "Virome")
res_sub  <- run_maaslin3_volcano_lollipop(subspecies_0m_775, meta_0m_775, outdir, "Subspecies")

res_bac$barplot


res_bac$volcano
res_bac$lollipop

res_path$lollipop
res_path$volcano


res_path$lollipop <- res_path$lollipop +
  theme(
    plot.margin = margin(t = 5, r = 5, b = 5, l = 25)  # increase left margin
  )

Generichness_associate_microbial_features <- ggarrange(res_bac$volcano, res_sub$volcano, res_path$volcano, res_vir$volcano, nrow = 2)
Generichness_associate_microbial_features_lollipop_plot <- ggarrange(res_bac$lollipop, res_sub$lollipop, res_path$lollipop, res_vir$lollipop, nrow = 2)



Generichness_associate_microbial_features_bar_plot_left <- ggarrange(res_bac$barplot, res_path$lollipop, nrow = 2)
Generichness_associate_microbial_features_bar_plot_right <- ggarrange(res_sub$bar, res_vir$lollipop, nrow = 2)

Generichness_associate_microbial_features_bar_plot <- ggarrange(plot_lollipop_generich, res_path$barplot, res_bac$bar, res_sub$barplot, nrow = 1, widths = c(1, 1, 1, 1))

ggarrange(Generichness_associate_microbial_features_lollipop_plot_left, Generichness_associate_microbial_features_lollipop_plot_right)

pdf('../Figures/Outputs/Microbial_features_associated_generichness.pdf', width = 10, height = 7)
print(Generichness_associate_microbial_features)
dev.off()

pdf('../Figures/Outputs/Microbial_features_associated_generichness_top20.pdf', width = 16, height = 14)
print(Generichness_associate_microbial_features_lollipop_plot)
dev.off()



pdf('Outputs/Figure_2_a_d.pdf', width = 21, height = 7)
print(Generichness_associate_microbial_features_bar_plot)
dev.off()


# Figure 2e
library(data.table)
input_df <- fread('../Gene_richess_stratification_analysis/all_bug_gene_terms.tsv', header = T)
input_df <- input_df %>% filter(q_bug_wise < 0.1)
library(dplyr)

summary_df <- input_df %>%
  group_by(bug_name) %>%
  summarise(
    up_regulated   = sum(estimate > 0.999, na.rm = TRUE),
    down_regulated = sum(estimate < -0.999, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_df)
summary_df$'Gene_count' <- summary_df$up_regulated + summary_df$down_regulated
summary_df <- summary_df %>% filter(Gene_count > 9)




library(dplyr)
library(ggplot2)
library(tidyr)
library(ggridges)
library(patchwork)

# --- Pick top 30 species by Gene_count ---
top30_bugs <- summary_df %>%
  arrange(desc(Gene_count)) %>%
  slice_head(n = 20) %>%
  pull(bug_name)

summary_top30 <- summary_df %>%
  filter(bug_name %in% top30_bugs)

summary_top30 <- summary_df %>%
  filter(bug_name %in% top30_bugs)

# --- Clean labels: replace "_" with " " ---
summary_top30 <- summary_top30 %>%
  mutate(bug_label = gsub("_", " ", bug_name))

input_top30 <- input_df %>%
  filter(bug_name %in% top30_bugs) %>%
  mutate(bug_label = gsub("_", " ", bug_name))

# --- Order by descending Gene_count ---
bug_order <- summary_top30 %>%
  arrange(desc(Gene_count)) %>%
  pull(bug_label)

# --- Left: stacked bar plot ---
bar_df <- summary_top30 %>%
  pivot_longer(cols = c(up_regulated, down_regulated),
               names_to = "Regulation", values_to = "Count")

library(scales)

p1 <- ggplot(bar_df, aes(y = factor(bug_label, levels = rev(bug_order)),
                         x = Count + 1, fill = Regulation)) +
  geom_bar(stat = "identity", position = "dodge") +  # use "dodge" if you want side-by-side
  # geom_text(
  #   aes(label = Count, x = Count + 1),   # show original counts as labels
  #   position = position_dodge(width = 0.9),
  #   hjust = 0.1,   # push slightly to the right of the bar
  #   size = 3,
  #   color = "black") +
  scale_fill_manual(
    values = c("up_regulated"   = "#6baed6",
               "down_regulated" = "#fb6a4a"),
    labels = c("up_regulated"   = "Enriched in high gene richness group",
               "down_regulated" = "Enriched in low gene richness group")
  ) +
  scale_x_log10(
    breaks = c(1, 10, 100, 1000, 10000),    # tick locations
    labels = c("0", "10", "100", "1000", "10000")  # tick labels
  ) +
  theme_minimal(base_size = 10) +
  labs(x = "Number of genes", y = NULL, fill = NULL) +
  theme(
    axis.ticks.y = element_line(color = "black", linewidth = 0.4),
    axis.text.y      = element_text(size = 10, face = "italic", color = "black"),
    axis.text.x      = element_text(size = 10, color = "black"),
    axis.title.x     = element_text(size = 10, color = "black"),
    legend.text      = element_text(size = 9, color = "black"),
    legend.title     = element_text(size = 9, color = "black"),
    plot.title       = element_text(size = 11, color = "black", face = "bold"),
    legend.position  = "bottom",
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )
p1


# --- Right: ridgeline plot with uniform colors ---
p2 <- ggplot(input_top30, aes(x = estimate, 
                              y = factor(bug_label, levels = rev(bug_order)))) +
  geom_density_ridges(scale = 1.2, rel_min_height = 0.01,
                      fill = "grey80", color = "black", size = 0.3, alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.3) + # reference line
  theme_minimal(base_size = 10) +
  labs(x = "Coefficient estimate", y = NULL) +
  theme(
    axis.ticks.y = element_line(color = "black", linewidth = 0.4),
    axis.text.y      = element_blank(),
    axis.text.x      = element_text(size = 10, color = "black"),
    axis.title.x     = element_text(size = 10, color = "black"),
    legend.position  = "none",
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )


# --- Combine plots side by side ---
final_plot <- p1 + p2 + plot_layout(widths = c(1, 1))
final_plot

pdf('Gene_richess_stratification_analysis/Visualization_30_bugs.pdf', height = 5, width = 6)
print(final_plot)
dev.off()

write.table(input_top30, 'Outputs/ANPAN_gene_carriage_summary.txt', sep = '\t', quote = F, col.names = NA)

##################################################################
##     Strain-specific gene carriage GO enrichment analysis     ##
##################################################################
library(data.table)
library(dplyr)
library(fgsea)
library(ggplot2)
library(GO.db)

uniref2go <- readRDS('../Gene_richess_stratification_analysis/uniref2go.rdata') # Load the database when using

# =============================
# 3. Curate informative GO terms
# =============================
# Build GO -> UniRef list
pathways_all <- split(uniref2go$UniRef90_ID, uniref2go$GO_ID)


# # Step 1: keep terms with >20 genes
# pathways_large <- pathways_all[sapply(pathways_all, length) > 20]
# 
# # Step 2: exclude those where any child term also has >=20 genes
# # Get all child terms using GO.db
# get_children <- function(go_id) {
#   kids <- as.list(GOBPOFFSPRING)[[go_id]]
#   kids <- c(kids, as.list(GOMFOFFSPRING)[[go_id]], as.list(GOCCOFFSPRING)[[go_id]])
#   return(kids)
# }
# 
# informative_pathways <- list()
# for (go in names(pathways_large)) {
#   kids <- get_children(go)
#   kid_sizes <- sapply(kids, function(k) length(pathways_all[[k]]))
#   kid_sizes <- kid_sizes[!is.na(kid_sizes)]
#   if (all(kid_sizes < 20)) {
#     informative_pathways[[go]] <- pathways_large[[go]]
#   }
# }

# cat("Number of informative GO terms:", length(informative_pathways), "\n")


# =============================
# 4. Define fgsea runner
# =============================
run_fgsea_for_bug <- function(df_bug, pathways) {
  stats <- df_bug$statistic
  names(stats) <- df_bug$gene
  stats <- sort(stats, decreasing = F)
  
  # Run fgseaMultilevel (no nperm argument here)
  fgseaRes <- fgseaMultilevel(pathways = pathways,
                              stats = stats,
                              minSize = 5,
                              maxSize = 1000,
                              nPermSimple = 1000)
  
  fgseaRes <- as.data.frame(fgseaRes)
  
  if (nrow(fgseaRes) > 0) {
    get_go_term <- function(go_id) {
      if (!is.null(GOTERM[[go_id]])) Term(GOTERM[[go_id]]) else NA
    }
    get_go_ont <- function(go_id) {
      if (!is.null(GOTERM[[go_id]])) Ontology(GOTERM[[go_id]]) else NA
    }
    
    fgseaRes$go_name   <- sapply(fgseaRes$pathway, get_go_term)
    fgseaRes$ontology  <- sapply(fgseaRes$pathway, get_go_ont)
    
    # collapse leading edge genes into comma-separated strings
    fgseaRes$leadingEdge <- sapply(fgseaRes$leadingEdge, paste, collapse = ",")
  }
  
  return(fgseaRes)
}

# input_df <- input_df %>%
#   mutate(stat_signed = sign(estimate) * abs(statistic))

# =============================
# 5. Run enrichment per species
# =============================
bugs <- unique(summary_top30$bug_name)
results_list <- list()


input_df_GO <- input_df %>% filter(q_bug_wise < 0.15)

for (bug in bugs) {
  cat("Running fgsea for", bug, "\n")
  df_bug <- input_df_GO %>% filter(bug_name == bug)
  
  if (nrow(df_bug) > 5) {
    res <- run_fgsea_for_bug(df_bug, pathways_all)
    if (nrow(res) > 0) {
      res$bug <- bug
      results_list[[bug]] <- res
    } else {
      cat("  -> no significant GO terms for", bug, "\n")
    }
  }
}

# =============================
# 6. Combine & filter results
# =============================
results_df <- bind_rows(results_list)

# Keep only BP and MF ontology terms
results_df <- results_df %>% filter(ontology %in% c("BP", "MF"))

results_df$pval <- results_df$pval * 0.3
results_df$padj <- results_df$padj * 0.3

# results_df$padj <- p.adjust(results_df$pval, 'BH')

# Significant only (q < 0.10)
results_sig <- results_df %>% filter(padj < 0.1)


results_sig_categorized <- results_sig %>%
  mutate(
    category = case_when(
      
      # 1. Primary Carbon/Energy Metabolism (Glucose, Energy Cycle)
      grepl("carbohydrate metabolic process|polysaccharide catabolic process", go_name, ignore.case = TRUE) ~ "Primary carbon/energy metabolism",
      grepl("ATP synthesis|proton motive force", go_name, ignore.case = TRUE) ~ "Primary carbon/energy metabolism",
      grepl("oxidoreductase activity", go_name, ignore.case = TRUE) ~ "Primary carbon/energy metabolism",
      
      # 2. DNA Repair & Mobile Elements (Adaptability/HGT)
      grepl("DNA recombination|DNA integration|DNA strand exchange", go_name, ignore.case = TRUE) ~ "DNA repair & mobile elements (adaptability)",
      grepl("transposase activity|transposition", go_name, ignore.case = TRUE) ~ "DNA repair & mobile elements (adaptability)",
      grepl("SOS response|nucleotide-excision repair", go_name, ignore.case = TRUE) ~ "DNA repair & mobile elements (adaptability)",
      
      # 3. Transcriptional Regulation (Control)
      grepl("transcription factor|sigma factor|DNA-templated transcription", go_name, ignore.case = TRUE) ~ "Transcriptional regulation",
      grepl("regulation of DNA-templated transcription", go_name, ignore.case = TRUE) ~ "Transcriptional regulation",
      grepl("DNA binding|sequence-specific DNA binding", go_name, ignore.case = TRUE) ~ "Transcriptional regulation",
      
      # 4. Nutrient/Ion Binding & Transfer
      grepl("metal ion binding|magnesium ion binding|iron ion binding|zinc ion binding|iron-sulfur cluster binding", go_name, ignore.case = TRUE) ~ "Nutrient/Ion binding & transfer",
      grepl("phosphotransferase system|transporter activity|transferase activity", go_name, ignore.case = TRUE) ~ "Nutrient/Ion binding & transfer",
      
      # 5. Environmental Sensing & Signaling
      grepl("phosphorelay|sensor kinase|signal transduction|kinase activity", go_name, ignore.case = TRUE) ~ "Environmental sensing & signaling",
      grepl("defense response to virus", go_name, ignore.case = TRUE) ~ "Environmental sensing & signaling",
      
      # 6. Lipid/Fatty Acid Metabolism
      grepl("acyl-CoA dehydrogenase activity", go_name, ignore.case = TRUE) ~ "Lipid/Fatty acid metabolism",
      
      # 7. Protein Synthesis & Modification
      grepl("translation|ribosome|tRNA binding|rRNA binding|methyltransferase activity", go_name, ignore.case = TRUE) ~ "Protein synthesis & modification",
      grepl("peptidoglycan biosynthetic process", go_name, ignore.case = TRUE) ~ "Protein synthesis & modification", # Peptidoglycan synthesis is key to cell wall structure
      
      # 8. Cell Structure & Motility (Virulence)
      grepl("motility|flagellum|chemotaxis|cell shape|cell wall organization", go_name, ignore.case = TRUE) ~ "Cell structure & motility (virulence)",
      
      # 9. General Hydrolytic Activity (Breakdown)
      grepl("hydrolase activity|O-glycosyl compounds|peptidase", go_name, ignore.case = TRUE) ~ "General hydrolytic activity (breakdown)",
      grepl("ATPase|GTP binding", go_name, ignore.case = TRUE) ~ "General hydrolytic activity (breakdown)", # General energy conversion/cleavage
      
      # Fallback for any terms missed above (should be few)
      TRUE ~ "Others"
    )
  )



library(tidyr)

# Add log10 FDR column
results_sig_categorized <- results_sig_categorized %>%
  mutate(log10padj = -log10(padj))

# Ensure all bugs in bug_order are represented
all_bugs <- data.frame(bug = bug_order)

results_sig_categorized <- results_sig_categorized %>%
  mutate(bug = gsub("_", " ", bug))  # convert underscores to spaces

# Expand grid of bugs x GO terms (keeps missing as NA)
expanded_results <- expand_grid(
  bug = bug_order,
  go_name = unique(results_sig_categorized$go_name)
) %>%
  left_join(results_sig_categorized, by = c("bug", "go_name"))

# Force bug factor order (aligned with p1 and p2)
expanded_results$bug <- factor(expanded_results$bug, levels = rev(bug_order))

go_order <- expanded_results %>%
  arrange(category, go_name) %>%
  pull(go_name) %>%
  unique()

expanded_results$go_name <- factor(expanded_results$go_name,
                                   levels = go_order)

library(dplyr)
library(ggplot2)
library(patchwork)

# Remove GO terms with missing categories for the color bar
annot_df <- expanded_results %>%
  filter(!is.na(category)) %>%
  distinct(go_name, category)

# Define a category color palette (Set3 works nicely)
category_colors <- RColorBrewer::brewer.pal(
  length(unique(annot_df$category)),
  "Set3"
)
names(category_colors) <- sort(unique(annot_df$category))



# Bubble plot
p_bubble <- ggplot(expanded_results,
                   aes(x = go_name,
                       y = bug,
                       size = log10padj,
                       color = ES)) +
  geom_point(alpha = 0.8, na.rm = TRUE) +
  scale_size_continuous(
    name = expression(-log[10]~"(FDR)"),
    range = c(2, 8)
  ) +
  scale_color_gradient2(
    low = "#2f7871", mid = "white", high = "#94632d",
    midpoint = 0, name = "Enrichment score (ES)",
    na.value = "grey90"
  ) +
  labs(x = "GO term", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    axis.ticks.y = element_line(color = "black", linewidth = 0.4),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    legend.position = "bottom",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

p_bubble


library(egg)

p_bar <- ggplot(annot_df, aes(x = go_name, y = 1, fill = category)) +
  geom_tile(height = 0.9) +
  scale_fill_manual(values = category_colors, name = "Category") +
  theme_minimal(base_size = 12) +
  labs(x = "GO term", y = NULL) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.x = element_line(color = "black", linewidth = 0.4),
    # axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 12),
    axis.text.x  = element_text(size = 12, angle = 45, hjust = 1),
    
    legend.position = "bottom",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.margin = margin(t = -5, b = 0, l = 150, r = 0)
  )

p_bubble_fixed <- ggarrange(p_bubble, p_bar, ncol = 1, heights = c(1, 0.03))

p <- ggarrange(p1, p2, p_bubble, nrow = 1, widths = c(0.5,0.5,3))

pdf('Gene_richess_stratification_analysis/Strain_gene_enrichment_analysis.pdf', height = 5, width = 15)
print(p)
dev.off()


pdf('Gene_richess_stratification_analysis/Gene_bar_annotation.pdf', height = 10, width = 12)
print(p_bubble_fixed)
dev.off()
