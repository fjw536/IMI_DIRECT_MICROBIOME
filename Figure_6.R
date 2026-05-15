# Figure 6a
# Validation in the prediabetic intervention study
inter_metab <- read.delim('Inputs/Prediabetes_intervention_validation/Metabolites.txt',
                          header = TRUE, check.names = F)

library(tidyverse)
library(ggpubr)
library(dplyr)

dat <- inter_metab

# Force overwrite without asking
options(ggsave.overwrite = TRUE)

# Identify metabolite columns (column 4 onward)
metab_cols <- colnames(dat)[4:ncol(dat)]

# Ensure Time is ordered
dat$`TimePoint` <- factor(dat$`Time Point`,
                           levels = c("Pre-intervention", "Post-intervention"))

# Output folder
out_dir <- "Outputs/paired_plots_intervention_validation"
if (!dir.exists(out_dir)) dir.create(out_dir)

# -----------------------------
# Create empty dataframe to store paired p values
# -----------------------------
all_pvals <- data.frame(
  Metabolite = character(),
  Diet = character(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# -----------------------------
# Loop for plotting and p-value extraction
# -----------------------------
for (metab in metab_cols) {
  
  df_plot <- dat %>%
    dplyr::select(Diet, Participant = `Participant ID`, Time = `Time Point`, value = all_of(metab))
  
  # Safe filename
  safe_name <- metab %>%
    str_replace_all("[^A-Za-z0-9_\\-]", "_") %>% 
    str_replace_all("_+", "_")
  
  # Paired test by Diet
  stat_df <- df_plot %>%
    group_by(Diet) %>%
    summarise(
      p = tryCatch(
        wilcox.test(value[Time == "Pre-intervention"],
                    value[Time == "Post-intervention"],
                    paired = TRUE)$p.value,
        error = function(e) NA
      ),
      .groups = "drop"
    ) %>%
    mutate(label = paste0("Paired Wilcox p = ", signif(p, 3)))
  
  # -------- Append p values to all_pvals --------
  all_pvals <- rbind(
    all_pvals,
    data.frame(
      Metabolite = metab,
      Diet = stat_df$Diet,
      p_value = stat_df$p,
      stringsAsFactors = FALSE
    )
  )
  
  # Plot
  p <- ggplot(df_plot, aes(Time, value, group = Participant)) +
    geom_line(alpha = 0.4) +
    geom_point(size = 2) +
    geom_boxplot(aes(group = Time), width = 0.4, alpha = 0.4, outlier.shape = NA) +
    geom_jitter(aes(color = Time), width = 0.15, size = 1.8, alpha = 0.7) +
    facet_wrap(~Diet, scales = "free_y") +
    geom_text(
      data = stat_df,
      inherit.aes = FALSE,
      aes(x = 1.5, y = Inf, label = label),
      vjust = 1.2,
      size = 4
    ) +
    labs(
      x = "Time",
      y = metab,
      title = metab
    ) +
    theme_bw(base_size = 12) +
    theme(
      strip.background = element_rect(fill = "grey95"),
      legend.position = "none"
    )
  
  # Save pdf
  outfile <- paste0(out_dir, "/", safe_name, ".pdf")
  ggsave(outfile, p, width = 6, height = 5)
}


pvals_wide <- all_pvals %>%
  pivot_wider(
    names_from = Diet,
    values_from = p_value
  )

pvals_wide$MED_q <- p.adjust(pvals_wide$`MED diet`, 'BH')
pvals_wide$PPT_q <- p.adjust(pvals_wide$`PPT diet`, 'BH')
# -----------------------------
# Save final p value table
# -----------------------------
write.table(pvals_wide,
            file = paste0(out_dir, "/paired_test_p_values.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Figure 6b
bac_stat_IGR_T2D <- read.delim('Outputs/paired_plots_intervention_validation/statics_for_plots/Bacterial_species.txt', header = T, row.names = 1, check.names = F)

library(tidyverse)
library(ggplot2)
library(patchwork)


# 1 Convert rownames into a real column
bac_stat_IGR_T2D <- bac_stat_IGR_T2D %>%
  rownames_to_column(var = "Bacterial species")

# 2 Now you can select columns safely
bac_stat_IGR_T2D %>%
  select(`Bacterial species`, MeanAbsSHAP, MED_q, PPT_q) %>%
  head()


# Convert to long format for dot plot
bac_stat_IGR_T2D_long <- bac_stat_IGR_T2D %>%
  dplyr::select(`Bacterial species`, MeanAbsSHAP, MED_q, PPT_q) %>%   # use the real column
  pivot_longer(
    cols = c(MED_q, PPT_q),
    names_to = "Diet",
    values_to = "q_value"
  ) %>%
  mutate(
    neglog10 = -log10(q_value),
    show_dot = q_value < 0.1
  )


# reorder metabolites by MeanAbsSHAP
# stat_metab_IGR_T2D$Metabolite <- rownames(stat_metab_IGR_T2D)
bac_stat_IGR_T2D <- bac_stat_IGR_T2D %>% arrange(MeanAbsSHAP)
bac_stat_IGR_T2D <- bac_stat_IGR_T2D %>%
  mutate(`Bacterial species` = gsub("_", " ", `Bacterial species`))
bac_stat_IGR_T2D$`Bacterial species` <- factor(bac_stat_IGR_T2D$`Bacterial species`, levels = bac_stat_IGR_T2D$`Bacterial species`)

bac_stat_IGR_T2D_long <- bac_stat_IGR_T2D_long %>%
  mutate(`Bacterial species` = gsub("_", " ", `Bacterial species`))
bac_stat_IGR_T2D_long$`Bacterial species` <- factor(bac_stat_IGR_T2D_long$`Bacterial species`, levels = bac_stat_IGR_T2D$`Bacterial species`)

bac_stat_IGR_T2D_long <- bac_stat_IGR_T2D_long %>%
  mutate(Diet = case_when(
    Diet == "MED_q" ~ "MED diet",
    Diet == "PPT_q" ~ "PPT diet",
    TRUE ~ Diet
  ))

# ===== Panel 1: Bar plot of MeanAbsSHAP =====
p1_bac <- ggplot(bac_stat_IGR_T2D, aes(x = `Bacterial species`, y = MeanAbsSHAP)) +
  geom_col(fill = "#16609a") +
  coord_flip() +
  labs(
    x = "Bacterial species",
    y = "Absolute SHAP value"
  ) +
  theme_bw(base_size = 12) +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(color = "black"),
    axis.text.y = element_text(size = 8, face = "italic")
  )

# ===== Panel 2: Dot plot for p values =====
p2_bac <- ggplot(bac_stat_IGR_T2D_long, aes(x = `Bacterial species`, y = Diet)) +
  
  # background empty points (keeps structure)
  geom_point(color = NA) +

  # actual significant dots
  geom_point(
    data = bac_stat_IGR_T2D_long %>% filter(q_value < 0.1),
    aes(size = neglog10, color = Diet),
    alpha = 0.8
  ) +
  guides(color = "none") +
  scale_color_manual(values = c("MED diet" = "#e3bc6a", "PPT diet" = "#00bfc4"), guide = "none") +
  coord_flip() +
  scale_size_continuous(name = "-log10(qvalue)") +
  labs(
    x = NULL,
    y = "Diet"
  ) +
  theme_bw(base_size = 12) +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(color = "black"),
    axis.text.y = element_blank(),       # <<< remove metabolite labels
    legend.position = "right"
  )

# Combine
library(egg)
p_bac_validation <- ggarrange(p1_bac, p2_bac, nrow = 1, widths = c(1, 0.8))

pdf('Outputs/paired_plots_intervention_validation/IGR_T2D_bacterial_species_validation.pdf', height = 7, width = 7)
print(p_bac_validation)
dev.off()

# Figure 6c
stat_metab_IGR_T2D <- read.delim('Outputs/paired_plots_intervention_validation/statics_for_plots/Metabolome.txt', header = T, check.names = F, row.names = 1)

library(tidyverse)
library(ggplot2)
library(patchwork)


# 1 Convert rownames into a real column
stat_metab_IGR_T2D <- stat_metab_IGR_T2D %>%
  rownames_to_column(var = "Metabolite")

# 2 Now you can select columns safely
stat_metab_IGR_T2D %>%
  select(Metabolite, MeanAbsSHAP, `MED diet`, `PPT diet`) %>%
  head()


# Convert to long format for dot plot
stat_metab_IGR_T2D_long <- stat_metab_IGR_T2D %>%
  dplyr::select(Metabolite, MeanAbsSHAP, MED_q, PPT_q) %>%   # use the real column
  pivot_longer(
    cols = c('MED_q', 'PPT_q'),
    names_to = "Diet",
    values_to = "q_value"
  ) %>%
  mutate(
    neglog10 = -log10(q_value),
    show_dot = q_value < 0.1
  )

stat_metab_IGR_T2D_long <- stat_metab_IGR_T2D_long %>%
  mutate(Diet = case_when(
    Diet == "MED_q" ~ "MED diet",
    Diet == "PPT_q" ~ "PPT diet",
    TRUE ~ Diet
  ))
# reorder metabolites by MeanAbsSHAP
# stat_metab_IGR_T2D$Metabolite <- rownames(stat_metab_IGR_T2D)
stat_metab_IGR_T2D <- stat_metab_IGR_T2D %>% arrange(MeanAbsSHAP)
stat_metab_IGR_T2D$Metabolite <- factor(stat_metab_IGR_T2D$Metabolite, levels = stat_metab_IGR_T2D$Metabolite)

stat_metab_IGR_T2D_long$Metabolite <- factor(stat_metab_IGR_T2D_long$Metabolite, levels = stat_metab_IGR_T2D$Metabolite)

# ===== Panel 1: Bar plot of MeanAbsSHAP =====
p1 <- ggplot(stat_metab_IGR_T2D, aes(x = Metabolite, y = MeanAbsSHAP)) +
  geom_col(fill = "#16609a") +
  coord_flip() +
  labs(
    x = "Metabolite",
    y = "Absolute SHAP value"
  ) +
  theme_bw(base_size = 12) +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(color = "black"),
    axis.text.y = element_text(size = 8)
  )

# ===== Panel 2: Dot plot for p values =====
p2 <- ggplot(stat_metab_IGR_T2D_long, aes(x = Metabolite, y = Diet)) +
  
  # background empty points (keeps structure)
  geom_point(color = NA) +

  # actual significant dots
  geom_point(
    data = stat_metab_IGR_T2D_long %>% filter(q_value < 0.1),
    aes(size = neglog10, color = Diet),
    alpha = 0.8
  ) +
  guides(color = "none") +
  scale_color_manual(values = c("MED diet" = "#e3bc6a", "PPT diet" = "#00bfc4"), guide = "none") +
  coord_flip() +
  scale_size_continuous(name = "-log10(qvalue)") +
  labs(
    x = NULL,
    y = "Diet"
  ) +
  theme_bw(base_size = 12) +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(color = "black"),
    axis.text.y = element_blank(),       # <<< remove metabolite labels
    legend.position = "right"
  )

# Combine
library(egg)
p_metab_validation <- ggarrange(p1, p2, nrow = 1, widths = c(1, 0.7))

pdf('Outputs/paired_plots_intervention_validation/IGR_T2D_metabolites_validation.pdf', height = 10, width = 8)
print(p_metab_validation)
dev.off()

nrow(stat_metab_IGR_T2D)
stat_metab_IGR_T2D %>% filter(MED_q < 0.1 | PPT_q < 0.1) %>% nrow()

stat_metab_IGR_T2D %>% filter(MED_q < 0.1) %>% nrow()

stat_metab_IGR_T2D %>% filter(PPT_q < 0.1) %>% nrow()

stat_metab_IGR_T2D$MDE_qval <- p.adjust(stat_metab_IGR_T2D$`MED diet`, 'BH')
stat_metab_IGR_T2D$PPT_qval <- p.adjust(stat_metab_IGR_T2D$`PPT diet`, 'BH')