###############################################################################
# Script: statistics_box_violin_plots.R
# Author: Koray Durak
# Contact: koray.durak@outlook.de
# Last updated: 2026-03-02
#
# Purpose:
#   Create reproducible violin + box plots from an online medical dataset,
#   with continuous biomarkers across three diagnosis groups.
###############################################################################

## ===============================
## 1) Dependencies
## ===============================
required_packages <- c("MedDataSets", "dplyr", "ggplot2", "tidyr", "forcats", "patchwork")

install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
}

install_if_missing(required_packages)

library(MedDataSets)
library(dplyr)
library(ggplot2)
library(tidyr)
library(forcats)
library(patchwork)


## ===============================
## 2) Load and prepare medical dataset
## ===============================
# Chosen dataset: Cushings_df
# - Diagnosis/groups: Type
# - Continuous biomarkers: Tetrahydrocortisone, Pregnanetriol

data("Cushings_df", package = "MedDataSets")

analysis_data <- Cushings_df %>%
  mutate(
    diagnosis_raw = as.character(Type)
  ) %>%
  filter(
    !is.na(diagnosis_raw),
    !is.na(Tetrahydrocortisone),
    !is.na(Pregnanetriol)
  )

# Keep exactly three diagnosis groups (most frequent), mimicking a 3-group study setup
diagnosis_top3 <- analysis_data %>%
  count(diagnosis_raw, sort = TRUE) %>%
  slice_head(n = 3) %>%
  pull(diagnosis_raw)

if (length(diagnosis_top3) < 3) {
  stop("Dataset does not contain at least three diagnosis groups after filtering.")
}

analysis_data <- analysis_data %>%
  filter(diagnosis_raw %in% diagnosis_top3) %>%
  mutate(
    diagnosis_group = factor(diagnosis_raw),
    diagnosis_group = fct_infreq(diagnosis_group)
  )


## ===============================
## 3) Long-format biomarker table
## ===============================
plot_data <- analysis_data %>%
  select(diagnosis_group, Tetrahydrocortisone, Pregnanetriol) %>%
  pivot_longer(
    cols = c(Tetrahydrocortisone, Pregnanetriol),
    names_to = "biomarker",
    values_to = "value"
  ) %>%
  mutate(
    biomarker_label = case_when(
      biomarker == "Tetrahydrocortisone" ~ "Tetrahydrocortisone (continuous biomarker)",
      biomarker == "Pregnanetriol" ~ "Pregnanetriol (continuous biomarker)",
      TRUE ~ biomarker
    )
  )


## ===============================
## 4) Violin + box plots
## ===============================
diagnosis_levels <- levels(analysis_data$diagnosis_group)
group_colors <- setNames(
  c("#2ca25f", "firebrick3", "#8800ff")[seq_along(diagnosis_levels)],
  diagnosis_levels
)

plot_template <- function(data, y_var, y_label) {
  ggplot(
    data,
    aes(x = .data$diagnosis_group, y = .data[[y_var]], color = .data$diagnosis_group)
  ) +
    geom_violin(
      aes(fill = .data$diagnosis_group),
      trim = FALSE,
      alpha = 0.2,
      linewidth = 0.8
    ) +
    geom_boxplot(
      width = 0.35,
      outlier.shape = NA,
      linewidth = 0.8,
      alpha = 0.85
    ) +
    geom_jitter(
      width = 0.08,
      size = 2.3,
      alpha = 0.75
    ) +
    scale_color_manual(values = group_colors) +
    scale_fill_manual(values = group_colors) +
    labs(
      x = "Diagnosis group",
      y = y_label
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      axis.title.x = element_text(size = 16, face = "bold", margin = margin(t = 10)),
      axis.title.y = element_text(size = 16, face = "bold", margin = margin(r = 10)),
      axis.text.x = element_text(size = 16, face = "bold", color = "black"),
      axis.text.y = element_text(size = 16, color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.6),
      plot.margin = margin(15, 20, 15, 15)
    )
}

p_tetra <- plot_template(
  analysis_data,
  y_var = "Tetrahydrocortisone",
  y_label = "Tetrahydrocortisone"
)

p_preg <- plot_template(
  analysis_data,
  y_var = "Pregnanetriol",
  y_label = "Pregnanetriol"
)

p_biomarkers <- p_tetra + p_preg +
  plot_layout(ncol = 2, widths = c(1, 1)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(size = 30, face = "bold"))
  )

p_biomarkers


## ===============================
## 5) Export outputs
## ===============================
output_plot <- "box_violin_biomarkers_cushings_top3.png"

ggsave(
  filename = output_plot,
  plot = p_biomarkers,
  width = 14,
  height = 7,
  dpi = 400
)

summary_table <- plot_data %>%
  group_by(biomarker_label, diagnosis_group) %>%
  summarise(
    n = n(),
    median = median(value, na.rm = TRUE),
    q25 = quantile(value, 0.25, na.rm = TRUE),
    q75 = quantile(value, 0.75, na.rm = TRUE),
    min = min(value, na.rm = TRUE),
    max = max(value, na.rm = TRUE),
    .groups = "drop"
  )

utils::write.csv(summary_table, "box_violin_biomarkers_cushings_top3_summary.csv", row.names = FALSE)
utils::write.csv(analysis_data, "box_violin_biomarkers_cushings_top3_dataset.csv", row.names = FALSE)

cat("Violin + box plot workflow complete.\n")
cat("Diagnosis groups used:\n")
print(levels(analysis_data$diagnosis_group))
cat("Outputs:\n")
cat("- box_violin_biomarkers_cushings_top3.png\n")
cat("- box_violin_biomarkers_cushings_top3_summary.csv\n")
cat("- box_violin_biomarkers_cushings_top3_dataset.csv\n")
