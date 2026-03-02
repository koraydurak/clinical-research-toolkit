###############################################################################
# Script: statistics_relative_true_positive_fraction.R
# Author: Koray Durak
# Contact: koray.durak@outlook.de
# Last updated: 2026-03-02
#
# Purpose:
#   Reproducible relative true positive fraction (rTPF) analysis with one
#   combined forest-style plot.
#   For studies comparing two diagnostic methods, rTPF is a suitable method in cases
#   where you do not want to miss a diagnosis
# 
# Data source:
#   MedDataSets::Pima_tr2_df
#
# Concept:
#   - Gold standard outcome: diabetes status (type)
#   - Method A (intervention-like): glucose-based screening (glu >= 140)
#   - Method B (control-like): BMI-based screening (bmi >= 30)
#   - rTPF (BMI vs glucose) = Sensitivity(BMI method) / Sensitivity(Glucose method)
#   - Non-inferiority margin for rTPF = 0.8
###############################################################################

## ===============================
## 1) Dependencies
## ===============================
required_packages <- c("MedDataSets", "dplyr", "ggplot2")

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


## ===============================
## 2) Load and prepare dataset
## ===============================
data("Pima_tr2_df", package = "MedDataSets")

as_binary_outcome <- function(x) {
  if (is.numeric(x) || is.integer(x)) {
    return(ifelse(x > 0, 1L, 0L))
  }

  x_chr <- tolower(trimws(as.character(x)))
  positive_labels <- c("1", "yes", "true", "diabetes", "positive")
  ifelse(x_chr %in% positive_labels, 1L, 0L)
}

analysis_data <- Pima_tr2_df %>%
  transmute(
    gold_standard = as_binary_outcome(type),
    glucose = as.numeric(glu),
    bmi = as.numeric(bmi)
  ) %>%
  filter(
    !is.na(gold_standard),
    !is.na(glucose),
    !is.na(bmi)
  ) %>%
  mutate(
    method_glucose = ifelse(glucose >= 140, 1L, 0L),
    method_bmi = ifelse(bmi >= 30, 1L, 0L)
  )

if (sum(analysis_data$gold_standard == 1L) == 0) {
  stop("No positive cases in gold standard after filtering.")
}


## ===============================
## 3) Core metrics
## ===============================
compute_metrics <- function(df) {
  diseased <- df$gold_standard == 1L
  n_diseased <- sum(diseased)

  tp_glucose <- sum(df$method_glucose == 1L & diseased)
  tp_bmi <- sum(df$method_bmi == 1L & diseased)

  sens_glucose <- tp_glucose / n_diseased
  sens_bmi <- tp_bmi / n_diseased
  rTPF <- sens_bmi / sens_glucose

  list(
    n_diseased = n_diseased,
    tp_glucose = tp_glucose,
    tp_bmi = tp_bmi,
    sens_glucose = sens_glucose,
    sens_bmi = sens_bmi,
    rTPF = rTPF
  )
}

observed <- compute_metrics(analysis_data)

if (observed$tp_glucose == 0) {
  stop("Glucose method has zero true positives, rTPF (BMI/glucose) is not estimable.")
}


## ===============================
## 4) Formula-based confidence intervals
## ===============================
# Sensitivity CIs: Wald approximation
# rTPF CI: log-Wald approximation for ratio of sensitivities
# Note: this is a simple formula-based approach without bootstrap.

alpha <- 0.05
z_crit <- qnorm(1 - alpha / 2)

ci_wald_prop <- function(successes, n, z = z_crit) {
  p <- successes / n
  se <- sqrt((p * (1 - p)) / n)
  lower <- max(0, p - z * se)
  upper <- min(1, p + z * se)
  c(lower = lower, upper = upper)
}

ci_sens_glucose <- ci_wald_prop(observed$tp_glucose, observed$n_diseased)
ci_sens_bmi <- ci_wald_prop(observed$tp_bmi, observed$n_diseased)

p1 <- observed$sens_bmi
p2 <- observed$sens_glucose

se_log_ratio <- sqrt(((1 - p1) / (observed$n_diseased * p1)) + ((1 - p2) / (observed$n_diseased * p2)))
log_ratio <- log(observed$rTPF)

ci_rTPF <- c(
  lower = exp(log_ratio - z_crit * se_log_ratio),
  upper = exp(log_ratio + z_crit * se_log_ratio)
)

non_inferiority_margin <- 0.8

if (observed$sens_glucose > observed$sens_bmi) {
  better_method <- "Glucose method has higher observed sensitivity than BMI method."
} else if (observed$sens_glucose < observed$sens_bmi) {
  better_method <- "BMI method has higher observed sensitivity than glucose method."
} else {
  better_method <- "Glucose and BMI methods have equal observed sensitivity."
}


## ===============================
## 5) Combined forest-style plot data
## ===============================
plot_df <- data.frame(
  Metric = factor(
    c("Relative true positive fraction (BMI / glucose)", "Sensitivity (glucose method)", "Sensitivity (BMI method)"),
    levels = c("Relative true positive fraction (BMI / glucose)", "Sensitivity (glucose method)", "Sensitivity (BMI method)")
  ),
  Estimate = c(observed$rTPF, observed$sens_glucose, observed$sens_bmi),
  CI_Lower = c(ci_rTPF["lower"], ci_sens_glucose["lower"], ci_sens_bmi["lower"]),
  CI_Upper = c(ci_rTPF["upper"], ci_sens_glucose["upper"], ci_sens_bmi["upper"])
)


## ===============================
## 6) Combined plot
## ===============================
combined_plot <- ggplot(plot_df, aes(x = Estimate, y = Metric, xmin = CI_Lower, xmax = CI_Upper)) +
  geom_point(size = 4, color = "#2c7fb8") +
  geom_errorbarh(height = 0.25, linewidth = 0.9, color = "#2c7fb8") +
  geom_vline(xintercept = non_inferiority_margin, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_vline(xintercept = 1.0, linetype = "dotted", color = "gray40", linewidth = 0.8) +
  labs(
    x = "Estimate with 95% CI",
    y = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    axis.title.x = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey80", linetype = "dotted")
  ) +
  coord_cartesian(xlim = c(0.5, 2))

combined_plot


## ===============================
## 7) Export outputs
## ===============================
ggsave(
  filename = "relative_tpf_combined_plot.png",
  plot = combined_plot,
  width = 10,
  height = 5.5,
  dpi = 400
)

summary_metrics <- data.frame(
  metric = c("rTPF_BMI_over_Glucose", "Sensitivity_glucose", "Sensitivity_BMI"),
  estimate = c(observed$rTPF, observed$sens_glucose, observed$sens_bmi),
  ci_lower = c(ci_rTPF["lower"], ci_sens_glucose["lower"], ci_sens_bmi["lower"]),
  ci_upper = c(ci_rTPF["upper"], ci_sens_glucose["upper"], ci_sens_bmi["upper"]),
  n_diseased = observed$n_diseased,
  tp_glucose = observed$tp_glucose,
  tp_bmi = observed$tp_bmi
)

utils::write.csv(summary_metrics, "relative_tpf_summary_metrics.csv", row.names = FALSE)
utils::write.csv(analysis_data, "relative_tpf_input_dataset.csv", row.names = FALSE)


## ===============================
## 8) Console summary
## ===============================
cat("Relative TPF workflow complete.\n")
cat("Dataset used: MedDataSets::Pima_tr2_df\n")
cat("Methods compared:\n")
cat("- Glucose method: glu >= 140\n")
cat("- BMI method: bmi >= 30\n")
cat("Definition: rTPF = Sensitivity(BMI) / Sensitivity(Glucose)\n")
cat("Non-inferiority margin for rTPF:", non_inferiority_margin, "\n")
cat("Observed rTPF (BMI / Glucose):", round(observed$rTPF, 3), "\n")
cat("95% CI for rTPF:", paste0("[", round(ci_rTPF["lower"], 3), ", ", round(ci_rTPF["upper"], 3), "]"), "\n")
cat(better_method, "\n")
cat("Outputs:\n")
cat("- relative_tpf_combined_plot.png\n")
cat("- relative_tpf_summary_metrics.csv\n")
cat("- relative_tpf_input_dataset.csv\n")
