###############################################################################
# Script: statistics_area_under_the_curve_timevarying.R
# Author: Koray Durak
# Contact: koray.durak@outlook.de
# Last updated: 2026-03-02
#
# Purpose:
#   Reproducible AUC workflow with exactly two analyses:
#   1) Standard ROC AUC at 1 year
#   2) Time-varying AUC over 1 year
#
# Dataset:
#   survival::pbc (public clinical dataset with biomarkers and time-to-event)
#
# Biomarkers compared in BOTH analyses:
#   - bilirubin (bili)
#   - albumin (albumin; direction-aligned as negative albumin)
###############################################################################

## ===============================
## 1) Dependencies
## ===============================
required_packages <- c("survival", "dplyr", "ggplot2", "pROC", "timeROC")

install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
}

install_if_missing(required_packages)

library(survival)
library(dplyr)
library(ggplot2)
library(pROC)
library(timeROC)


## ===============================
## 2) Load and prepare medical data
## ===============================
data("pbc", package = "survival")

analysis_data <- pbc %>%
  select(time, status, bili, albumin) %>%
  filter(
    !is.na(time),
    !is.na(status),
    !is.na(bili),
    !is.na(albumin)
  ) %>%
  mutate(
    # Event of interest: death (status == 2)
    event_death = ifelse(status == 2, 1, 0),

    # One-year endpoint for standard ROC
    event_death_1y = ifelse(time <= 365 & event_death == 1, 1, 0),

    # Censor follow-up at 1 year for time-varying AUC
    time_1y = pmin(time, 365),
    event_death_1y_time = ifelse(time <= 365 & event_death == 1, 1, 0),

    # Direction alignment: higher value = higher risk for both biomarkers
    neg_albumin = -albumin
  )

if (sum(analysis_data$event_death_1y) == 0) {
  stop("No one-year death events found after filtering.")
}


## ===============================
## 3) Standard ROC AUC (1-year death)
## ===============================
roc_bili <- pROC::roc(
  response = analysis_data$event_death_1y,
  predictor = analysis_data$bili,
  ci = TRUE,
  quiet = TRUE,
  direction = "<"
)

roc_albumin <- pROC::roc(
  response = analysis_data$event_death_1y,
  predictor = analysis_data$neg_albumin,
  ci = TRUE,
  quiet = TRUE,
  direction = "<"
)

roc_comparison <- pROC::roc.test(roc_bili, roc_albumin)

auc_bili <- as.numeric(pROC::auc(roc_bili))
auc_albumin <- as.numeric(pROC::auc(roc_albumin))
ci_bili <- roc_bili$ci
ci_albumin <- roc_albumin$ci

label_bili <- sprintf(
  "Bilirubin\nAUC %.3f (95%% CI %.3f-%.3f)",
  auc_bili, ci_bili[1], ci_bili[3]
)

label_albumin <- sprintf(
  "Albumin (inverse)\nAUC %.3f (95%% CI %.3f-%.3f)",
  auc_albumin, ci_albumin[1], ci_albumin[3]
)

roc_plot <- pROC::ggroc(
  list(label_bili = roc_bili, label_albumin = roc_albumin),
  size = 1.0
) +
  geom_abline(intercept = 1, slope = 1, linetype = "dashed", color = "gray") +
  scale_color_manual(
    values = c("label_bili" = "firebrick3", "label_albumin" = "#0072B2"),
    labels = c(label_bili, label_albumin),
    name = "Biomarker"
  ) +
  labs(
    x = "Specificity",
    y = "Sensitivity"
  ) +
  theme_bw() +
  theme(
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 11)
  )

roc_plot

## ===============================
## 4) Time-varying AUC (up to 1 year)
## ===============================
times_1y <- seq(1, 365, by = 1)

td_auc_bili <- timeROC(
  T = analysis_data$time_1y,
  delta = analysis_data$event_death_1y_time,
  marker = analysis_data$bili,
  cause = 1,
  times = times_1y,
  iid = TRUE
)

td_auc_albumin <- timeROC(
  T = analysis_data$time_1y,
  delta = analysis_data$event_death_1y_time,
  marker = analysis_data$neg_albumin,
  cause = 1,
  times = times_1y,
  iid = TRUE
)

td_auc_df <- data.frame(
  Time = rep(times_1y, 2),
  AUC = c(td_auc_bili$AUC, td_auc_albumin$AUC),
  SD = c(td_auc_bili$inference$vect_sd_1, td_auc_albumin$inference$vect_sd_1),
  Biomarker = factor(
    rep(c("Bilirubin", "Albumin (inverse)"), each = length(times_1y)),
    levels = c("Bilirubin", "Albumin (inverse)")
  )
)

td_auc_df <- td_auc_df %>%
  mutate(
    Lower = AUC - 1.96 * SD,
    Upper = AUC + 1.96 * SD
  )

# Keep raw AUC values for export, but for plotting fill only leading missing
# values with the first estimable AUC so curves visibly start at day 1.
fill_leading_auc <- function(df) {
  idx <- which(!is.na(df$AUC))[1]
  if (!is.na(idx) && idx > 1) {
    df$AUC[1:(idx - 1)] <- df$AUC[idx]
  }
  df
}

td_auc_df_plot <- td_auc_df %>%
  group_by(Biomarker) %>%
  group_modify(~ fill_leading_auc(.x)) %>%
  ungroup()

bili_365 <- td_auc_df %>% filter(Biomarker == "Bilirubin", Time == 365)
alb_365 <- td_auc_df %>% filter(Biomarker == "Albumin (inverse)", Time == 365)

label_td_bili <- sprintf(
  "Bilirubin",
  bili_365$AUC[1], bili_365$Lower[1], bili_365$Upper[1]
)

label_td_albumin <- sprintf(
  "Albumin (inverse)",
  alb_365$AUC[1], alb_365$Lower[1], alb_365$Upper[1]
)

td_auc_plot <- ggplot(td_auc_df_plot, aes(x = Time, y = AUC, color = Biomarker)) +
  geom_line(size = 1.0) +
  scale_color_manual(
    values = c("Bilirubin" = "firebrick3", "Albumin (inverse)" = "#0072B2"),
    labels = c("Bilirubin" = label_td_bili, "Albumin (inverse)" = label_td_albumin)
  ) +
  labs(
    x = "Days from baseline",
    y = "Time-dependent AUC"
  ) +
  coord_cartesian(ylim = c(0.5, 1.0)) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 11)
  )

td_auc_plot

## ===============================
## 5) Export outputs
## ===============================
ggsave(
  filename = "ROC_AUC_1year_death_two_biomarkers.png",
  plot = roc_plot,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  filename = "TimeVarying_AUC_1year_death_two_biomarkers.png",
  plot = td_auc_plot,
  width = 8,
  height = 6,
  dpi = 300
)

summary_auc <- data.frame(
  analysis = c("Standard ROC (1y death)", "Standard ROC (1y death)", "ROC comparison"),
  biomarker = c("Bilirubin", "Albumin (inverse)", "Bilirubin vs Albumin (inverse)"),
  auc_or_stat = c(
    round(auc_bili, 4),
    round(auc_albumin, 4),
    round(as.numeric(roc_comparison$statistic), 4)
  ),
  ci_lower = c(round(ci_bili[1], 4), round(ci_albumin[1], 4), NA),
  ci_upper = c(round(ci_bili[3], 4), round(ci_albumin[3], 4), NA),
  p_value = c(NA, NA, round(roc_comparison$p.value, 6))
)

utils::write.csv(summary_auc, "AUC_summary_1year_two_biomarkers.csv", row.names = FALSE)
utils::write.csv(td_auc_df, "TimeVarying_AUC_values_1year_two_biomarkers.csv", row.names = FALSE)
utils::write.csv(analysis_data, "AUC_input_dataset_pbc_1year.csv", row.names = FALSE)


## ===============================
## 6) Console summary
## ===============================
cat("AUC workflow complete.\n")
cat("Dataset used: survival::pbc\n")
cat("Biomarkers compared: bilirubin vs inverse albumin\n")
cat("Outputs:\n")
cat("- ROC_AUC_1year_death_two_biomarkers.png\n")
cat("- TimeVarying_AUC_1year_death_two_biomarkers.png\n")
cat("- AUC_summary_1year_two_biomarkers.csv\n")
cat("- TimeVarying_AUC_values_1year_two_biomarkers.csv\n")
cat("- AUC_input_dataset_pbc_1year.csv\n")
