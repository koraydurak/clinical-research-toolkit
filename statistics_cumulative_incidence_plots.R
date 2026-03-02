###############################################################################
# Script: statistics_cumulative_incidence_plots.R
# Author: Koray Durak
# Contact: koray.durak@outlook.de
# Last updated: 2026-03-02
#
# Purpose:
#   Reproducible cumulative incidence workflow with two methods:
#   A) All-cause mortality (KM shown as cumulative incidence)
#   B) Melanoma-related mortality with competing risk of non-melanoma death
#
# Data source:
#   MedDataSets::Melanoma_df
###############################################################################

## ===============================
## 1) Dependencies
## ===============================
required_packages <- c(
  "MedDataSets", "dplyr", "ggplot2", "survival", "survminer",
  "tidycmprsk", "ggsurvfit", "scales", "magick"
)

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
library(survival)
library(survminer)
library(tidycmprsk)
library(ggsurvfit)
library(scales)
library(magick)


## ===============================
## 2) Data preparation
## ===============================
data("Melanoma_df", package = "MedDataSets")

analysis_data <- Melanoma_df %>%
  transmute(
    time_days = as.numeric(time),
    status = as.numeric(status),
    thickness = as.numeric(thickness)
  ) %>%
  filter(!is.na(time_days), !is.na(status), !is.na(thickness)) %>%
  mutate(
    biomarker_group = ifelse(thickness >= median(thickness, na.rm = TRUE), "High thickness", "Low thickness"),
    biomarker_group = factor(biomarker_group, levels = c("Low thickness", "High thickness")),

    time_1y = pmin(time_days, 365),

    # All-cause mortality (status 1 or 3) within 1 year
    event_allcause_1y = ifelse(time_days <= 365 & status %in% c(1, 3), 1, 0),

    # Competing-risk endpoint:
    # 0 = censored/alive, 1 = melanoma death (event of interest), 2 = non-melanoma death (competing)
    event_competing_code = case_when(
      time_days <= 365 & status == 1 ~ 1,
      time_days <= 365 & status == 3 ~ 2,
      TRUE ~ 0
    ),
    event_competing = factor(
      event_competing_code,
      levels = c(0, 1, 2),
      labels = c("Censored", "Melanoma death", "Non-melanoma death")
    )
  )

if (nrow(analysis_data) == 0) {
  stop("No observations available after filtering.")
}


## ===============================
## 3) Global plotting settings
## ===============================
group_colors <- c("Low thickness" = "#2ca25f", "High thickness" = "firebrick3")
x_breaks <- c(0, 90, 180, 270, 365)

# Y-axis scaling option:
# - "full_100": always 0% to 100% (default)
# - "auto": data-driven upper limit
y_axis_mode <- c("full_100", "auto")[1]

if (y_axis_mode == "full_100") {
  y_lim <- 1
  y_breaks <- seq(0, 1, by = 0.1)
} else {
  # conservative upper limit for this dataset
  y_lim <- 0.35
  y_breaks <- seq(0, y_lim, by = 0.05)
}


## ===============================
## 4) A - All-cause mortality (KM cumulative incidence)
## ===============================
surv_obj <- with(analysis_data, survival::Surv(time_1y, event_allcause_1y))
fit_allcause <- survival::survfit(surv_obj ~ biomarker_group, data = analysis_data)

plot_allcause <- survminer::ggsurvplot(
  fit_allcause,
  data = analysis_data,
  fun = "event",
  axes.offset = TRUE,
  ylim = c(0, y_lim),
  ylab = "Cumulative incidence",
  break.y.by = ifelse(y_axis_mode == "full_100", 0.1, 0.05),
  xlim = c(0, 365),
  xlab = "Follow up (days)",
  break.time.by = 90,
  censor = FALSE,
  risk.table = TRUE,
  risk.table.title = "At risk",
  tables.height = 0.14,
  tables.theme = theme_cleantable(),
  tables.y.text = FALSE,
  fontsize = 4,
  size = 0.75,
  linetype = 1,
  palette = c("#2ca25f", "firebrick3"),
  legend = "top",
  legend.title = "",
  legend.labs = c("Low thickness", "High thickness")
)

plot_allcause$plot <- plot_allcause$plot +
  labs(title = "A - All-cause mortality") +
  scale_y_continuous(
    limits = c(0, y_lim),
    breaks = y_breaks,
    labels = percent_format(accuracy = 1)
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 13, face = "bold", margin = margin(r = 10)),
    axis.text.x = element_text(size = 14, face = "bold", color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.6),
    plot.title = element_text(size = 14, face = "bold", hjust = 0)
  )

## ===============================
## 5) B - Melanoma-related mortality with competing risk
## ===============================
cif_fit <- tidycmprsk::cuminc(
  survival::Surv(time_1y, event_competing) ~ biomarker_group,
  data = analysis_data
)

plot_competing <- cif_fit %>%
  ggsurvfit::ggcuminc(outcome = "Melanoma death", linewidth = 0.75) +
  labs(
    title = "B - Melanoma-related mortality",
    x = "Follow up (days)",
    y = "Cumulative incidence"
  ) +
  scale_x_continuous(
    limits = c(0, 365),
    breaks = x_breaks
  ) +
  scale_y_continuous(
    limits = c(0, y_lim),
    breaks = y_breaks,
    labels = percent_format(accuracy = 1)
  ) +
  scale_color_manual(
    values = c("Low thickness" = "#2ca25f", "High thickness" = "firebrick3"),
    labels = c("Low thickness", "High thickness")
  ) +
  ggsurvfit::add_risktable(
    risktable_stats = "n.risk",
    risktable_title = "At risk",
    size = 4
  ) +
  ggsurvfit::add_risktable_strata_symbol() +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 13, face = "bold", margin = margin(r = 10)),
    axis.text.x = element_text(size = 14, face = "bold", color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.6),
    plot.title = element_text(size = 14, face = "bold", hjust = 0),
    ggsurvfit.risktable.title = element_text(size = 11, face = "bold", color = "black"),
    ggsurvfit.risktable.text = element_text(size = 11, color = "black"),
    ggsurvfit.risktable.y.text = element_blank()
  )


## ===============================
## 6) Export single plots (with risk tables)
## ===============================
out_allcause <- "CIF_allcause_1year_km.png"
out_competing <- "CIF_competingrisk_1year_melanoma_vs_otherdeath.png"
out_combined <- "Figure_combined_CIF_1year.png"

# Save ggsurvplot (includes risk table)
png(filename = out_allcause, width = 8, height = 5.8, units = "in", res = 300, type = "cairo")
print(plot_allcause)
dev.off()

# Save competing-risk plot (includes risk table)
ggsave(
  filename = out_competing,
  plot = plot_competing,
  width = 8,
  height = 5.8,
  dpi = 300
)


## ===============================
## 7) Combine A and B like original workflow
## ===============================
img_A <- magick::image_read(out_allcause)
img_B <- magick::image_read(out_competing)
combined_img <- magick::image_append(c(img_A, img_B), stack = TRUE)
magick::image_write(combined_img, path = out_combined, format = "png")


## ===============================
## 8) Optional summary export
## ===============================
km_summary <- summary(fit_allcause, times = 365)

summary_allcause <- data.frame(
  analysis = "All-cause mortality (KM cumulative incidence)",
  group = sub("biomarker_group=", "", km_summary$strata),
  time_days = km_summary$time,
  cumulative_incidence = 1 - km_summary$surv,
  stringsAsFactors = FALSE
)

# extract CIF estimates at nearest time <=365
extract_cif_365 <- function(group_name, fit_obj) {
  key <- paste0(group_name, " 1")
  curve <- fit_obj[[key]]
  if (is.null(curve)) return(NA_real_)

  valid_idx <- which(curve$time <= 365)
  if (length(valid_idx) == 0) return(0)
  curve$est[max(valid_idx)]
}

summary_competing <- data.frame(
  analysis = "Melanoma-related mortality (competing risk CIF)",
  group = levels(analysis_data$biomarker_group),
  time_days = 365,
  cumulative_incidence = c(
    extract_cif_365("Low thickness", cif_fit),
    extract_cif_365("High thickness", cif_fit)
  ),
  stringsAsFactors = FALSE
)

summary_out <- bind_rows(summary_allcause, summary_competing)
utils::write.csv(summary_out, "CIF_summary_365days.csv", row.names = FALSE)
utils::write.csv(analysis_data, "CIF_input_dataset_melanoma.csv", row.names = FALSE)


## ===============================
## 9) Console summary
## ===============================
cat("Cumulative incidence workflow complete.\n")
cat("Dataset used: MedDataSets::Melanoma_df\n")
cat("Analyses:\n")
cat("- A: All-cause mortality (KM, cumulative incidence)\n")
cat("- B: Melanoma-related mortality (competing-risk CIF)\n")
cat("Outputs:\n")
cat("-", out_allcause, "\n")
cat("-", out_competing, "\n")
cat("-", out_combined, "\n")
cat("- CIF_summary_365days.csv\n")
cat("- CIF_input_dataset_melanoma.csv\n")
