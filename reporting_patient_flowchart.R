###############################################################################
# Script: reporting_patient_flowchart.R
# Author: Koray Durak
# Contact: koray.durak@outlook.de
# Last updated: 2026-03-02
#
# Purpose:
#   Build a reproducible patient flow chart from a synthetic test dataset.
#   This script demonstrates transparent exclusion logic that others can test and adapt.
###############################################################################

## ===============================
## 1) Dependencies
## ===============================
required_packages <- c("consort", "dplyr", "glue")

install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
}

install_if_missing(required_packages)
library(consort)
library(dplyr)
library(glue)


## ===============================
## 2) Create synthetic test dataset
## ===============================
set.seed(20260302)
n_screened <- 720

patients <- data.frame(
  patient_id = sprintf("PT_%04d", seq_len(n_screened)),
  age = round(rnorm(n_screened, mean = 67, sd = 11)),
  sex = sample(c("Female", "Male"), n_screened, replace = TRUE),
  consent_available = sample(c(TRUE, FALSE), n_screened, replace = TRUE, prob = c(0.93, 0.07)),
  consent_withdrawn = sample(c(TRUE, FALSE), n_screened, replace = TRUE, prob = c(0.05, 0.95)),
  surgery_performed = sample(c(TRUE, FALSE), n_screened, replace = TRUE, prob = c(0.9, 0.1)),
  hospital_stay_hours = round(rlnorm(n_screened, meanlog = log(70), sdlog = 0.45)),
  screening_failure = sample(c(TRUE, FALSE), n_screened, replace = TRUE, prob = c(0.08, 0.92)),
  heart_surgery = sample(c(TRUE, FALSE), n_screened, replace = TRUE, prob = c(0.12, 0.88)),
  hs_ctnt_gen6_measurements = sample(0:4, n_screened, replace = TRUE, prob = c(0.08, 0.2, 0.34, 0.24, 0.14)),
  prior_inclusion_within_1y = sample(c(TRUE, FALSE), n_screened, replace = TRUE, prob = c(0.09, 0.91))
)


## ===============================
## 3) Define exclusion logic
## ===============================
patients <- patients %>%
  mutate(
    reason_general = case_when(
      !consent_available ~ "No consent available",
      consent_withdrawn ~ "Consent withdrawn",
      !surgery_performed | hospital_stay_hours < 48 ~ "Surgery not performed or hospital stay < 48 h",
      screening_failure ~ "Screening failure",
      heart_surgery ~ "Involvement of heart surgery",
      TRUE ~ "Eligible"
    )
  )

general_eligible <- patients %>% filter(reason_general == "Eligible")

general_eligible <- general_eligible %>%
  mutate(
    reason_final = case_when(
      hs_ctnt_gen6_measurements <= 1 ~ "No or only one hs-cTnT-gen6 measurement",
      prior_inclusion_within_1y ~ "Multiple inclusion within one year",
      TRUE ~ "Included"
    )
  )

dataset_analysis <- general_eligible

## ===============================
## 4) Count each flow-chart step
## ===============================
count_screened <- nrow(patients)

excluded_noconsent <- sum(patients$reason_general == "No consent available")
excluded_consentwithdrawn <- sum(patients$reason_general == "Consent withdrawn")
excluded_nosurgery_or_shortstay <- sum(patients$reason_general == "Surgery not performed or hospital stay < 48 h")
excluded_screening_failure <- sum(patients$reason_general == "Screening failure")
excluded_heartsurgery <- sum(patients$reason_general == "Involvement of heart surgery")
excluded_general_all <- count_screened - nrow(general_eligible)

count_general_eligible <- nrow(general_eligible)

excluded_missing_measurement <- sum(general_eligible$reason_final == "No or only one hs-cTnT-gen6 measurement")
excluded_multiple_inclusion <- sum(general_eligible$reason_final == "Multiple inclusion within one year")
excluded_final_all <- excluded_missing_measurement + excluded_multiple_inclusion

count_final_included <- sum(general_eligible$reason_final == "Included")


## ===============================
## 5) Build text boxes
## ===============================
txt1 <- glue("Patients screened in test cohort (n={count_screened})")

txt1_side <- glue(
  "Excluded (n={excluded_general_all}):\n",
  "• No consent available (n={excluded_noconsent})\n",
  "• Consent withdrawn (n={excluded_consentwithdrawn})\n",
  "• Surgery not performed or hospital stay < 48 h (n={excluded_nosurgery_or_shortstay})\n",
  "• Screening failure (n={excluded_screening_failure})\n",
  "• Involvement of heart surgery (n={excluded_heartsurgery})"
)

txt2 <- glue("Patients meeting general eligibility criteria (n={count_general_eligible})")

txt2_side <- glue(
  "Excluded (n={excluded_final_all}):\n",
  "• No or only one cardiac troponin measurement (n={excluded_missing_measurement})\n",
  "• Multiple inclusion within one year (n={excluded_multiple_inclusion})"
)

txt3 <- glue("Patients included in final analysis (n={count_final_included})")


## ===============================
## 6) Draw and export flow chart
## ===============================
flow_chart <- add_box(txt = txt1) %>%
  add_side_box(txt = txt1_side, side = "right") %>%
  add_box(txt = txt2) %>%
  add_side_box(txt = txt2_side, side = "left") %>%
  add_box(txt = txt3)

plot(flow_chart)

output_figure <- "patient_flowchart_demo.png"
grDevices::png(
  filename = output_figure,
  width = 30,
  height = 25,
  res = 300,
  units = "cm",
  type = "cairo"
)
plot(flow_chart)
dev.off()

output_data <- "patient_flowchart_demo_dataset.csv"
utils::write.csv(patients, output_data, row.names = FALSE)

cat("Flow chart created and saved to:", output_figure, "\n")
cat("Synthetic test dataset saved to:", output_data, "\n")
cat("Final included patients:", count_final_included, "\n")
