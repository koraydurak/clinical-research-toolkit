###############################################################################
# Script: planning_calendar_randomization.R
# Author: Koray Durak
# Contact: koray.durak@outlook.de
# Last updated: 2026-03-02
#
# Purpose:
#   Create reproducible day-level randomization schedules for study planning
#   of a randomized trial where the randomization is performed on working-day level.
#   This script does NOT require sensitive patient-level data.
#
# Reproducibility notes:
#   - Fixed random seeds are used per schedule.
#   - Weekday-only allocation (Mon-Fri) is enforced.
#   - Output is saved to Excel for sharing/audit.
#
# Public data note:
#   - An optional demo section shows how to load an openly available dataset
#     from the MedDataSets package:
#     https://github.com/lightbluetitan/meddatasets
###############################################################################

## ===============================
## 1) Dependencies
## ===============================
required_packages <- c("writexl")

install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
}

install_if_missing(required_packages)
library(writexl)


## ===============================
## 2) Helper functions
## ===============================
create_workday_randomization <- function(start_date,
                                         end_date,
                                         seed,
                                         group_labels = c("RPA", "SOC")) {
  set.seed(seed)

  all_days <- seq(as.Date(start_date), as.Date(end_date), by = "day")
  isoweek <- as.integer(format(all_days, "%u"))
  work_days <- all_days[isoweek %in% 1:5]

  weekday_en <- c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday")[
    as.integer(format(work_days, "%u"))
  ]
  weekday_local <- weekdays(work_days)

  groups <- rep(group_labels, length.out = length(work_days))
  assignment <- sample(groups)

  schedule <- data.frame(
    Date = work_days,
    Weekday_EN = weekday_en,
    Weekday_Local = weekday_local,
    Group = assignment,
    stringsAsFactors = FALSE
  )

  schedule
}

save_schedule <- function(schedule, output_file) {
  write_xlsx(schedule, output_file)
  invisible(output_file)
}

print_schedule_summary <- function(schedule, label = "Schedule") {
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat(label, "\n")
  cat(strrep("-", 70), "\n", sep = "")
  cat("Date range:", as.character(min(schedule$Date)), "to", as.character(max(schedule$Date)), "\n")
  cat("Number of workdays:", nrow(schedule), "\n")
  cat("Group balance:\n")
  print(table(schedule$Group))
  cat("Preview (first 10 rows):\n")
  print(utils::head(schedule, 10))
}


## ===============================
## 3) Randomization plan (single one-year example)
## ===============================
# Edit this block for future studies.
randomization_plan <- list(
  plan_name = "Example plan: Full year 2026",
  start_date = "2026-01-01",
  end_date = "2026-12-31",
  seed = 2026,
  output_file = "RPA_randomization_schedule_2026_full_year.xlsx"
)


## ===============================
## 4) Run plan
## ===============================
schedule <- create_workday_randomization(
  start_date = randomization_plan$start_date,
  end_date = randomization_plan$end_date,
  seed = randomization_plan$seed,
  group_labels = c("RPA", "SOC")
)

print_schedule_summary(schedule, randomization_plan$plan_name)
save_schedule(schedule, randomization_plan$output_file)

cat("\nRandomization planning complete. Output Excel files saved in the project folder.\n")
