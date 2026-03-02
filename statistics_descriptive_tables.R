###############################################################################
# Script: statistics_descriptive_tables.R
# Author: Koray Durak
# Contact: koray.durak@outlook.de
# Last updated: 2026-03-02
#
# Purpose:
#   Create reproducible descriptive tables using an online dataset from the
#   MedDataSets package.
#
# What this script does:
#   1) Lists available datasets in MedDataSets (catalog export).
#   2) Loads a stable demo dataset (VA_df).
#   3) Creates one descriptive table (stratified by treatment) with tableone.
#   4) Exports outputs as CSV and DOCX for easy sharing on GitHub.
###############################################################################

## ===============================
## 1) Dependencies
## ===============================
required_packages <- c("MedDataSets", "dplyr", "tableone", "flextable", "officer", "tibble")

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
library(tableone)
library(flextable)
library(officer)
library(tibble)


## ===============================
## Helper functions for publishable tables
## ===============================
variable_labels <- c(
  age = "Age, years",
  Karn = "Karnofsky performance status",
  `diag.time` = "Time since diagnosis, days",
  treat = "Treatment arm",
  cell = "Tumor cell type",
  prior = "Prior therapy",
  status_group = "Vital status"
)

format_tableone_output <- function(table_df, labels_map) {
  df <- table_df

  raw_var <- df$Variable
  raw_clean <- gsub(" \\(%\\)$", "", raw_var)
  raw_clean <- gsub(" \\(median \\[IQR\\]\\)$", "", raw_clean)
  raw_clean <- gsub(" \\(mean \\(SD\\)\\)$", "", raw_clean)

  has_level <- grepl(" = ", raw_clean, fixed = TRUE)
  base_var <- ifelse(has_level, sub(" = .*", "", raw_clean), raw_clean)
  level_label <- ifelse(has_level, sub(".* = ", "", raw_clean), NA_character_)

  pretty_var <- labels_map[base_var]
  pretty_var[is.na(pretty_var)] <- base_var[is.na(pretty_var)]

  vars_with_parent_row <- unique(base_var[!has_level])
  inserted_headers <- character(0)
  out_rows <- vector("list", length = 0)

  for (i in seq_len(nrow(df))) {
    if (has_level[i]) {
      needs_header <- !(base_var[i] %in% vars_with_parent_row) && !(base_var[i] %in% inserted_headers)

      if (needs_header) {
        header_row <- df[i, , drop = FALSE]
        header_row[1, ] <- ""
        header_row$Variable <- pretty_var[i]
        out_rows[[length(out_rows) + 1]] <- header_row
        inserted_headers <- c(inserted_headers, base_var[i])
      }

      level_row <- df[i, , drop = FALSE]
      level_row$Variable <- paste0("  - ", level_label[i])
      out_rows[[length(out_rows) + 1]] <- level_row
    } else {
      parent_row <- df[i, , drop = FALSE]
      parent_row$Variable <- pretty_var[i]
      out_rows[[length(out_rows) + 1]] <- parent_row
    }
  }

  dplyr::bind_rows(out_rows)
}

rename_columns_with_n <- function(df, n_lookup, overall_n = NULL) {
  out <- df
  current_names <- names(out)

  if (!is.null(overall_n) && "Overall" %in% current_names) {
    names(out)[names(out) == "Overall"] <- paste0("Overall (n=", overall_n, ")")
  }

  for (grp in names(n_lookup)) {
    if (grp %in% names(out)) {
      names(out)[names(out) == grp] <- paste0(grp, " (n=", n_lookup[[grp]], ")")
    }
  }

  out
}

make_publishable_flextable <- function(df, caption_text) {
  qflextable(df) %>%
    set_caption(caption = caption_text) %>%
    theme_booktabs() %>%
    fontsize(size = 10, part = "all") %>%
    align(align = "left", j = 1, part = "all") %>%
    align(align = "center", j = 2:ncol(df), part = "all") %>%
    bold(part = "header") %>%
    bold(i = ~ !grepl("^  - ", Variable), j = 1, part = "body") %>%
    bg(bg = "#F0F0F0", part = "header") %>%
    autofit()
}


## ===============================
## 2) List available MedDataSets data
## ===============================
meddatasets_catalog <- as.data.frame(data(package = "MedDataSets")$results)

if (nrow(meddatasets_catalog) == 0) {
  stop("No datasets found in MedDataSets. Please check package installation.")
}

meddatasets_catalog <- meddatasets_catalog %>%
  select(Package, Item, Title)

utils::write.csv(meddatasets_catalog, "meddatasets_catalog.csv", row.names = FALSE)

cat("Number of MedDataSets datasets found:", nrow(meddatasets_catalog), "\n")
cat("Catalog exported to: meddatasets_catalog.csv\n")
cat("First 12 datasets:\n")
print(utils::head(meddatasets_catalog$Item, 12))


## ===============================
## 3) Load online demo dataset
## ===============================
# Chosen dataset:
#   VA_df (Veteran's Administration Lung Cancer Trial)
# Why:
#   - Medical context
#   - Mix of categorical and continuous variables
#   - Good for Table 1 style summaries

data("VA_df", package = "MedDataSets")

va_data <- VA_df %>%
  mutate(
    status_group = ifelse(status == 1, "Died", "Alive"),
    status_group = factor(status_group, levels = c("Alive", "Died")),
    treat = as.factor(treat),
    cell = as.factor(cell),
    prior = ifelse(as.character(prior) %in% c("10", "1", "yes", "Yes", "TRUE"), "Yes", "No"),
    prior = factor(prior, levels = c("No", "Yes"))
  )


## ===============================
## 4) Define variables for summaries
## ===============================
vars_for_tables <- c("age", "Karn", "diag.time", "treat", "cell", "prior")
categorical_vars <- c("treat", "cell", "prior")
continuous_nonnormal <- c("age", "Karn", "diag.time")


## ===============================
## 5) Table: Baseline by treatment arm
## ===============================
table_treatment <- CreateTableOne(
  vars = c("age", "Karn", "diag.time", "status_group", "cell", "prior"),
  strata = "treat",
  factorVars = c("status_group", "cell", "prior"),
  addOverall = TRUE,
  test = FALSE,
  data = va_data
)

table_treatment_df <- print(
  table_treatment,
  noSpaces = TRUE,
  nonnormal = continuous_nonnormal,
  catDigits = 1,
  contDigits = 1,
  quote = FALSE,
  printToggle = FALSE
) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Variable")

treat_counts <- table(va_data$treat)
table_treatment_df <- format_tableone_output(table_treatment_df, variable_labels)
table_treatment_df <- rename_columns_with_n(
  table_treatment_df,
  n_lookup = treat_counts,
  overall_n = nrow(va_data)
)

utils::write.csv(table_treatment_df, "descriptive_table_va_by_treatment.csv", row.names = FALSE)


## ===============================
## 6) Export table to Word
## ===============================
ft_treatment <- make_publishable_flextable(
  table_treatment_df,
  "Table 1. Baseline characteristics by treatment arm (VA_df)"
)

output_doc <- "descriptive_tables_meddatasets_demo.docx"

read_docx() %>%
  body_add_par("Table 1. Baseline characteristics by treatment arm (VA_df)", style = "heading 1") %>%
  body_add_flextable(ft_treatment) %>%
  print(target = output_doc)


## ===============================
## 8) Console summary
## ===============================
cat("\nDescriptive table workflow complete.\n")
cat("Outputs:\n")
cat("- meddatasets_catalog.csv\n")
cat("- descriptive_table_va_by_treatment.csv\n")
cat("- descriptive_tables_meddatasets_demo.docx\n")
