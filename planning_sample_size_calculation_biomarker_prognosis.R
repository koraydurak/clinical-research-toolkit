## =============================================================================
## Sample size calculation for planning of a study for a prognostic biomarker
## =============================================================================
## Purpose:
##   Estimate the number of required events for study planning of a prognostic biomarker 
##   
## Reference:
##   Schmoor, C., Sauerbrei, W. and Schumacher, M. (2000),
##   Sample size considerations for the evaluation of prognostic factors in
##   survival analysis. Stat Med, 19: 441-452.
##   https://doi.org/10.1002/(SICI)1097-0258(20000229)19:4<441::AID-SIM349>3.0.CO;2-N
## =============================================================================

## ---- 1) User inputs ---------------------------------------------------------
settings <- list(
  # Set to "meddataset_test" for a public/demo run (no internal data required).
  # Set to "local_study_data" to use your prepared project data.
  data_source = "meddataset_test",

  alpha = 0.05,
  power = 0.90,
  assumed_or = 1.1, # Odds Ratio / Risk ratio for endpoint of interest, e.g. mortality, MACE
                    # Assuming e.g. 10% increase per unit of biomarker
  
  # Core exposure used for event planning calculation (SD taken from this variable).
  exposure_var = "continous_exposure",

  # VIF options.
  use_empirical_vif = TRUE,
  fixed_vif = 1
)

meddataset_settings <- list(
  package = "MedDataSets",
  dataset = "Melanoma_df",
  biomarker_var = "thickness",
  demographic_vars = c("age", "sex")
)


## ---- 2) Small helpers -------------------------------------------------------
install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
}

estimate_exposure_vif <- function(data_input, biomarker_var, demographic_vars) {
  required_vars <- c(biomarker_var, demographic_vars)

  missing_vars <- setdiff(required_vars, names(data_input))
  if (length(missing_vars) > 0) {
    stop(paste("Missing variables for VIF estimation:", paste(missing_vars, collapse = ", ")))
  }

  data_model <- data_input[, required_vars, drop = FALSE]
  data_model <- stats::na.omit(data_model)

  if (nrow(data_model) < 20) {
    stop("Too few complete rows for VIF estimation after NA removal.")
  }

  for (variable in demographic_vars) {
    if (is.character(data_model[[variable]])) {
      data_model[[variable]] <- as.factor(data_model[[variable]])
    }
  }

  biomarker <- data_model[[biomarker_var]]
  rhs <- paste(demographic_vars, collapse = " + ")
  formula_vif <- stats::as.formula(paste(biomarker_var, "~", rhs))

  if (is.factor(biomarker) || is.character(biomarker) || length(unique(biomarker)) == 2) {
    if (is.character(data_model[[biomarker_var]])) {
      data_model[[biomarker_var]] <- as.factor(data_model[[biomarker_var]])
    }
    model_r2 <- stats::glm(formula_vif, family = stats::binomial, data = data_model)
    r_squared <- 1 - (model_r2$deviance / model_r2$null.deviance)
  } else {
    model_r2 <- stats::lm(formula_vif, data = data_model)
    r_squared <- summary(model_r2)$r.squared
  }

  if (is.na(r_squared) || r_squared >= 1) {
    stop("Invalid R-squared for VIF estimation (NA or >= 1).")
  }

  1 / (1 - r_squared)
}


## ---- 3) Load analysis data --------------------------------------------------
if (identical(settings$data_source, "local_study_data")) {
  source("3_Code/1_data_preparation.R")
  analysis_data <- data_prepared
  exposure_var_used <- settings$exposure_var

  vif_data <- data_patients_postopctn
  vif_biomarker_var <- "cardiac_troponin_continous"
  vif_demographic_vars <- c("age", "history_of_mi_bin")
} else if (identical(settings$data_source, "meddataset_test")) {
  install_if_missing(c(meddataset_settings$package))
  suppressPackageStartupMessages(
    library(meddataset_settings$package, character.only = TRUE)
  )

  utils::data(list = meddataset_settings$dataset, package = meddataset_settings$package, envir = environment())
  analysis_data <- get(meddataset_settings$dataset, envir = environment())
  exposure_var_used <- meddataset_settings$biomarker_var

  vif_data <- analysis_data
  vif_biomarker_var <- meddataset_settings$biomarker_var
  vif_demographic_vars <- meddataset_settings$demographic_vars
} else {
  stop("settings$data_source must be either 'local_study_data' or 'meddataset_test'.")
}


## ---- 4) VIF estimation ------------------------------------------------------
exposure_vif <- if (isTRUE(settings$use_empirical_vif)) {
  estimate_exposure_vif(
    data_input = vif_data,
    biomarker_var = vif_biomarker_var,
    demographic_vars = vif_demographic_vars
  )
} else {
  settings$fixed_vif
}


## ---- 5) Core calculation ----------------------------------------------------
if (!exposure_var_used %in% names(analysis_data)) {
  stop(paste("Exposure variable not found:", exposure_var_used))
}

if (settings$assumed_or <= 0) {
  stop("assumed_or must be > 0.")
}

exposure_sd <- stats::sd(analysis_data[[exposure_var_used]], na.rm = TRUE)
if (is.na(exposure_sd) || exposure_sd <= 0) {
  stop("Exposure SD is missing or <= 0. Check exposure variable values.")
}

exposure_logodds_assumed <- log(settings$assumed_or)
z_alpha <- stats::qnorm(1 - settings$alpha / 2)
z_power <- stats::qnorm(settings$power)

event_number_needed <- (z_alpha + z_power)^2 /
  (exposure_sd * exposure_logodds_assumed)^2

event_number_needed_vif_adjusted <- event_number_needed * exposure_vif


## ---- 6) Structured output ---------------------------------------------------
sample_size_result <- list(
  method = "Sample size calculation formula according to Schmoor et al.",
  data_source = settings$data_source,
  vif_model = list(
    biomarker = vif_biomarker_var,
    demographics = vif_demographic_vars
  ),
  assumptions = list(
    alpha = settings$alpha,
    power = settings$power,
    assumed_or = settings$assumed_or,
    exposure_variable = exposure_var_used,
    exposure_sd = exposure_sd,
    exposure_vif = exposure_vif
  ),
  results = list(
    events_required_unadjusted = event_number_needed,
    events_required_vif_adjusted = event_number_needed_vif_adjusted
  )
)

print(sample_size_result)







