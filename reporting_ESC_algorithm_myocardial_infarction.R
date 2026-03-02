###############################################################################
# Script: reporting_ESC_algorithm_myocardial_infarction.R
# Author: Koray Durak
# Contact: koray.durak@outlook.de
# Last updated: 2026-03-02
#
# Purpose:
#   Reproducible ESC 0/3h algorithm reporting workflow using a synthetic cohort.
#
# Design targets requested:
#   - n = 400 patients
#   - Rule-out group ~40% with only 1-2 NSTEMI
#   - Observe zone with intermediate NSTEMI proportion
#   - Rule-in group with highest NSTEMI proportion
###############################################################################

## ===============================
## 1) Dependencies
## ===============================
required_packages <- c("dplyr", "glue", "DiagrammeR")

install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
}

install_if_missing(required_packages)

library(dplyr)
library(glue)
library(DiagrammeR)


## ===============================
## 2) Synthetic ESC 0/3h algorithm cohort (n = 400)
## ===============================
set.seed(20260302)

n_total <- 400
n_rule_out <- 160   # 40.0%
n_observe <- 160    # 40.0%
n_rule_in <- 80     # 20.0%

# Controlled NSTEMI counts by triage group
# - Rule-out: very low (2)
# - Observe: intermediate (20)
# - Rule-in: highest (48)
n_nstemi_rule_out <- 2
n_nstemi_observe <- 20
n_nstemi_rule_in <- 48

n_nstemi_total <- n_nstemi_rule_out + n_nstemi_observe + n_nstemi_rule_in

triage <- c(
  rep("Rule-out", n_rule_out),
  rep("Observe", n_observe),
  rep("Rule-in", n_rule_in)
)

nstemi <- c(
  c(rep(1, n_nstemi_rule_out), rep(0, n_rule_out - n_nstemi_rule_out)),
  c(rep(1, n_nstemi_observe), rep(0, n_observe - n_nstemi_observe)),
  c(rep(1, n_nstemi_rule_in), rep(0, n_rule_in - n_nstemi_rule_in))
)

synthetic_data <- data.frame(
  patient_id = sprintf("ESC_%03d", seq_len(n_total)),
  triage_group = triage,
  nstemi = nstemi,
  stringsAsFactors = FALSE
)

# Shuffle rows for realistic patient ordering
synthetic_data <- synthetic_data[sample.int(n_total), ]

# Add synthetic clinical/biomarker fields consistent with triage logic
synthetic_data <- synthetic_data %>%
  mutate(
    hs_ctni_0h = case_when(
      triage_group == "Rule-out" ~ runif(n(), 1, 10.5),
      triage_group == "Observe" ~ runif(n(), 8, 45),
      triage_group == "Rule-in" ~ runif(n(), 45, 220),
      TRUE ~ runif(n(), 1, 220)
    ),
    hs_ctni_3h = case_when(
      triage_group == "Rule-out" ~ hs_ctni_0h + runif(n(), -2, 2),
      triage_group == "Observe" ~ hs_ctni_0h + runif(n(), -6, 14),
      triage_group == "Rule-in" ~ hs_ctni_0h + runif(n(), 5, 35),
      TRUE ~ hs_ctni_0h
    ),
    delta_0_3h = abs(hs_ctni_3h - hs_ctni_0h),
    chest_pain_over_6h = ifelse(triage_group == "Rule-out", rbinom(n(), 1, 0.45), rbinom(n(), 1, 0.2)),
    grace_lt_140 = ifelse(triage_group == "Rule-out", rbinom(n(), 1, 0.92), ifelse(triage_group == "Observe", rbinom(n(), 1, 0.7), rbinom(n(), 1, 0.5))),
    symptoms_resolved = ifelse(triage_group == "Rule-out", rbinom(n(), 1, 0.9), ifelse(triage_group == "Observe", rbinom(n(), 1, 0.55), rbinom(n(), 1, 0.25)))
  ) %>%
  mutate(
    triage_binary_rule_out = ifelse(triage_group == "Rule-out", 1, 0),
    triage_binary_rule_in = ifelse(triage_group == "Rule-in", 1, 0)
  )


## ===============================
## 3) Performance metrics + CIs
## ===============================
format_pct <- function(x) sprintf("%.1f%%", 100 * x)
format_ci <- function(x) sprintf("%.1f-%.1f", 100 * x[1], 100 * x[2])

# Rule-out safety: treat "not rule-out" as test positive for NSTEMI detection
fn_rule_out <- sum(synthetic_data$triage_group == "Rule-out" & synthetic_data$nstemi == 1)
tp_not_rule_out <- sum(synthetic_data$triage_group != "Rule-out" & synthetic_data$nstemi == 1)

tn_rule_out <- sum(synthetic_data$triage_group == "Rule-out" & synthetic_data$nstemi == 0)
all_rule_out <- sum(synthetic_data$triage_group == "Rule-out")

sens_rule_out <- tp_not_rule_out / (tp_not_rule_out + fn_rule_out)
npv_rule_out <- tn_rule_out / all_rule_out

ci_sens_rule_out <- stats::binom.test(tp_not_rule_out, tp_not_rule_out + fn_rule_out)$conf.int
ci_npv_rule_out <- stats::binom.test(tn_rule_out, all_rule_out)$conf.int

# Rule-in efficacy: test positive = rule-in
non_nstemi_total <- sum(synthetic_data$nstemi == 0)
tn_not_rule_in <- sum(synthetic_data$triage_group != "Rule-in" & synthetic_data$nstemi == 0)

spec_rule_in <- tn_not_rule_in / non_nstemi_total
ppv_rule_in <- sum(synthetic_data$triage_group == "Rule-in" & synthetic_data$nstemi == 1) / sum(synthetic_data$triage_group == "Rule-in")

ci_spec_rule_in <- stats::binom.test(tn_not_rule_in, non_nstemi_total)$conf.int
ci_ppv_rule_in <- stats::binom.test(
  sum(synthetic_data$triage_group == "Rule-in" & synthetic_data$nstemi == 1),
  sum(synthetic_data$triage_group == "Rule-in")
)$conf.int

# Group summaries
summary_by_group <- synthetic_data %>%
  group_by(triage_group) %>%
  summarise(
    n = n(),
    nstemi_n = sum(nstemi == 1),
    nstemi_pct = 100 * mean(nstemi == 1),
    .groups = "drop"
  )

n_rule_out_actual <- summary_by_group$n[summary_by_group$triage_group == "Rule-out"]
n_observe_actual <- summary_by_group$n[summary_by_group$triage_group == "Observe"]
n_rule_in_actual <- summary_by_group$n[summary_by_group$triage_group == "Rule-in"]

nstemi_rule_out_actual <- summary_by_group$nstemi_n[summary_by_group$triage_group == "Rule-out"]
nstemi_observe_actual <- summary_by_group$nstemi_n[summary_by_group$triage_group == "Observe"]
nstemi_rule_in_actual <- summary_by_group$nstemi_n[summary_by_group$triage_group == "Rule-in"]


## ===============================
## 4) Report-style summary sentence
## ===============================
summary_sentence <- glue(
  "In the synthetic main cohort (n={n_total}; {format_pct(n_nstemi_total / n_total)} with NSTEMI), ",
  "the ESC 0/3h algorithm rule-out criteria classified {n_rule_out_actual} patients ({format_pct(n_rule_out_actual / n_total)}) as rule-out, ",
  "with sensitivity {format_pct(sens_rule_out)} (95% CI {format_ci(ci_sens_rule_out)}) and NPV {format_pct(npv_rule_out)} ",
  "(95% CI {format_ci(ci_npv_rule_out)}). ",
  "Rule-in criteria classified {n_rule_in_actual} patients ({format_pct(n_rule_in_actual / n_total)}) as rule-in, ",
  "with specificity {format_pct(spec_rule_in)} (95% CI {format_ci(ci_spec_rule_in)}) and PPV {format_pct(ppv_rule_in)} ",
  "(95% CI {format_ci(ci_ppv_rule_in)}). ",
  "A total of {n_observe_actual} patients ({format_pct(n_observe_actual / n_total)}) remained in the observe zone, ",
  "of whom {nstemi_observe_actual} ({format_pct(nstemi_observe_actual / n_observe_actual)}) had adjudicated NSTEMI."
)

cat(summary_sentence, "\n\n")


## ===============================
## 5) ESC 0/3h algorithm diagram
## ===============================
uln <- 11
rule_in_abs_cutoff <- 55
rule_in_delta_cutoff <- 11

n_lowuln <- sum(synthetic_data$hs_ctni_0h < uln, na.rm = TRUE)
n_highuln <- sum(synthetic_data$hs_ctni_0h >= uln, na.rm = TRUE)
lowuln_pp <- sprintf("Patients: %s (%s)", n_lowuln, format_pct(n_lowuln / n_total))
highuln_pp <- sprintf("Patients: %s (%s)", n_highuln, format_pct(n_highuln / n_total))

cohort_label <- sprintf("Main cohort\\n(n = %d)", n_total)

rule_out_label <- glue(
  "Rule-out\\n",
  "Patients: {n_rule_out_actual} ({format_pct(n_rule_out_actual / n_total)})\\n",
  "NSTEMI: {nstemi_rule_out_actual} ({format_pct(nstemi_rule_out_actual / n_rule_out_actual)})\\n",
  "Sensitivity: {format_pct(sens_rule_out)}\\n",
  "NPV: {format_pct(npv_rule_out)}"
)

observe_label <- glue(
  "Observe zone\\n",
  "Patients: {n_observe_actual} ({format_pct(n_observe_actual / n_total)})\\n",
  "NSTEMI: {nstemi_observe_actual} ({format_pct(nstemi_observe_actual / n_observe_actual)})"
)

rule_in_label <- glue(
  "Rule-in\\n",
  "Patients: {n_rule_in_actual} ({format_pct(n_rule_in_actual / n_total)})\\n",
  "NSTEMI: {nstemi_rule_in_actual} ({format_pct(nstemi_rule_in_actual / n_rule_in_actual)})\\n",
  "Specificity: {format_pct(spec_rule_in)}\\n",
  "PPV: {format_pct(ppv_rule_in)}"
)

esc_diagram <- DiagrammeR::grViz(sprintf(
  "digraph esc_03h_algorithm {
     graph [layout = dot, rankdir = TB, splines = ortho]

     node [shape = box, style = filled, fontname = Helvetica, fontsize = 20]

     cohort [
       label = '%s',
       fillcolor = '#004C97',
       fontcolor = 'white',
       width = 3
     ]

     re3h_low [
       label = 'hs-cTnI < %s ng/L at 0h\\n%s',
       fillcolor = '#AECDE8',
       width = 3
     ]

     re3h_high [
       label = 'hs-cTnI >= %s ng/L at 0h\\n%s',
       fillcolor = '#AECDE8',
       width = 3
     ]

     { rank = same; re3h_low; re3h_high }

     pre_ruleout [
       label = 'Rule-out if\\n(CPO > 6h or\\nhs-cTn at 3h < %s ng/L)\\nAND\\nGRACE score <140\\nAND\\npainfree',
       fillcolor = '#A5D6A7',
       width = 3
     ]

     pre_observe [
       label = '',
       shape = box,
       style = invis,
       width = 3
     ]

     pre_rulein [
       label = 'Rule-in if\\nhs-cTnI >= %s ng/L or\\nDelta0/3h >= %s ng/L',
       fillcolor = '#EF9A9A',
       width = 3
     ]

     { rank = same; pre_ruleout; pre_observe; pre_rulein }

     ruleout [
       label = '%s',
       fillcolor = '#4CAF50',
       width = 3
     ]

     observe [
       label = '%s',
       fillcolor = '#F2C94C',
       width = 3
     ]

     rulein [
       label = '%s',
       fillcolor = '#E63946',
       width = 3
     ]

     { rank = same; ruleout; observe; rulein }

     cohort -> re3h_low
     cohort -> re3h_high

     re3h_low  -> pre_ruleout
     re3h_low  -> pre_observe [style = invis]
     re3h_low  -> pre_rulein

     re3h_high -> pre_rulein
     re3h_high -> pre_observe [style = invis]
     re3h_high -> pre_ruleout [style = invis]

     pre_ruleout -> ruleout [style = invis]
     pre_rulein  -> rulein [style = invis]
     pre_observe -> observe [style = invis]
   }",
  cohort_label,
  uln,
  lowuln_pp,
  uln,
  highuln_pp,
  uln,
  rule_in_abs_cutoff,
  rule_in_delta_cutoff,
  rule_out_label,
  observe_label,
  rule_in_label
))

esc_diagram

# Optional export if SVG tooling is available
if (requireNamespace("DiagrammeRsvg", quietly = TRUE) && requireNamespace("rsvg", quietly = TRUE)) {
  svg_txt <- DiagrammeRsvg::export_svg(esc_diagram)
  rsvg::rsvg_png(charToRaw(svg_txt), file = "esc_03h_algorithm_reporting_workflow.png", width = 2600, height = 1800)
}


## ===============================
## 6) Export outputs
## ===============================
metrics_table <- data.frame(
  metric = c(
    "Cohort NSTEMI prevalence",
    "Rule-out proportion",
    "Rule-out sensitivity",
    "Rule-out NPV",
    "Rule-in proportion",
    "Rule-in specificity",
    "Rule-in PPV",
    "Observe proportion",
    "Observe NSTEMI proportion"
  ),
  estimate = c(
    format_pct(n_nstemi_total / n_total),
    format_pct(n_rule_out_actual / n_total),
    format_pct(sens_rule_out),
    format_pct(npv_rule_out),
    format_pct(n_rule_in_actual / n_total),
    format_pct(spec_rule_in),
    format_pct(ppv_rule_in),
    format_pct(n_observe_actual / n_total),
    format_pct(nstemi_observe_actual / n_observe_actual)
  ),
  ci_95 = c(
    "-",
    "-",
    sprintf("%s%%-%s%%", round(100 * ci_sens_rule_out[1], 1), round(100 * ci_sens_rule_out[2], 1)),
    sprintf("%s%%-%s%%", round(100 * ci_npv_rule_out[1], 1), round(100 * ci_npv_rule_out[2], 1)),
    "-",
    sprintf("%s%%-%s%%", round(100 * ci_spec_rule_in[1], 1), round(100 * ci_spec_rule_in[2], 1)),
    sprintf("%s%%-%s%%", round(100 * ci_ppv_rule_in[1], 1), round(100 * ci_ppv_rule_in[2], 1)),
    "-",
    "-"
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(synthetic_data, "esc_03h_algorithm_synthetic_dataset_n400.csv", row.names = FALSE)
utils::write.csv(metrics_table, "esc_03h_algorithm_reporting_metrics.csv", row.names = FALSE)

writeLines(summary_sentence, con = "esc_03h_algorithm_summary_sentence.txt")


## ===============================
## 7) Console summary
## ===============================
cat("ESC 0/3h algorithm reporting workflow complete.\n")
cat("Synthetic cohort size:", n_total, "\n")
cat("Triage distribution:\n")
cat("- Rule-out:", n_rule_out_actual, "patients, NSTEMI:", nstemi_rule_out_actual, "\n")
cat("- Observe:", n_observe_actual, "patients, NSTEMI:", nstemi_observe_actual, "\n")
cat("- Rule-in:", n_rule_in_actual, "patients, NSTEMI:", nstemi_rule_in_actual, "\n")
cat("Outputs:\n")
cat("- esc_03h_algorithm_synthetic_dataset_n400.csv\n")
cat("- esc_03h_algorithm_reporting_metrics.csv\n")
cat("- esc_03h_algorithm_summary_sentence.txt\n")
cat("- esc_03h_algorithm_reporting_workflow.png (if DiagrammeRsvg + rsvg available)\n")
