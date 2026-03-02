# clinical-research-toolkit
A toolkit of various R scripts that can be useful within clinical research, ranging from data cleaning to statistical analysis and reporting. All scripts were adjusted to free datasets which are available as R package or synthetic datasets. 

## Planning: 
- planning_calendar_randomization.R: Generates a weekday-only randomization schedule (RPA vs SOC) for a full study period using fixed seeds for reproducibility. It prints balance summaries and exports the final schedule to Excel for planning/audit use.

- planning_sample_size_calculation_biomarker_prognosis.R: Estimates required event counts for a prognostic biomarker study using the Schmoor et al. approach with user-defined alpha, power, and assumed effect size. It can run on demo data or local study data, optionally estimates VIF empirically, and returns VIF-adjusted event requirements.

## Reporting: 
- reporting_ESC_algorithm_myocardial_infarction.R: Builds a synthetic 400-patient ESC 0/3h hs-cTnI cohort with rule-out/observe/rule-in groups and controlled NSTEMI rates. It computes sensitivity/NPV/specificity/PPV with CIs, creates a publication-style algorithm diagram, and exports dataset, metrics table, and summary text.

- reporting_patient_flowchart.R: Simulates a screened cohort and applies explicit stepwise exclusion criteria (consent, surgery/stay, screening failure, heart surgery, measurement availability, repeat inclusion). It produces a CONSORT-style flowchart PNG plus a CSV of the synthetic patient-level data.

## Statistics: 
- statistics_area_under_the_curve_timevarying.R: Uses survival::pbc to compare bilirubin vs inverse albumin with both standard 1-year ROC AUC and day-by-day time-dependent AUC up to 365 days. It exports ROC/time-varying plots and CSV files for summary statistics, AUC trajectories, and analysis input data.

- statistics_box_violin_plots.R: Loads Cushings data, keeps the three most frequent diagnosis groups, and visualizes two biomarkers with violin+box+jitter panels. It saves a combined figure and exports grouped distribution summaries and the filtered dataset.

- statistics_cumulative_incidence_plots.R: Creates two 1-year cumulative incidence analyses from Melanoma data: KM-based all-cause mortality and competing-risk CIF for melanoma death vs non-melanoma death. It exports each panel (with risk tables), a combined figure, and CSV summaries at 365 days plus analysis data.

- statistics_descriptive_tables.R: Builds a reproducible Table 1 workflow from MedDataSets (VA_df), including variable relabeling and publication-style formatting. It exports a dataset catalog CSV, treatment-stratified descriptive table CSV, and a formatted DOCX table.

- statistics_relative_true_positive_fraction.R: Performs an rTPF analysis on Pima data comparing BMI-threshold screening against glucose-threshold screening using diabetes status as gold standard. It computes sensitivities and log-Wald CIs, checks against a non-inferiority margin, and exports a combined forest-style plot with metrics/data CSVs.
