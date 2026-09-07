#!/usr/bin/env Rscript

###############################################################################
# 02_predict_GReX.R
#
# Apply Korean- and GTEx-based whole-blood GReX models and evaluate prediction
# performance.
#
# GTEx model resource
# -------------------
# The GTEx-based GReX model weights were obtained from the GTEx v8 Elastic-Net
# eQTL models distributed through PredictDB:
#
#   https://predictdb.org/post/2021/07/21/gtex-v8-models-on-eqtl-and-sqtl/elastic_net_eqtl.tar
#
# The Whole Blood model extracted from this archive was used for GTEx-based
# GReX prediction.
#
# Evaluation strategy
# -------------------
# 1. Korean GReX model -> Asan+Chosun
#    Prediction performance is taken from the 5-fold outer cross-validation
#    performed during Korean GReX model training (test_R2_avg).
#
# 2. GTEx GReX model -> GTEx
#    Prediction performance is taken from the corresponding cross-validation
#    results generated during GTEx model training.
#
# 3. GTEx GReX model -> Asan+Chosun
# 4. Korean GReX model -> CODA
# 5. GTEx GReX model -> CODA
#    For these external applications, predicted expression is compared with
#    observed expression using gene-wise Pearson correlation, gene-wise R2,
#    and sample-wise Pearson correlation.
#
# 6. Korean and GTEx GReX models -> KoGES+GENIE
#    Predicted expression is generated for downstream TWAS.
#
# PrediXcan implementation
# ------------------------
# Genotype dosage conversion and GReX prediction were performed using scripts
# obtained from the original PrediXcan repository:
#
#   https://github.com/hakyimlab/PrediXcan
#
# Specifically:
#   - convert_plink_to_dosage.py was used to generate PrediXcan-ready dosage
#     files from chromosome-specific PLINK genotype files.
#   - PrediXcan.py was used with the --predict option to generate predicted
#     gene expression from the Korean- and GTEx-based GReX model databases.
#
# These upstream scripts are not redistributed in this repository.
#
# PrediXcan dosage inputs are assumed to have already been generated before
# running this script and to contain chromosome-specific chr*.txt.gz files
# plus samples.txt.
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(DBI)
  library(RSQLite)
})

# =============================================================================
# 1. Configuration
# =============================================================================

PYTHON_BIN <- "python3"

# PrediXcan.py from:
# https://github.com/hakyimlab/PrediXcan
PREDIXCAN_SCRIPT <- "/path/to/PrediXcan/Software/PrediXcan.py"

# Databases used for prediction.
#
# Korean:
#   Korean whole-blood GReX model trained in this study.
#
# GTEx:
#   GTEx v8 Elastic-Net eQTL model downloaded from PredictDB:
#   https://predictdb.org/post/2021/07/21/gtex-v8-models-on-eqtl-and-sqtl/elastic_net_eqtl.tar
#
#   The Whole Blood prediction database extracted from the downloaded archive
#   was used for GTEx-based GReX prediction.
PREDICTION_MODELS <- list(
  Korean = "/path/to/Korean_Whole_Blood_models_filtered_signif.db",
  GTEx = "/path/to/GTEx_v8_ElasticNet_eQTL/en_Whole_Blood.db"
)

# Unfiltered/model-training databases used to extract cross-validation
# performance. These databases are expected to contain a model_summaries table
# with test_R2_avg, test_R2_sd, zscore_pval, and rho_avg.
CV_DATABASES <- list(
  Korean = "/path/to/Korean_Whole_Blood_models.db",
  GTEx = "/path/to/GTEx_Whole_Blood_models.db"
)

# PrediXcan-ready dosage inputs are assumed to be available in each dosage_dir.
# These files were generated beforehand using convert_plink_to_dosage.py from:
# https://github.com/hakyimlab/PrediXcan
#
# Required structure:
#   dosage_dir/
#     chr1.txt.gz ... chr22.txt.gz
#     samples.txt
#
# observed_expression_file is required only for Asan+Chosun and CODA.
DATASETS <- list(
  Asan_Chosun = list(
    dosage_dir = "/path/to/Asan_Chosun/genotype",
    samples_file = "/path/to/Asan_Chosun/genotype/samples.txt",
    observed_expression_file = "/path/to/Asan_Chosun/residuals_expression.txt",
    output_dir = "/path/to/output/Asan_Chosun"
  ),
  CODA = list(
    dosage_dir = "/path/to/CODA/genotype",
    samples_file = "/path/to/CODA/genotype/samples.txt",
    observed_expression_file = "/path/to/CODA/residuals_expression.txt",
    output_dir = "/path/to/output/CODA"
  ),
  KoGES_GENIE = list(
    dosage_dir = "/path/to/KoGES_GENIE/genotype",
    samples_file = "/path/to/KoGES_GENIE/genotype/samples.txt",
    observed_expression_file = NA,
    output_dir = "/path/to/output/KoGES_GENIE"
  )
)

MODEL_FILTER <- "zscore_pval < 0.05 AND rho_avg > 0.1"

# =============================================================================
# 2. Utility functions
# =============================================================================

check_file <- function(path) {
  if (!file.exists(path)) stop("Required file not found: ", path, call. = FALSE)
}

check_dosage_input <- function(dataset) {
  check_file(dataset$samples_file)
  for (chr in 1:22) check_file(file.path(dataset$dosage_dir, paste0("chr", chr, ".txt.gz")))
}

run_command <- function(command, args) {
  status <- system2(command, args = args)
  if (!identical(status, 0L)) stop("Command failed: ", command, call. = FALSE)
}

# =============================================================================
# 3. Cross-validation performance from model training
# =============================================================================

extract_cv_performance <- function(model_name, db_file, output_dir) {
  check_file(db_file)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  conn <- dbConnect(SQLite(), db_file)
  on.exit(if (dbIsValid(conn)) dbDisconnect(conn), add = TRUE)

  if (!"model_summaries" %in% dbListTables(conn)) {
    stop(model_name, " database does not contain a model_summaries table.", call. = FALSE)
  }

  model_summaries <- dbGetQuery(
    conn,
    paste0(
      "SELECT * FROM model_summaries WHERE ",
      MODEL_FILTER
    )
  )

  required_cols <- c("gene", "test_R2_avg")
  if (!all(required_cols %in% colnames(model_summaries))) {
    stop(
      model_name,
      " model_summaries must contain: ",
      paste(required_cols, collapse = ", "),
      call. = FALSE
    )
  }

  keep_cols <- intersect(
    c("gene", "test_R2_avg", "test_R2_sd", "rho_avg", "rho_avg_squared", "zscore_pval"),
    colnames(model_summaries)
  )
  cv_results <- model_summaries[, keep_cols, drop = FALSE]

  cv_summary <- data.frame(
    model = model_name,
    n_models = nrow(cv_results),
    mean_test_R2 = mean(cv_results$test_R2_avg, na.rm = TRUE),
    median_test_R2 = median(cv_results$test_R2_avg, na.rm = TRUE),
    mean_test_R2_sd = if ("test_R2_sd" %in% colnames(cv_results)) mean(cv_results$test_R2_sd, na.rm = TRUE) else NA_real_,
    stringsAsFactors = FALSE
  )

  fwrite(cv_results, file.path(output_dir, paste0(model_name, "_cross_validation_R2.tsv")), sep = "\t")
  fwrite(cv_summary, file.path(output_dir, paste0(model_name, "_cross_validation_summary.tsv")), sep = "\t")

  message(
    "[", model_name, " training CV] ",
    nrow(cv_results), " models; mean test R2 = ",
    signif(cv_summary$mean_test_R2, 4)
  )
  invisible(cv_summary)
}

# =============================================================================
# 4. Apply GReX models with PrediXcan
# =============================================================================

run_predixcan <- function(dataset_name, model_name) {
  dataset <- DATASETS[[dataset_name]]
  model_db <- PREDICTION_MODELS[[model_name]]

  check_file(PREDIXCAN_SCRIPT)
  check_file(model_db)
  check_dosage_input(dataset)

  result_dir <- file.path(dataset$output_dir, "predicted_expression")
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  output_prefix <- file.path(result_dir, paste0(dataset_name, "_", model_name, "_Whole_Blood"))
  predicted_file <- paste0(output_prefix, "_predicted_expression.txt")

  if (!file.exists(predicted_file)) {
    message("[", dataset_name, "] Applying ", model_name, " GReX model...")
    run_command(
      PYTHON_BIN,
      c(
        PREDIXCAN_SCRIPT,
        "--predict",
        "--dosages", dataset$dosage_dir,
        "--dosages_prefix", "chr",
        "--samples", dataset$samples_file,
        "--weights", model_db,
        "--output_prefix", output_prefix
      )
    )
  } else {
    message("[", dataset_name, " / ", model_name, "] Existing prediction found; skipping.")
  }

  check_file(predicted_file)
  predicted_file
}

# =============================================================================
# 5. External prediction performance
# =============================================================================
#
# Based on the comparison strategy used in Comparing.R:
#   - match common samples and genes
#   - gene-wise Pearson correlation
#   - gene-wise R2 from observed expression ~ predicted expression
#   - sample-wise Pearson correlation
# =============================================================================

read_predicted_expression <- function(path) {
  x <- fread(path, data.table = FALSE, check.names = FALSE)
  if (ncol(x) < 3) stop("Unexpected PrediXcan output format: ", path, call. = FALSE)
  sample_ids <- as.character(x[[1]])
  pred <- x[, -(1:2), drop = FALSE]
  rownames(pred) <- sample_ids
  pred
}

read_observed_expression <- function(path) {
  x <- fread(path, data.table = FALSE, check.names = FALSE)
  if (ncol(x) < 2) stop("Unexpected observed-expression format: ", path, call. = FALSE)
  sample_ids <- as.character(x[[1]])
  obs <- x[, -1, drop = FALSE]
  rownames(obs) <- sample_ids
  obs
}

safe_cor <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3 || sd(x[keep]) == 0 || sd(y[keep]) == 0) return(NA_real_)
  cor(x[keep], y[keep], method = "pearson")
}

safe_r2 <- function(observed, predicted) {
  keep <- is.finite(observed) & is.finite(predicted)
  if (sum(keep) < 3 || sd(predicted[keep]) == 0) return(NA_real_)
  summary(lm(observed[keep] ~ predicted[keep]))$r.squared
}

compare_external_prediction <- function(dataset_name, model_name, predicted_file) {
  dataset <- DATASETS[[dataset_name]]
  check_file(dataset$observed_expression_file)

  predicted <- read_predicted_expression(predicted_file)
  observed <- read_observed_expression(dataset$observed_expression_file)

  common_samples <- intersect(rownames(observed), rownames(predicted))
  common_genes <- intersect(colnames(observed), colnames(predicted))
  if (length(common_samples) < 3) stop("Fewer than three overlapping samples for ", dataset_name, ".", call. = FALSE)
  if (length(common_genes) == 0) stop("No overlapping genes for ", dataset_name, ".", call. = FALSE)

  observed <- observed[common_samples, common_genes, drop = FALSE]
  predicted <- predicted[common_samples, common_genes, drop = FALSE]

  gene_cor <- vapply(common_genes, function(g) safe_cor(observed[[g]], predicted[[g]]), numeric(1))
  gene_r2 <- vapply(common_genes, function(g) safe_r2(observed[[g]], predicted[[g]]), numeric(1))
  gene_results <- data.frame(
    gene = common_genes,
    correlation = gene_cor,
    R2 = gene_r2,
    n_samples = length(common_samples),
    stringsAsFactors = FALSE
  )

  sample_cor <- vapply(common_samples, function(id) {
    safe_cor(as.numeric(observed[id, ]), as.numeric(predicted[id, ]))
  }, numeric(1))
  sample_results <- data.frame(
    sample_id = common_samples,
    correlation = sample_cor,
    n_genes = length(common_genes),
    stringsAsFactors = FALSE
  )

  summary_results <- data.frame(
    dataset = dataset_name,
    model = model_name,
    n_samples = length(common_samples),
    n_genes = length(common_genes),
    mean_gene_correlation = mean(gene_cor, na.rm = TRUE),
    median_gene_correlation = median(gene_cor, na.rm = TRUE),
    mean_gene_R2 = mean(gene_r2, na.rm = TRUE),
    median_gene_R2 = median(gene_r2, na.rm = TRUE),
    mean_sample_correlation = mean(sample_cor, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  comparison_dir <- file.path(dataset$output_dir, "prediction_performance", model_name)
  dir.create(comparison_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(gene_results, file.path(comparison_dir, "gene_level_prediction_performance.tsv"), sep = "\t")
  fwrite(sample_results, file.path(comparison_dir, "sample_level_prediction_correlation.tsv"), sep = "\t")
  fwrite(summary_results, file.path(comparison_dir, "prediction_performance_summary.tsv"), sep = "\t")

  message(
    "[", dataset_name, " / ", model_name, "] ",
    length(common_genes), " genes; mean correlation = ",
    signif(summary_results$mean_gene_correlation, 4),
    "; mean R2 = ", signif(summary_results$mean_gene_R2, 4)
  )
  invisible(summary_results)
}

# =============================================================================
# 6. Analysis workflow
# =============================================================================

run_workflow <- function() {
  # A. Internal model-training performance
  extract_cv_performance(
    model_name = "Korean_Asan_Chosun",
    db_file = CV_DATABASES$Korean,
    output_dir = file.path(DATASETS$Asan_Chosun$output_dir, "cross_validation")
  )

  extract_cv_performance(
    model_name = "GTEx_GTEx",
    db_file = CV_DATABASES$GTEx,
    output_dir = file.path(dirname(PREDICTION_MODELS$GTEx), "cross_validation")
  )

  # B. External prediction performance
  pred <- run_predixcan("Asan_Chosun", "GTEx")
  compare_external_prediction("Asan_Chosun", "GTEx", pred)

  pred <- run_predixcan("CODA", "Korean")
  compare_external_prediction("CODA", "Korean", pred)

  pred <- run_predixcan("CODA", "GTEx")
  compare_external_prediction("CODA", "GTEx", pred)

  # C. Predicted expression for downstream TWAS
  run_predixcan("KoGES_GENIE", "Korean")
  run_predixcan("KoGES_GENIE", "GTEx")

  message("All GReX prediction and performance-evaluation steps completed.")
}

# =============================================================================
# 7. Command-line interface
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0 || args[1] == "all") {
  run_workflow()
} else if (args[1] == "cv") {
  extract_cv_performance(
    "Korean_Asan_Chosun",
    CV_DATABASES$Korean,
    file.path(DATASETS$Asan_Chosun$output_dir, "cross_validation")
  )
  extract_cv_performance(
    "GTEx_GTEx",
    CV_DATABASES$GTEx,
    file.path(dirname(PREDICTION_MODELS$GTEx), "cross_validation")
  )
} else if (args[1] == "external") {
  p <- run_predixcan("Asan_Chosun", "GTEx")
  compare_external_prediction("Asan_Chosun", "GTEx", p)
  p <- run_predixcan("CODA", "Korean")
  compare_external_prediction("CODA", "Korean", p)
  p <- run_predixcan("CODA", "GTEx")
  compare_external_prediction("CODA", "GTEx", p)
} else if (args[1] == "twas") {
  run_predixcan("KoGES_GENIE", "Korean")
  run_predixcan("KoGES_GENIE", "GTEx")
} else {
  stop("Usage: Rscript 02_predict_GReX.R [all|cv|external|twas]", call. = FALSE)
}
