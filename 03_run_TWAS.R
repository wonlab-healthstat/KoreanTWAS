#!/usr/bin/env Rscript

###############################################################################
# 03_run_TWAS.R
#
# Individual-level TWAS in KoGES+GENIE using Korean- and GTEx-based GReX.
#
# Phenotype columns are assumed to use the study abbreviations directly
# (e.g. BMI, HDL, BRCA, PROCA) rather than the original "B_" coding names.
#
# Analysis:
#   Continuous traits -> linear regression
#   Binary traits     -> Firth bias-reduced logistic regression
#   SMOKE             -> analyzed separately
#   PROCA, BPH        -> male-only analyses
#   BRCA, UTCA        -> female-only analyses
#
# For both Korean and GTEx models, genes are restricted to:
#   model_summaries$rho_avg_squared > 0.05
#
# The number of GReX models satisfying R2 > 0.05 is used as the
# multiple-testing denominator for each model.
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(DBI)
  library(RSQLite)
  library(logistf)
})

# =============================================================================
# 1. Input files
# =============================================================================

PHENOTYPE_FILE <- "/path/to/KoGES_GENIE_phenotypes.csv"
PHENOTYPE_ID_COL <- "IID"

PREDICTED_EXPRESSION <- list(
  Korean = "/path/to/KoGES_GENIE_Korean_Whole_Blood_predicted_expression.txt",
  GTEx = "/path/to/KoGES_GENIE_GTEx_Whole_Blood_predicted_expression.txt"
)

MODEL_DATABASES <- list(
  Korean = "/path/to/Korean_Whole_Blood_models.db",
  GTEx = "/path/to/GTEx_Whole_Blood_models.db"
)

OUTPUT_DIR <- "/path/to/TWAS_results"
R2_THRESHOLD <- 0.05
MAX_ZERO_RATIO <- 0.90

# =============================================================================
# 2. Trait information
# =============================================================================

CONTINUOUS_TRAITS <- c(
  "BMI", "HEIGHT", "HIP", "WAIST", "DBP", "SBP", "MAP",
  "HCT", "HB", "MCH", "MCHC", "MCV", "PLAT", "RBC", "WBC",
  "CREATININE", "eGFR", "BUN", "ALT", "AST", "R_GTP", "T_BIL",
  "ALBUMIN", "GLU", "HbA1C", "HDL", "LDL", "TCHL", "TG"
)

BINARY_TRAITS <- c(
  "SMOKE", "BRCA", "UTCA", "COLCA", "GCA", "HCCCA", "LCA", "PROCA",
  "THYCA", "ANG", "ARRHY", "CVA", "CAD", "CHF", "HTN", "MI", "PV",
  "STR", "TIA", "VD", "LIV", "GASTRO", "HEPATITISB", "HEPATITISC",
  "FLIV", "GB", "STOMUL", "CIRRHOSIS", "ULCER", "THY", "CATA",
  "GLAU", "KD", "LIP", "DM", "ARTH", "FRAC", "OSTE", "NOI", "PARK",
  "ALLER", "GT", "POL", "PER", "ASTH", "BRON", "COPD", "CLD", "TB",
  "UB", "BPH", "UT"
)

TRAIT_INFO <- data.frame(
  trait = c(CONTINUOUS_TRAITS, BINARY_TRAITS),
  type = c(rep("continuous", length(CONTINUOUS_TRAITS)), rep("binary", length(BINARY_TRAITS))),
  analysis_group = c(
    rep("all", length(CONTINUOUS_TRAITS)),
    ifelse(BINARY_TRAITS == "SMOKE", "smoking",
      ifelse(BINARY_TRAITS %in% c("PROCA", "BPH"), "male",
        ifelse(BINARY_TRAITS %in% c("BRCA", "UTCA"), "female", "all")
      )
    )
  ),
  exclude_bmi_covariate = c(
    CONTINUOUS_TRAITS %in% c("BMI", "HEIGHT", "HIP", "WAIST"),
    rep(FALSE, length(BINARY_TRAITS))
  ),
  stringsAsFactors = FALSE
)

SEX_VAR <- "SEX"
AGE_VAR <- "AGE"
SMOKING_VAR <- "SMOKE"
BMI_VAR <- "BMI"

# =============================================================================
# 3. Utility functions
# =============================================================================

check_file <- function(path) {
  if (!file.exists(path)) stop("Required file not found: ", path, call. = FALSE)
}

check_columns <- function(data, columns, label) {
  missing <- setdiff(columns, colnames(data))
  if (length(missing) > 0) stop(label, " columns not found: ", paste(missing, collapse = ", "), call. = FALSE)
}

read_predicted_expression <- function(path) {
  check_file(path)
  x <- fread(path, data.table = FALSE, check.names = FALSE)
  if (ncol(x) < 3) stop("Unexpected PrediXcan output format: ", path, call. = FALSE)
  sample_ids <- as.character(x[[1]])
  pred <- x[, -(1:2), drop = FALSE]
  rownames(pred) <- sample_ids
  pred
}

# =============================================================================
# 4. R2 filtering of GReX models
# =============================================================================

get_r2_filtered_genes <- function(db_file, r2_threshold = R2_THRESHOLD) {
  check_file(db_file)
  conn <- dbConnect(SQLite(), db_file)
  on.exit(if (dbIsValid(conn)) dbDisconnect(conn), add = TRUE)
  if (!"model_summaries" %in% dbListTables(conn)) stop("model_summaries table not found in: ", db_file, call. = FALSE)

  fields <- dbListFields(conn, "model_summaries")
  if (!all(c("gene", "rho_avg_squared") %in% fields)) {
    stop("model_summaries must contain gene and rho_avg_squared: ", db_file, call. = FALSE)
  }

  x <- dbGetQuery(
    conn,
    paste0(
      "SELECT gene, rho_avg_squared AS R2 FROM model_summaries ",
      "WHERE rho_avg_squared > ", r2_threshold
    )
  )
  x <- x[!is.na(x$R2) & !duplicated(x$gene), , drop = FALSE]
  x
}

select_twas_genes <- function(grex, r2_info) {
  model_genes <- r2_info$gene
  available_genes <- intersect(model_genes, colnames(grex))
  if (length(available_genes) == 0) stop("No R2-filtered genes found in predicted expression.", call. = FALSE)

  zero_ratio <- colSums(grex[, available_genes, drop = FALSE] == 0, na.rm = TRUE) / nrow(grex)
  analysis_genes <- available_genes[which(zero_ratio <= MAX_ZERO_RATIO)]

  list(
    model_genes = model_genes,
    available_genes = available_genes,
    analysis_genes = analysis_genes,
    n_model_genes = length(model_genes),
    n_available_genes = length(available_genes),
    n_analysis_genes = length(analysis_genes)
  )
}

# =============================================================================
# 5. Phenotype preparation
# =============================================================================

load_phenotypes <- function() {
  check_file(PHENOTYPE_FILE)
  pheno <- fread(PHENOTYPE_FILE, data.table = FALSE, check.names = FALSE)
  required <- unique(c(PHENOTYPE_ID_COL, SEX_VAR, AGE_VAR, TRAIT_INFO$trait))
  check_columns(pheno, required, "Phenotype")

  # Original study coding:
  # SEX   : 1 = male, 2 = female -> 0 = male, 1 = female
  # SMOKE : 1 = never, 2/3 = former/current -> 0 = never, 1 = ever
  # Binary disease traits: 1 = no, 2 = yes -> 0 = control, 1 = case
  pheno[[SEX_VAR]] <- ifelse(pheno[[SEX_VAR]] == 1, 0, 1)
  pheno[[SMOKING_VAR]] <- ifelse(pheno[[SMOKING_VAR]] == 1, 0, 1)

  disease_traits <- BINARY_TRAITS[BINARY_TRAITS != SMOKING_VAR]
  for (trait in disease_traits) pheno[[trait]] <- ifelse(pheno[[trait]] == 1, 0, 1)
  pheno
}

merge_grex_with_phenotypes <- function(pheno, grex) {
  common_ids <- intersect(as.character(pheno[[PHENOTYPE_ID_COL]]), rownames(grex))
  if (length(common_ids) == 0) stop("No overlapping samples between phenotype and GReX data.", call. = FALSE)

  pheno <- pheno[match(common_ids, pheno[[PHENOTYPE_ID_COL]]), , drop = FALSE]
  grex <- grex[common_ids, , drop = FALSE]
  cbind(pheno, grex)
}

# =============================================================================
# 6. Continuous TWAS
# =============================================================================

run_linear_trait <- function(data, trait, genes, covariates, n_tests) {
  data <- data[!is.na(data[[trait]]), , drop = FALSE]

  results <- lapply(genes, function(gene) {
    fit <- lm(reformulate(c(gene, covariates), response = trait), data = data)
    coef_table <- summary(fit)$coefficients
    if (!gene %in% rownames(coef_table)) return(NULL)

    data.frame(
      trait = trait,
      gene = gene,
      beta = coef_table[gene, "Estimate"],
      SE = coef_table[gene, "Std. Error"],
      Zscore = coef_table[gene, "Estimate"] / coef_table[gene, "Std. Error"],
      statistic = coef_table[gene, "t value"],
      p = coef_table[gene, "Pr(>|t|)"],
      N = nobs(fit),
      stringsAsFactors = FALSE
    )
  })

  out <- bind_rows(results)
  out$p_bonferroni <- p.adjust(out$p, method = "bonferroni", n = n_tests)
  out$p_fdr <- p.adjust(out$p, method = "fdr", n = n_tests)
  out$bonferroni_threshold <- 0.05 / n_tests
  out$significant_bonferroni <- out$p < out$bonferroni_threshold
  out
}

run_continuous_twas <- function(data, genes, n_tests, model_name) {
  traits <- TRAIT_INFO$trait[which(TRAIT_INFO$type == "continuous")]
  out_dir <- file.path(OUTPUT_DIR, model_name, "continuous")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  for (trait in traits) {
    omit_bmi <- TRAIT_INFO$exclude_bmi_covariate[which(TRAIT_INFO$trait == trait)]
    covariates <- c(SEX_VAR, AGE_VAR, SMOKING_VAR, BMI_VAR)
    if (isTRUE(omit_bmi)) covariates <- c(SEX_VAR, AGE_VAR, SMOKING_VAR)

    message("[", model_name, "] Continuous: ", trait)
    result <- run_linear_trait(data, trait, genes, covariates, n_tests)
    fwrite(result, file.path(out_dir, paste0(trait, ".tsv")), sep = "\t")
  }
}

# =============================================================================
# 7. Binary TWAS using Firth logistic regression
# =============================================================================

run_firth_trait <- function(data, trait, genes, covariates, n_tests) {
  data <- data[!is.na(data[[trait]]), , drop = FALSE]

  results <- lapply(genes, function(gene) {
    fit <- tryCatch(
      logistf::logistf(reformulate(c(gene, covariates), response = trait), data = data),
      error = function(e) NULL
    )
    if (is.null(fit) || !gene %in% names(fit$coefficients)) return(NULL)

    dropped <- tryCatch(drop1(fit), error = function(e) NULL)
    chisq <- NA_real_
    p_value <- NA_real_

    if (!is.null(dropped)) {
      row_id <- if (gene %in% rownames(dropped)) gene else rownames(dropped)[1]
      chisq <- suppressWarnings(as.numeric(dropped[row_id, 1]))
      p_value <- suppressWarnings(as.numeric(dropped[row_id, 3]))
    }

    beta_value <- fit$coefficients[gene]
    se_value <- sqrt(diag(vcov(fit)))[gene]

    data.frame(
      trait = trait,
      gene = gene,
      beta = beta_value,
      SE = se_value,
      Zscore = beta_value / se_value,
      Chisq = chisq,
      p = p_value,
      N = nrow(data),
      N_case = sum(data[[trait]] == 1, na.rm = TRUE),
      N_control = sum(data[[trait]] == 0, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })

  out <- bind_rows(results)
  out$p_bonferroni <- p.adjust(out$p, method = "bonferroni", n = n_tests)
  out$p_fdr <- p.adjust(out$p, method = "fdr", n = n_tests)
  out$bonferroni_threshold <- 0.05 / n_tests
  out$significant_bonferroni <- out$p < out$bonferroni_threshold
  out
}

run_binary_twas <- function(data, genes, n_tests, model_name) {
  smoking_traits <- TRAIT_INFO$trait[which(TRAIT_INFO$type == "binary" & TRAIT_INFO$analysis_group == "smoking")]
  all_traits <- TRAIT_INFO$trait[which(TRAIT_INFO$type == "binary" & TRAIT_INFO$analysis_group == "all")]
  male_traits <- TRAIT_INFO$trait[which(TRAIT_INFO$type == "binary" & TRAIT_INFO$analysis_group == "male")]
  female_traits <- TRAIT_INFO$trait[which(TRAIT_INFO$type == "binary" & TRAIT_INFO$analysis_group == "female")]

  for (trait in smoking_traits) {
    out_dir <- file.path(OUTPUT_DIR, model_name, "binary", "smoking")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    message("[", model_name, "] Binary (smoking): ", trait)
    result <- run_firth_trait(data, trait, genes, c(SEX_VAR, AGE_VAR), n_tests)
    fwrite(result, file.path(out_dir, paste0(trait, ".tsv")), sep = "\t")
  }

  for (trait in all_traits) {
    out_dir <- file.path(OUTPUT_DIR, model_name, "binary", "all")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    message("[", model_name, "] Binary (all): ", trait)
    result <- run_firth_trait(data, trait, genes, c(SEX_VAR, AGE_VAR, SMOKING_VAR, BMI_VAR), n_tests)
    fwrite(result, file.path(out_dir, paste0(trait, ".tsv")), sep = "\t")
  }

  male_data <- data[data[[SEX_VAR]] == 0, , drop = FALSE]
  for (trait in male_traits) {
    out_dir <- file.path(OUTPUT_DIR, model_name, "binary", "male")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    message("[", model_name, "] Binary (male): ", trait)
    result <- run_firth_trait(male_data, trait, genes, c(AGE_VAR, SMOKING_VAR, BMI_VAR), n_tests)
    fwrite(result, file.path(out_dir, paste0(trait, ".tsv")), sep = "\t")
  }

  female_data <- data[data[[SEX_VAR]] == 1, , drop = FALSE]
  for (trait in female_traits) {
    out_dir <- file.path(OUTPUT_DIR, model_name, "binary", "female")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    message("[", model_name, "] Binary (female): ", trait)
    result <- run_firth_trait(female_data, trait, genes, c(AGE_VAR, SMOKING_VAR, BMI_VAR), n_tests)
    fwrite(result, file.path(out_dir, paste0(trait, ".tsv")), sep = "\t")
  }
}

# =============================================================================
# 8. Run Korean- and GTEx-based TWAS
# =============================================================================

run_model_twas <- function(model_name, phenotypes) {
  message("============================================================")
  message("Running ", model_name, "-based TWAS")
  message("============================================================")

  grex <- read_predicted_expression(PREDICTED_EXPRESSION[[model_name]])
  r2_info <- get_r2_filtered_genes(MODEL_DATABASES[[model_name]])
  gene_info <- select_twas_genes(grex, r2_info)
  n_tests <- gene_info$n_model_genes

  message(
    "[", model_name, "] R2 > ", R2_THRESHOLD, ": ", gene_info$n_model_genes,
    " genes; available GReX: ", gene_info$n_available_genes,
    "; after zero-ratio filter: ", gene_info$n_analysis_genes,
    "; Bonferroni threshold: ", signif(0.05 / n_tests, 5)
  )

  model_dir <- file.path(OUTPUT_DIR, model_name)
  dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

  fwrite(
    data.frame(
      model = model_name,
      R2_threshold = R2_THRESHOLD,
      n_R2_filtered_models = gene_info$n_model_genes,
      n_models_in_predicted_expression = gene_info$n_available_genes,
      n_models_after_zero_ratio_filter = gene_info$n_analysis_genes,
      multiple_testing_n = n_tests,
      bonferroni_threshold = 0.05 / n_tests
    ),
    file.path(model_dir, "TWAS_model_summary.tsv"),
    sep = "\t"
  )
  fwrite(r2_info, file.path(model_dir, "R2_filtered_GReX_models.tsv"), sep = "\t")

  analysis_data <- merge_grex_with_phenotypes(
    phenotypes,
    grex[, gene_info$analysis_genes, drop = FALSE]
  )

  run_continuous_twas(analysis_data, gene_info$analysis_genes, n_tests, model_name)
  run_binary_twas(analysis_data, gene_info$analysis_genes, n_tests, model_name)
}

# =============================================================================
# 9. Main
# =============================================================================

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
phenotypes <- load_phenotypes()

args <- commandArgs(trailingOnly = TRUE)
models_to_run <- if (length(args) == 0 || tolower(args[1]) == "all") {
  c("Korean", "GTEx")
} else {
  args[1]
}

if (!all(models_to_run %in% names(MODEL_DATABASES))) {
  stop("Usage: Rscript 03_run_TWAS.R [all|Korean|GTEx]", call. = FALSE)
}

for (model_name in models_to_run) run_model_twas(model_name, phenotypes)
message("TWAS analysis completed.")
