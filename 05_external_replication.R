#!/usr/bin/env Rscript

###############################################################################
# 05_external_replication.R
#
# External replication of Korean GReX-based TWAS using BBJ and CKB GWAS
# summary statistics with S-PrediXcan.
#
# GWAS summary statistics from BioBank Japan (BBJ) and the China Kadoorie Biobank (CKB), used for external replication, 
# are publicly available through their respective PheWeb repositories (https://pheweb.jp/downloads, and (https://pheweb.ckbiobank.org/, respectively). 
# S-PrediXcan implementation:
#   MetaXcan
#   https://github.com/hakyimlab/MetaXcan
#
# Specifically:
#   MetaXcan/software/SPrediXcan.py
#
# This script does not redistribute MetaXcan. It documents the study-specific
# model, covariance, GWAS input schemas, and execution settings.
#
# Analyses:
#   1. Korean GReX model -> BBJ continuous traits
#   2. Korean GReX model -> BBJ binary traits
#   3. Korean GReX model -> CKB continuous and binary traits
#
# Usage:
#   Rscript 05_external_replication.R all
#   Rscript 05_external_replication.R BBJ_continuous
#   Rscript 05_external_replication.R BBJ_binary
#   Rscript 05_external_replication.R CKB
#   Rscript 05_external_replication.R collect
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

# =============================================================================
# 1. Configuration
# =============================================================================

PYTHON_BIN <- "python3"
SPREDIXCAN_SCRIPT <- "/path/to/MetaXcan/software/SPrediXcan.py"

# Korean whole-blood GReX model and covariance generated during model training.
MODEL_DB <- "/path/to/Korean_Whole_Blood_models_filtered_signif.db"
COVARIANCE_FILE <- "/path/to/Model_training_All_covariances.txt.gz"

# External GWAS summary statistics.
BBJ_CONTINUOUS_DIR <- "/path/to/BBJ/continuous"
BBJ_BINARY_DIR <- "/path/to/BBJ/binary"
CKB_DIR <- "/path/to/CKB"

OUTPUT_DIR <- "/path/to/external_replication"
VERBOSITY <- "10"

# BBJ binary summary statistics retained in the original analysis.
BBJ_BINARY_PATTERN <- "\\.AF\\.|\\.CP\\.|\\.GP\\.|\\.InH\\.|\\.OP\\."

# =============================================================================
# 2. GWAS column schemas
# =============================================================================

SCHEMAS <- list(
  BBJ_continuous = list(
    snp = "SNPID_new",
    effect_allele = "ALLELE1",
    non_effect_allele = "ALLELE0",
    chromosome = "CHR_hg38",
    position = "POS_hg38",
    frequency = "A1FREQ",
    beta = "BETA",
    se = "SE",
    p = "P_BOLT_LMM_INF"
  ),
  BBJ_binary = list(
    snp = "SNPID_new",
    effect_allele = "Allele2",
    non_effect_allele = "Allele1",
    chromosome = "CHR_hg38",
    position = "POS_hg38",
    frequency = "AF_Allele2",
    beta = "BETA",
    se = "SE",
    p = "p.value"
  ),
  CKB = list(
    snp = "SNP",
    effect_allele = "effect_allele",
    non_effect_allele = "other_allele",
    chromosome = "CHR_hg38",
    position = "POS_hg38",
    frequency = "effect_allele_frequency",
    beta = "beta",
    se = "standard_error",
    p = "p_value"
  )
)

# =============================================================================
# 3. Utility functions
# =============================================================================

check_file <- function(path) {
  if (!file.exists(path)) stop("Required file not found: ", path, call. = FALSE)
}

check_dir <- function(path) {
  if (!dir.exists(path)) stop("Required directory not found: ", path, call. = FALSE)
}

run_spredixcan <- function(gwas_folder, gwas_pattern, schema, output_file) {
  check_file(SPREDIXCAN_SCRIPT)
  check_file(MODEL_DB)
  check_file(COVARIANCE_FILE)
  check_dir(gwas_folder)

  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(output_file)) {
    message("Existing result found; skipping: ", output_file)
    return(invisible(output_file))
  }

  args <- c(
    SPREDIXCAN_SCRIPT,
    "--model_db_path", MODEL_DB,
    "--covariance", COVARIANCE_FILE,
    "--gwas_folder", gwas_folder,
    "--gwas_file_pattern", gwas_pattern,
    "--snp_column", schema$snp,
    "--effect_allele_column", schema$effect_allele,
    "--non_effect_allele_column", schema$non_effect_allele,
    "--chromosome_column", schema$chromosome,
    "--position_column", schema$position,
    "--freq_column", schema$frequency,
    "--beta_column", schema$beta,
    "--se_column", schema$se,
    "--pvalue_column", schema$p,
    "--verbosity", VERBOSITY,
    "--additional_output",
    "--output_file", output_file
  )

  status <- system2(PYTHON_BIN, args = args)
  if (!identical(status, 0L)) stop("S-PrediXcan failed: ", output_file, call. = FALSE)
  invisible(output_file)
}

extract_bbj_trait <- function(path) {
  parts <- str_split(basename(path), "\\.")[[1]]
  if (length(parts) < 4) stop("Unexpected BBJ filename: ", basename(path), call. = FALSE)
  parts[4]
}

extract_ckb_trait <- function(path) {
  str_split(basename(path), "\\.")[[1]][1]
}

escape_regex <- function(x) {
  str_replace_all(x, "([][{}()+*^$|\\\\?.])", "\\\\\\1")
}

# =============================================================================
# 4. BBJ continuous traits
# =============================================================================

run_bbj_continuous <- function() {
  check_dir(BBJ_CONTINUOUS_DIR)
  files <- list.files(
    BBJ_CONTINUOUS_DIR,
    pattern = "\\.auto\\.txt$",
    full.names = TRUE,
    recursive = TRUE
  )
  if (length(files) == 0) stop("No BBJ continuous GWAS files found.", call. = FALSE)

  traits <- vapply(files, extract_bbj_trait, character(1))
  out_dir <- file.path(OUTPUT_DIR, "BBJ", "continuous")

  for (i in seq_along(files)) {
    trait <- traits[i]
    message("[BBJ continuous] ", trait, " ", i, "/", length(files))
    run_spredixcan(
      gwas_folder = dirname(files[i]),
      gwas_pattern = ".*auto\\.txt$",
      schema = SCHEMAS$BBJ_continuous,
      output_file = file.path(out_dir, paste0(trait, ".csv"))
    )
  }
}

# =============================================================================
# 5. BBJ binary traits
# =============================================================================

run_bbj_binary <- function() {
  check_dir(BBJ_BINARY_DIR)
  files <- list.files(
    BBJ_BINARY_DIR,
    pattern = "\\.auto\\.txt$",
    full.names = TRUE,
    recursive = TRUE
  )
  files <- files[str_detect(files, BBJ_BINARY_PATTERN)]
  if (length(files) == 0) stop("No selected BBJ binary GWAS files found.", call. = FALSE)

  traits <- vapply(files, extract_bbj_trait, character(1))
  out_dir <- file.path(OUTPUT_DIR, "BBJ", "binary")

  for (i in seq_along(files)) {
    trait <- traits[i]
    message("[BBJ binary] ", trait, " ", i, "/", length(files))
    run_spredixcan(
      gwas_folder = dirname(files[i]),
      gwas_pattern = ".*auto\\.txt$",
      schema = SCHEMAS$BBJ_binary,
      output_file = file.path(out_dir, paste0(trait, ".csv"))
    )
  }
}

# =============================================================================
# 6. CKB continuous and binary traits
# =============================================================================

run_ckb <- function() {
  check_dir(CKB_DIR)
  files <- list.files(
    CKB_DIR,
    pattern = "\\.tsv\\.gz$",
    full.names = TRUE,
    recursive = TRUE
  )
  files <- files[!str_detect(basename(files), "chr")]
  if (length(files) == 0) stop("No CKB GWAS files found.", call. = FALSE)

  traits <- vapply(files, extract_ckb_trait, character(1))
  out_dir <- file.path(OUTPUT_DIR, "CKB")

  for (i in seq_along(files)) {
    trait <- traits[i]
    message("[CKB] ", trait, " ", i, "/", length(files))
    run_spredixcan(
      gwas_folder = dirname(files[i]),
      gwas_pattern = paste0("^", escape_regex(basename(files[i])), "$"),
      schema = SCHEMAS$CKB,
      output_file = file.path(out_dir, paste0(trait, ".csv"))
    )
  }
}

# =============================================================================
# 7. Combine S-PrediXcan outputs
# =============================================================================

collect_results <- function() {
  result_files <- list.files(
    OUTPUT_DIR,
    pattern = "\\.csv$",
    full.names = TRUE,
    recursive = TRUE
  )
  if (length(result_files) == 0) {
    message("No S-PrediXcan result files found.")
    return(invisible(NULL))
  }

  rows <- lapply(result_files, function(f) {
    x <- fread(f, data.table = FALSE, check.names = FALSE)
    normalized <- normalizePath(f, winslash = "/", mustWork = TRUE)
    out_root <- normalizePath(OUTPUT_DIR, winslash = "/", mustWork = TRUE)
    rel <- sub(paste0("^", escape_regex(out_root), "/"), "", normalized)
    parts <- str_split(rel, "/")[[1]]

    if (parts[1] == "BBJ") {
      x$cohort <- "BBJ"
      x$trait_type <- parts[2]
    } else {
      x$cohort <- "CKB"
      x$trait_type <- NA_character_
    }

    x$trait <- str_remove(basename(f), "\\.csv$")
    x
  })

  combined <- as.data.frame(rbindlist(rows, fill = TRUE))
  out_file <- file.path(OUTPUT_DIR, "S_PrediXcan_external_replication_ALL.tsv")
  fwrite(combined, out_file, sep = "\t")
  message("Combined results saved: ", out_file)
  invisible(combined)
}

# =============================================================================
# 8. Command-line interface
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) == 0) "all" else tolower(args[1])

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (mode == "all") {
  run_bbj_continuous()
  run_bbj_binary()
  run_ckb()
} else if (mode == "bbj_continuous") {
  run_bbj_continuous()
} else if (mode == "bbj_binary") {
  run_bbj_binary()
} else if (mode == "ckb") {
  run_ckb()
} else if (mode == "collect") {
  collect_results()
} else {
  stop(
    "Usage: Rscript 05_external_replication.R [all|BBJ_continuous|BBJ_binary|CKB|collect]",
    call. = FALSE
  )
}
