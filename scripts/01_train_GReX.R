#!/usr/bin/env Rscript

###############################################################################
# 01_train_GReX.R
#
# Korean whole-blood GReX model training using PredictDB-Tutorial.
#
# Reference implementation:
# https://github.com/hakyimlab/PredictDB-Tutorial
#
# The original model-training code is not redistributed here. This script
# documents how the PredictDB-Tutorial training function was executed with the
# study-specific preprocessed inputs and parameters.
#
# Required inputs:
#   <input_dir>/gene_annot.parsed.txt
#   <input_dir>/genotype/snp_annot.chr{1..22}.txt
#   <input_dir>/genotype/genotype.chr{1..22}.txt
#   <input_dir>/expression/transformed_expression.txt
#   <input_dir>/PEER/covariates.txt
#
# Training settings:
#   MAF >= 0.01
#   cis-window = +/- 1 Mb
#   Elastic-Net alpha = 0.5
#   inner CV = 10 folds
#   outer nested CV = 5 folds
#
# Usage:
#   Rscript 01_train_GReX.R <chr|all> <input_dir> <output_dir> <PredictDB-Tutorial_dir> [seed]
#
# Examples:
#   Rscript 01_train_GReX.R 1 ./input ./model_output /path/to/PredictDB-Tutorial
#   Rscript 01_train_GReX.R all ./input ./model_output /path/to/PredictDB-Tutorial
###############################################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Usage: Rscript 01_train_GReX.R <chr|all> <input_dir> <output_dir> <PredictDB-Tutorial_dir> [seed]")

chr_arg <- args[1]
input_dir <- normalizePath(args[2], mustWork = TRUE)
output_dir <- args[3]
predictdb_dir <- normalizePath(args[4], mustWork = TRUE)
seed <- if (length(args) >= 5) as.numeric(args[5]) else NA

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)
for (d in c("code", "summary", "weights", "covariances")) dir.create(file.path(output_dir, d), recursive = TRUE, showWarnings = FALSE)

training_code <- file.path(predictdb_dir, "code", "gtex_v7_nested_cv_elnet.R")
if (!file.exists(training_code)) stop("PredictDB training code not found: ", training_code)
source(training_code)

gene_annot_file <- file.path(input_dir, "gene_annot.parsed.txt")
expression_file <- file.path(input_dir, "expression", "transformed_expression.txt")
covariates_file <- file.path(input_dir, "PEER", "covariates.txt")
for (f in c(gene_annot_file, expression_file, covariates_file)) if (!file.exists(f)) stop("Required input not found: ", f)

chromosomes <- if (tolower(chr_arg) == "all") 1:22 else as.numeric(chr_arg)
if (any(is.na(chromosomes)) || any(!chromosomes %in% 1:22)) stop("<chr> must be 1-22 or 'all'.")

old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(file.path(output_dir, "code"))

for (chrom in chromosomes) {
  snp_annot_file <- file.path(input_dir, "genotype", paste0("snp_annot.chr", chrom, ".txt"))
  genotype_file <- file.path(input_dir, "genotype", paste0("genotype.chr", chrom, ".txt"))
  for (f in c(snp_annot_file, genotype_file)) if (!file.exists(f)) stop("Required input not found: ", f)

  message("Training chromosome ", chrom, "...")
  main(
    snp_annot_file = snp_annot_file,
    gene_annot_file = gene_annot_file,
    genotype_file = genotype_file,
    expression_file = expression_file,
    covariates_file = covariates_file,
    chrom = chrom,
    prefix = "Model_training",
    maf = 0.01,
    n_folds = 10,
    n_train_test_folds = 5,
    seed = seed,
    cis_window = 1000000,
    alpha = 0.5,
    null_testing = FALSE
  )
}

message("GReX model training completed.")

###############################################################################
# PredictDB database construction
#
# Database construction and filtering were performed following the
# "Make a database" and "Filter the database" sections of PredictDB-Tutorial:
# https://github.com/hakyimlab/PredictDB-Tutorial
#
# The chromosome-specific outputs generated above were used as input:
#   summary/Model_training_chr*_model_summaries.txt
#   summary/Model_training_chr*_summary.txt
#   weights/Model_training_chr*_weights.txt
#
# The tutorial database procedure was applied with study-specific changes:
#   - population in sample_info: "KOR" instead of the tutorial example "EUR"
#   - database paths and filenames were changed for the Korean whole-blood model
#   - the tutorial filtering criterion was retained:
#       zscore_pval < 0.05 and rho_avg > 0.10
#
# The database-construction code itself is not duplicated in this repository
# because it is provided directly in the PredictDB-Tutorial.
###############################################################################
