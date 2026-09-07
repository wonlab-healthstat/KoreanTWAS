#!/usr/bin/env Rscript

###############################################################################
# 04_run_GIFT.R
#
# Conditional fine-mapping of Bonferroni-significant Korean-based TWAS loci
# using GIFT.
#
# Workflow:
#   1. Load Bonferroni-significant Korean TWAS associations from 03_run_TWAS.R.
#   2. Define the local candidate-gene set using EAS LD blocks.
#   3. Harmonize KoGES+GENIE alleles to the Asan+Chosun reference genotype.
#   4. Generate gene-specific PLINK files for Asan+Chosun and KoGES+GENIE.
#   5. Run GIFT_individual() for each significant gene-trait association.
#   6. Summarize the conditional result for the focal TWAS gene and calculate
#      regional and transcriptome-wide Bonferroni-adjusted P values.
#
# GIFT settings used in the study:
#   maxiter = 100
#   tol     = 1e-3
#   pleio   = 0
#   filter  = TRUE
#
# Local cis-genotype window:
#   +/- 100 kb around each candidate gene.
#
# Candidate genes:
#   GReX models with rho_avg_squared > 0.05 that overlap the LD block
#   containing the focal TWAS gene.
#
# A minimum of 10 common cis-SNPs between Asan+Chosun and KoGES+GENIE is
# required for a candidate gene.
#
# LD blocks:
#   GRCh38 East Asian LD blocks from:
#   https://github.com/jmacdon/LDblocks_GRCh38
#
# Expression preprocessing:
#   GIFT uses covariate-adjusted Asan+Chosun gene expression prepared before
#   running this script. In the study, expression was adjusted for sex, age,
#   three genotype PCs, and 60 PEER factors.
#
# Usage:
#   Rscript 04_run_GIFT.R prepare
#   Rscript 04_run_GIFT.R run
#   Rscript 04_run_GIFT.R summarize
#   Rscript 04_run_GIFT.R all
#
# Optional single-trait execution:
#   Rscript 04_run_GIFT.R prepare MCHC
#   Rscript 04_run_GIFT.R run MCHC
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(genio)
  library(GIFT)
})

# =============================================================================
# 1. Configuration
# =============================================================================

PLINK_BIN <- "/path/to/plink"

# Output from 03_run_TWAS.R for the Korean GReX model.
TWAS_DIR <- "/path/to/TWAS_results/Korean"
TWAS_MODEL_SUMMARY <- file.path(TWAS_DIR, "TWAS_model_summary.tsv")
R2_FILTERED_MODELS <- file.path(TWAS_DIR, "R2_filtered_GReX_models.tsv")

# KoGES+GENIE phenotype file. Trait columns are assumed to use abbreviations
# directly (e.g. MCHC, HDL, KD, LIP) and sample ID is IID.
PHENOTYPE_FILE <- "/path/to/KoGES_GENIE_phenotypes.csv"
PHENOTYPE_ID_COL <- "IID"
SEX_COL <- "SEX"

# Asan+Chosun GReX-reference data.
ASAN_GENOTYPE_FILE <- "/path/to/Asan_Chosun/Final_genotype.csv"
ASAN_SNP_LOCATION_FILE <- "/path/to/Asan_Chosun/Final_SNP_loc.txt"
ASAN_SNP_ANNOTATION_FILE <- "/path/to/Asan_Chosun/snp_annot.all.txt"
ASAN_EXPRESSION_FILE <- "/path/to/Asan_Chosun/expression_adjusted.txt"
ASAN_COVARIATE_FILE <- "/path/to/Asan_Chosun/Final_cvrt.csv"

# Gene annotation used for GReX/TWAS.
GENE_ANNOTATION_FILE <- "/path/to/Position_gene.txt"

# PrediXcan-ready KoGES+GENIE dosage files generated before this analysis:
#   chr1.txt.gz ... chr22.txt.gz
# Expected leading columns:
#   chr, SNP, pos, major, minor, maf
KOGES_GENOTYPE_DIR <- "/path/to/KoGES_GENIE/genotype"

# GRCh38 EAS LD blocks.
LD_BLOCK_FILE <- "/path/to/pyrho_EAS_LD_blocks.bed"

WORK_DIR <- "/path/to/GIFT"
PLINK_INPUT_DIR <- file.path(WORK_DIR, "input", "gene_plink")
GIFT_OUTPUT_DIR <- file.path(WORK_DIR, "output")
SUMMARY_DIR <- file.path(WORK_DIR, "summary")

CIS_WINDOW <- 100000
MIN_CIS_SNPS <- 10
N_CORES <- 1

# =============================================================================
# 2. Utility functions
# =============================================================================

check_file <- function(path) {
  if (!file.exists(path)) stop("Required file not found: ", path, call. = FALSE)
}

check_dir <- function(path) {
  if (!dir.exists(path)) stop("Required directory not found: ", path, call. = FALSE)
}

normalize_chr <- function(x) {
  as.character(str_remove(as.character(x), "^chr"))
}

load_korean_twas <- function(target_trait = NULL) {
  dirs <- c(file.path(TWAS_DIR, "continuous"), file.path(TWAS_DIR, "binary"))
  files <- unlist(lapply(dirs[dir.exists(dirs)], function(x) {
    list.files(x, pattern = "\\.tsv$", recursive = TRUE, full.names = TRUE)
  }))
  if (length(files) == 0) stop("No Korean TWAS result files found under: ", TWAS_DIR, call. = FALSE)

  out <- bind_rows(lapply(files, function(f) {
    x <- fread(f, data.table = FALSE, check.names = FALSE)
    if (!all(c("trait", "gene", "p") %in% colnames(x))) return(NULL)
    x$trait_type <- if (grepl(paste0(.Platform$file.sep, "continuous", .Platform$file.sep), f)) "CONTINUOUS" else "BINARY"
    x
  }))
  if (!is.null(target_trait)) out <- out[out$trait == target_trait, , drop = FALSE]

  if ("significant_bonferroni" %in% colnames(out)) {
    sig <- out[which(as.logical(out$significant_bonferroni)), , drop = FALSE]
  } else if ("bonferroni_threshold" %in% colnames(out)) {
    sig <- out[which(out$p < out$bonferroni_threshold), , drop = FALSE]
  } else {
    stop("TWAS results must contain significant_bonferroni or bonferroni_threshold.", call. = FALSE)
  }

  sig <- sig[!duplicated(paste(sig$trait, sig$gene, sep = "::")), , drop = FALSE]
  sig
}

load_twas_model_info <- function() {
  check_file(TWAS_MODEL_SUMMARY)
  check_file(R2_FILTERED_MODELS)

  model_summary <- fread(TWAS_MODEL_SUMMARY, data.table = FALSE)
  if (!"multiple_testing_n" %in% colnames(model_summary)) {
    stop("TWAS_model_summary.tsv must contain multiple_testing_n.", call. = FALSE)
  }

  r2_models <- fread(R2_FILTERED_MODELS, data.table = FALSE)
  if (!"gene" %in% colnames(r2_models)) stop("R2_filtered_GReX_models.tsv must contain gene.", call. = FALSE)

  list(
    n_transcriptome_tests = as.integer(model_summary$multiple_testing_n[1]),
    genes = unique(as.character(r2_models$gene))
  )
}

# =============================================================================
# 3. Load Asan+Chosun reference data
# =============================================================================

load_asan_reference <- function(grex_genes) {
  check_file(ASAN_GENOTYPE_FILE)
  check_file(ASAN_SNP_LOCATION_FILE)
  check_file(ASAN_SNP_ANNOTATION_FILE)
  check_file(ASAN_EXPRESSION_FILE)
  check_file(ASAN_COVARIATE_FILE)
  check_file(GENE_ANNOTATION_FILE)

  raw_geno <- fread(ASAN_GENOTYPE_FILE, data.table = FALSE, check.names = FALSE)
  variant_ids <- as.character(raw_geno[[1]])
  G <- t(data.matrix(raw_geno[, -1, drop = FALSE]))
  colnames(G) <- variant_ids
  sample_ids <- rownames(G)

  snp_loc <- fread(ASAN_SNP_LOCATION_FILE, data.table = FALSE, check.names = FALSE)
  snp_annot <- fread(ASAN_SNP_ANNOTATION_FILE, data.table = FALSE, check.names = FALSE)
  required_loc <- c("chr", "snpid", "pos")
  required_annot <- c("varID", "ref_vcf", "alt_vcf", "rsid")
  if (!all(required_loc %in% colnames(snp_loc))) stop("Unexpected Asan SNP-location format.", call. = FALSE)
  if (!all(required_annot %in% colnames(snp_annot))) stop("Unexpected Asan SNP-annotation format.", call. = FALSE)

  snp_loc$CHR <- normalize_chr(snp_loc$chr)
  idx <- match(snp_loc$snpid, snp_annot$varID)
  snp_loc$ref_vcf <- snp_annot$ref_vcf[idx]
  snp_loc$alt_vcf <- snp_annot$alt_vcf[idx]
  snp_loc$rsid <- snp_annot$rsid[idx]
  snp_loc$SNP_GIFT <- paste0(snp_loc$rsid, "_", snp_loc$alt_vcf)

  map <- snp_loc$SNP_GIFT[match(colnames(G), snp_loc$snpid)]
  keep <- !is.na(map) & map != "NA_NA" & !duplicated(map)
  G <- G[, keep, drop = FALSE]
  colnames(G) <- map[keep]
  storage.mode(G) <- "integer"

  snp_loc <- snp_loc[snp_loc$SNP_GIFT %in% colnames(G), , drop = FALSE]
  snp_loc <- snp_loc[!duplicated(snp_loc$SNP_GIFT), , drop = FALSE]

  expression <- fread(ASAN_EXPRESSION_FILE, data.table = FALSE, check.names = FALSE)
  colnames(expression)[1] <- "ID"
  expression$ID <- as.character(expression$ID)
  expression <- expression[match(sample_ids, expression$ID), , drop = FALSE]
  if (any(is.na(expression$ID))) stop("Some Asan+Chosun genotype samples are missing from expression data.", call. = FALSE)

  raw_cov <- fread(ASAN_COVARIATE_FILE, data.table = FALSE, check.names = FALSE)
  covariates <- as.data.frame(t(raw_cov[, -1, drop = FALSE]), check.names = FALSE)
  colnames(covariates) <- as.character(raw_cov[[1]])
  covariates$ID <- rownames(covariates)
  covariates <- covariates[match(sample_ids, covariates$ID), , drop = FALSE]
  if (!"sex" %in% colnames(covariates)) stop("Asan covariate file must contain sex.", call. = FALSE)
  asan_sex <- ifelse(as.numeric(covariates$sex) == 0, 1, 2)

  genes <- fread(GENE_ANNOTATION_FILE, data.table = FALSE, check.names = FALSE)
  required_gene <- c("chr", "gene_id", "gene_name", "start", "end", "gene_type")
  if (!all(required_gene %in% colnames(genes))) stop("Unexpected gene-annotation format.", call. = FALSE)

  genes$chr <- normalize_chr(genes$chr)
  genes <- genes[
    genes$chr %in% as.character(1:22) &
      genes$gene_type %in% c("protein_coding", "lincRNA") &
      genes$gene_id %in% grex_genes &
      genes$gene_id %in% colnames(expression),
    ,
    drop = FALSE
  ]

  bim <- data.frame(
    chr = snp_loc$CHR,
    id = snp_loc$SNP_GIFT,
    posg = 0,
    pos = snp_loc$pos,
    alt = snp_loc$alt_vcf,
    ref = snp_loc$ref_vcf,
    stringsAsFactors = FALSE
  )

  list(
    G = G,
    sample_ids = sample_ids,
    snp = snp_loc,
    bim = bim,
    expression = expression,
    sex = asan_sex,
    genes = genes
  )
}

# =============================================================================
# 4. Load and harmonize KoGES+GENIE genotype
# =============================================================================

load_koges_chr <- function(chr, asan_snp, phenotype_ids) {
  f <- file.path(KOGES_GENOTYPE_DIR, paste0("chr", chr, ".txt.gz"))
  check_file(f)

  geno <- fread(f, data.table = FALSE, check.names = FALSE)
  required <- c("chr", "SNP", "pos", "major", "minor", "maf")
  if (!all(required %in% colnames(geno))) stop("Unexpected KoGES dosage format: ", f, call. = FALSE)

  target_alt <- asan_snp$alt_vcf[match(geno$SNP, asan_snp$rsid)]
  target_ref <- asan_snp$ref_vcf[match(geno$SNP, asan_snp$rsid)]
  matched <- !is.na(target_alt) & !is.na(target_ref)
  aligned <- matched & geno$minor == target_alt & geno$major == target_ref
  flipped <- matched & geno$major == target_alt & geno$minor == target_ref
  keep <- aligned | flipped

  sample_cols <- phenotype_ids[phenotype_ids %in% colnames(geno)]
  if (length(sample_cols) == 0) stop("No KoGES phenotype samples found in genotype file: ", f, call. = FALSE)

  geno <- geno[keep, c(required, sample_cols), drop = FALSE]
  flip_idx <- which(flipped[keep])
  if (length(flip_idx) > 0) {
    old_major <- geno$major[flip_idx]
    geno$major[flip_idx] <- geno$minor[flip_idx]
    geno$minor[flip_idx] <- old_major

    G <- data.matrix(geno[, sample_cols, drop = FALSE])
    G[flip_idx, ] <- 2 - G[flip_idx, ]
    geno[, sample_cols] <- G
  }

  geno$SNP_GIFT <- paste0(geno$SNP, "_", geno$minor)
  geno <- geno[!duplicated(geno$SNP_GIFT), , drop = FALSE]
  geno
}

# =============================================================================
# 5. Prepare gene-specific PLINK inputs
# =============================================================================

prepare_gift_inputs <- function(target_trait = NULL) {
  check_file(PHENOTYPE_FILE)
  check_file(LD_BLOCK_FILE)
  check_dir(KOGES_GENOTYPE_DIR)

  sig_twas <- load_korean_twas(target_trait)
  if (nrow(sig_twas) == 0) {
    message("No Bonferroni-significant Korean TWAS associations to prepare.")
    return(invisible(NULL))
  }

  model_info <- load_twas_model_info()
  ref <- load_asan_reference(model_info$genes)

  pheno <- fread(PHENOTYPE_FILE, data.table = FALSE, check.names = FALSE)
  if (!all(c(PHENOTYPE_ID_COL, SEX_COL) %in% colnames(pheno))) {
    stop("Phenotype file must contain ", PHENOTYPE_ID_COL, " and ", SEX_COL, ".", call. = FALSE)
  }
  pheno[[PHENOTYPE_ID_COL]] <- as.character(pheno[[PHENOTYPE_ID_COL]])

  ld <- fread(LD_BLOCK_FILE, data.table = FALSE, check.names = FALSE)
  colnames(ld)[1:3] <- c("chr", "start", "end")
  ld$chr <- normalize_chr(ld$chr)

  gene_dt <- as.data.table(ref$genes)
  snp_dt <- as.data.table(ref$snp)
  ld_dt <- as.data.table(ld)
  setkey(gene_dt, chr, start, end)
  setkey(snp_dt, CHR, pos)
  setkey(ld_dt, chr, start, end)

  dir.create(file.path(PLINK_INPUT_DIR, "Asan_Chosun"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(PLINK_INPUT_DIR, "KoGES_GENIE"), recursive = TRUE, showWarnings = FALSE)

  for (chr_now in unique(gene_dt$chr[gene_dt$gene_id %in% sig_twas$gene])) {
    chr_targets <- sig_twas[sig_twas$gene %in% gene_dt$gene_id[gene_dt$chr == chr_now], , drop = FALSE]
    if (nrow(chr_targets) == 0) next

    message("Preparing chromosome ", chr_now)
    koges_chr <- load_koges_chr(chr_now, ref$snp, pheno[[PHENOTYPE_ID_COL]])
    koges_sample_ids <- pheno[[PHENOTYPE_ID_COL]][pheno[[PHENOTYPE_ID_COL]] %in% colnames(koges_chr)]
    G2 <- data.matrix(koges_chr[, koges_sample_ids, drop = FALSE])
    rownames(G2) <- koges_chr$SNP_GIFT
    colnames(G2) <- koges_sample_ids
    storage.mode(G2) <- "integer"

    bim2 <- data.frame(
      chr = normalize_chr(koges_chr$chr),
      id = koges_chr$SNP_GIFT,
      posg = 0,
      pos = koges_chr$pos,
      alt = koges_chr$minor,
      ref = koges_chr$major,
      stringsAsFactors = FALSE
    )

    for (r in seq_len(nrow(chr_targets))) {
      trait <- chr_targets$trait[r]
      focal_gene <- chr_targets$gene[r]
      if (!trait %in% colnames(pheno)) {
        message("Skipping ", trait, ": phenotype column not found.")
        next
      }

      target <- gene_dt[gene_id == focal_gene]
      if (nrow(target) == 0) next
      target <- target[1]

      blocks <- ld_dt[
        chr == target$chr &
          start <= target$end &
          end >= target$start
      ]
      if (nrow(blocks) == 0) {
        message(trait, " / ", focal_gene, ": no overlapping LD block.")
        next
      }

      region_start <- min(blocks$start)
      region_end <- max(blocks$end)
      candidates <- gene_dt[
        chr == target$chr &
          start <= region_end &
          end >= region_start
      ]
      if (nrow(candidates) == 0) next

      asan_out <- file.path(PLINK_INPUT_DIR, "Asan_Chosun", trait, focal_gene)
      koges_out <- file.path(PLINK_INPUT_DIR, "KoGES_GENIE", trait, focal_gene)
      dir.create(asan_out, recursive = TRUE, showWarnings = FALSE)
      dir.create(koges_out, recursive = TRUE, showWarnings = FALSE)

      for (j in seq_len(nrow(candidates))) {
        gene <- candidates$gene_id[j]
        cis_start <- max(candidates$start[j] - CIS_WINDOW, 1)
        cis_end <- candidates$end[j] + CIS_WINDOW

        local_snps <- snp_dt[
          CHR == target$chr &
            pos >= cis_start &
            pos <= cis_end &
            SNP_GIFT %in% rownames(G2)
        ]
        local_snps <- local_snps[!duplicated(SNP_GIFT)]
        if (nrow(local_snps) < MIN_CIS_SNPS) next

        snps <- local_snps$SNP_GIFT
        snps <- snps[snps %in% colnames(ref$G) & snps %in% rownames(G2)]
        if (length(snps) < MIN_CIS_SNPS) next

        asan_bim <- ref$bim[match(snps, ref$bim$id), , drop = FALSE]
        asan_X <- t(ref$G[, snps, drop = FALSE])
        asan_expression <- ref$expression[[gene]]
        asan_fam <- data.frame(
          fam = ref$sample_ids,
          id = ref$sample_ids,
          pat = 0,
          mat = 0,
          sex = ref$sex,
          pheno = asan_expression,
          stringsAsFactors = FALSE
        )
        genio::write_plink(
          file = file.path(asan_out, gene),
          X = asan_X,
          bim = asan_bim,
          fam = asan_fam
        )

        koges_bim <- bim2[match(snps, bim2$id), , drop = FALSE]
        pheno_idx <- match(koges_sample_ids, pheno[[PHENOTYPE_ID_COL]])
        koges_fam <- data.frame(
          fam = koges_sample_ids,
          id = koges_sample_ids,
          pat = 0,
          mat = 0,
          sex = pheno[[SEX_COL]][pheno_idx],
          pheno = pheno[[trait]][pheno_idx],
          stringsAsFactors = FALSE
        )
        genio::write_plink(
          file = file.path(koges_out, gene),
          X = G2[snps, , drop = FALSE],
          bim = koges_bim,
          fam = koges_fam
        )
      }

      n_candidates <- length(list.files(asan_out, pattern = "\\.bed$"))
      message(trait, " / ", focal_gene, ": ", n_candidates, " local GIFT candidates prepared.")
    }

    rm(koges_chr, G2, bim2)
    gc()
  }

  invisible(sig_twas)
}

# =============================================================================
# 6. Run GIFT
# =============================================================================

load_expression_for_gift <- function() {
  x <- fread(ASAN_EXPRESSION_FILE, data.table = FALSE, check.names = FALSE)
  colnames(x)[1] <- "ID"
  x$ID <- paste0(x$ID, "_", x$ID)
  rownames(x) <- x$ID
  x
}

load_phenotype_for_gift <- function() {
  y <- fread(PHENOTYPE_FILE, data.table = FALSE, check.names = FALSE)
  y[[PHENOTYPE_ID_COL]] <- as.character(y[[PHENOTYPE_ID_COL]])
  y$GIFT_ID <- paste0(y[[PHENOTYPE_ID_COL]], "_", y[[PHENOTYPE_ID_COL]])
  rownames(y) <- y$GIFT_ID
  y
}

run_one_gift <- function(trait, focal_gene, Gene, Pheno) {
  asan_dir <- file.path(PLINK_INPUT_DIR, "Asan_Chosun", trait, focal_gene)
  koges_dir <- file.path(PLINK_INPUT_DIR, "KoGES_GENIE", trait, focal_gene)
  out_dir <- file.path(GIFT_OUTPUT_DIR, trait, focal_gene)
  out_file <- file.path(out_dir, "GIFT_results.RData")

  if (file.exists(out_file)) {
    return(data.frame(trait = trait, gene = focal_gene, status = "ALREADY_DONE"))
  }
  if (!dir.exists(asan_dir) || !dir.exists(koges_dir)) {
    return(data.frame(trait = trait, gene = focal_gene, status = "NO_PLINK_INPUT"))
  }

  n_asan <- length(list.files(asan_dir, pattern = "\\.bed$"))
  n_koges <- length(list.files(koges_dir, pattern = "\\.bed$"))
  if (n_asan < 2 || n_koges < 2) {
    return(data.frame(trait = trait, gene = focal_gene, status = "LESS_THAN_2_CANDIDATES"))
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  GReX_convert <- tryCatch(
    GIFT::pre_process_individual(asan_dir, PLINK_BIN),
    error = function(e) NULL
  )
  TWAS_convert <- tryCatch(
    GIFT::pre_process_individual(koges_dir, PLINK_BIN),
    error = function(e) NULL
  )
  if (is.null(GReX_convert) || is.null(TWAS_convert)) {
    return(data.frame(trait = trait, gene = focal_gene, status = "PREPROCESS_ERROR"))
  }

  genes <- GReX_convert$gene
  Zx <- GReX_convert$Z
  pindex <- GReX_convert$pindex
  Zy <- TWAS_convert$Z

  missing_genes <- setdiff(genes, colnames(Gene))
  if (length(missing_genes) > 0) {
    return(data.frame(trait = trait, gene = focal_gene, status = "EXPRESSION_GENE_MISSING"))
  }
  if (!trait %in% colnames(Pheno)) {
    return(data.frame(trait = trait, gene = focal_gene, status = "PHENOTYPE_MISSING"))
  }

  X <- Gene[match(rownames(Zx), rownames(Gene)), genes, drop = FALSE]
  Y <- Pheno[match(rownames(Zy), rownames(Pheno)), trait, drop = FALSE]
  X <- as.matrix(X)
  Y <- as.matrix(Y)

  result <- tryCatch(
    GIFT::GIFT_individual(
      X, Y, Zx, Zy, genes, pindex,
      maxiter = 100,
      tol = 1e-3,
      pleio = 0,
      ncores = 1,
      filter = TRUE
    ),
    error = function(e) e
  )

  if (inherits(result, "error")) {
    return(data.frame(trait = trait, gene = focal_gene, status = paste0("GIFT_ERROR: ", conditionMessage(result))))
  }

  GIFT_results_gene <- list(filterT = result)
  save(GIFT_results_gene, file = out_file)

  data.frame(
    trait = trait,
    gene = focal_gene,
    status = "DONE",
    n_candidates = length(genes),
    stringsAsFactors = FALSE
  )
}

run_gift <- function(target_trait = NULL) {
  sig_twas <- load_korean_twas(target_trait)
  if (nrow(sig_twas) == 0) {
    message("No Bonferroni-significant Korean TWAS associations to run.")
    return(invisible(NULL))
  }

  Gene <- load_expression_for_gift()
  Pheno <- load_phenotype_for_gift()

  jobs <- unique(sig_twas[, c("trait", "gene"), drop = FALSE])
  worker <- function(i) run_one_gift(jobs$trait[i], jobs$gene[i], Gene, Pheno)

  if (N_CORES > 1 && .Platform$OS.type != "windows") {
    statuses <- parallel::mclapply(seq_len(nrow(jobs)), worker, mc.cores = N_CORES)
  } else {
    statuses <- lapply(seq_len(nrow(jobs)), worker)
  }

  statuses <- bind_rows(statuses)
  dir.create(SUMMARY_DIR, recursive = TRUE, showWarnings = FALSE)
  fwrite(statuses, file.path(SUMMARY_DIR, "GIFT_run_status.tsv"), sep = "\t")
  print(table(statuses$status))
  invisible(statuses)
}

# =============================================================================
# 7. Summarize GIFT results
# =============================================================================

extract_focal_result <- function(res, focal_gene) {
  if (is.null(res) || nrow(res) == 0 || !all(c("gene", "causal_effect", "p") %in% colnames(res))) {
    return(data.frame(gene_in_table = NA_character_, causal_effect = NA_real_, p = NA_real_))
  }

  idx <- which(res$gene == focal_gene)
  if (length(idx) == 0) {
    return(data.frame(gene_in_table = NA_character_, causal_effect = NA_real_, p = NA_real_))
  }

  data.frame(
    gene_in_table = as.character(res$gene[idx[1]]),
    causal_effect = as.numeric(res$causal_effect[idx[1]]),
    p = as.numeric(res$p[idx[1]])
  )
}

summarize_gift <- function(target_trait = NULL) {
  sig_twas <- load_korean_twas(target_trait)
  model_info <- load_twas_model_info()

  rows <- lapply(seq_len(nrow(sig_twas)), function(i) {
    trait <- sig_twas$trait[i]
    focal_gene <- sig_twas$gene[i]
    f <- file.path(GIFT_OUTPUT_DIR, trait, focal_gene, "GIFT_results.RData")

    base <- data.frame(
      trait = trait,
      trait_type = sig_twas$trait_type[i],
      target_gene = focal_gene,
      twas_beta = if ("beta" %in% colnames(sig_twas)) sig_twas$beta[i] else NA_real_,
      twas_p = sig_twas$p[i],
      stringsAsFactors = FALSE
    )

    if (!file.exists(f)) {
      return(cbind(
        base,
        gene_in_table = NA_character_,
        causal_effect = NA_real_,
        gift_p = NA_real_,
        n_candidates = NA_integer_,
        lowest_p = NA
      ))
    }

    e <- new.env(parent = emptyenv())
    load(f, envir = e)
    if (!exists("GIFT_results_gene", envir = e)) {
      return(cbind(
        base,
        gene_in_table = NA_character_,
        causal_effect = NA_real_,
        gift_p = NA_real_,
        n_candidates = NA_integer_,
        lowest_p = NA
      ))
    }

    obj <- get("GIFT_results_gene", envir = e)
    res <- if (is.list(obj) && "filterT" %in% names(obj)) as.data.frame(obj$filterT) else NULL
    picked <- extract_focal_result(res, focal_gene)

    min_p <- if (!is.null(res) && nrow(res) > 0 && !all(is.na(res$p))) min(res$p, na.rm = TRUE) else NA_real_
    lowest <- if (!is.na(picked$p) && !is.na(min_p)) isTRUE(all.equal(picked$p, min_p)) else NA

    cbind(
      base,
      gene_in_table = picked$gene_in_table,
      causal_effect = picked$causal_effect,
      gift_p = picked$p,
      n_candidates = if (!is.null(res)) nrow(res) else NA_integer_,
      lowest_p = lowest
    )
  })

  summary <- bind_rows(rows)
  summary$gene_name <- NA_character_

  if (file.exists(GENE_ANNOTATION_FILE)) {
    annot <- fread(GENE_ANNOTATION_FILE, data.table = FALSE, check.names = FALSE)
    if (all(c("gene_id", "gene_name") %in% colnames(annot))) {
      summary$gene_name <- annot$gene_name[match(summary$target_gene, annot$gene_id)]
    }
  }

  summary$p_regional <- pmin(summary$gift_p * summary$n_candidates, 1)
  summary$p_transcriptome <- pmin(summary$gift_p * model_info$n_transcriptome_tests, 1)
  summary$regional_significant <- summary$p_regional < 0.05
  summary$transcriptome_significant <- summary$p_transcriptome < 0.05
  summary <- summary %>%
    select(
      trait, trait_type, target_gene, gene_name,
      twas_beta, twas_p, causal_effect, gift_p,
      n_candidates, lowest_p,
      p_regional, regional_significant,
      p_transcriptome, transcriptome_significant
    ) %>%
    arrange(trait_type, trait, target_gene)

  dir.create(SUMMARY_DIR, recursive = TRUE, showWarnings = FALSE)
  fwrite(summary, file.path(SUMMARY_DIR, "GIFT_finemap_summary_ALL.tsv"), sep = "\t")
  fwrite(summary[summary$trait_type == "BINARY", ], file.path(SUMMARY_DIR, "GIFT_finemap_summary_BINARY.tsv"), sep = "\t")
  fwrite(summary[summary$trait_type == "CONTINUOUS", ], file.path(SUMMARY_DIR, "GIFT_finemap_summary_CONTINUOUS.tsv"), sep = "\t")

  message(
    "GIFT summary: ", nrow(summary), " focal associations; ",
    sum(summary$regional_significant, na.rm = TRUE), " regionally significant; ",
    sum(summary$transcriptome_significant, na.rm = TRUE), " transcriptome-wide significant."
  )
  invisible(summary)
}

# =============================================================================
# 8. Command-line interface
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) >= 1) tolower(args[1]) else "all"
target_trait <- if (length(args) >= 2) args[2] else NULL

dir.create(WORK_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(GIFT_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SUMMARY_DIR, recursive = TRUE, showWarnings = FALSE)

if (mode == "prepare") {
  prepare_gift_inputs(target_trait)
} else if (mode == "run") {
  run_gift(target_trait)
} else if (mode == "summarize") {
  summarize_gift(target_trait)
} else if (mode == "all") {
  prepare_gift_inputs(target_trait)
  run_gift(target_trait)
  summarize_gift(target_trait)
} else {
  stop("Usage: Rscript 04_run_GIFT.R [prepare|run|summarize|all] [optional_trait]", call. = FALSE)
}
