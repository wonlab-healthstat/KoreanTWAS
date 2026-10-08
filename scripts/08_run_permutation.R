#!/usr/bin/env Rscript

###############################################################################
# 08_run_permutation.R | Korean GReX-based TWAS permutation analysis
#
# Modes (run separately):
#   null        - all 1,630 TWAS genes x 81 traits; 1,000 permutations for
#                 continuous traits, 100 for binary traits. Summarize null
#                 P-value distributions and genomic inflation (lambda_GC).
#   selected348 - preserve the manifest of 348 Bonferroni-significant
#                 gene-trait associations; 100,000 ADDITIONAL permutations
#                 per pair (indices 1001:101000), with Monte Carlo P-values.
#
# Same permutation schemes as the original study code:
#   Continuous: Freedman-Lane residual permutations from covariate-only lm;
#               residuals shuffled within SEX x SMOKE strata. Compare |t|.
#   Binary: stratified label permutations, with adaptive strata selection
#           across SEX/SMOKE/age-quantile/BMI-quantile combinations.
#           Compare the PLR chi-square statistic from Firth logistic regression.
#   Selected-pair Monte Carlo P = (K + 1) / (B_valid + 1).
#   No generalized Pareto tail extrapolation is performed.
#
# Input can be a combined analysis-ready file, as in the original permutation
# script, OR separate phenotype + predicted-expression files from 03_run_TWAS.R.
# Trait columns must use abbreviations (e.g. BMI, HTN, SMOKE); prediction
# columns must use ENSG gene IDs. No genotype preparation is done here.
#
# Required selected-pair manifest: 348 rows with columns trait_id and Gene,
# optional analysis_group; an original B_-coded trait_id can be converted below.
# The script NEVER replaces the selected 348 pairs with a new list.
#
###############################################################################

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
suppressPackageStartupMessages({ library(data.table); library(logistf) })
setDTthreads(1L)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

# 1. Configuration ------------------------------------------------------------
# Set DATA_FILE to a combined, analysis-ready table to use the original format.
# Alternatively, leave DATA_FILE as NA and supply the two separate files.
DATA_FILE <- NA_character_
PHENOTYPE_FILE <- "/path/to/KoGES_GENIE_phenotypes.csv"
KOREAN_GREX_FILE <- "/path/to/KoGES_GENIE_Korean_predicted_expression.txt"
PHENOTYPE_ID_COL <- "IID"
GREX_ID_COL <- "FID"        # PrediXcan predicted-expression sample column
GENE_LIST_FILE <- "/path/to/TWAS_results/Korean/R2_filtered_GReX_models.tsv"
SELECTED_PAIRS_FILE <- "/path/to/selected_348_gene_trait_pairs.tsv"
OUTPUT_DIR <- "/path/to/permutation_results"

SEED <- 20260428L
N_NULL_CONTINUOUS <- 1000L
N_NULL_BINARY <- 100L
SELECTED_START <- 1001L
SELECTED_END <- 101000L        # 100,000 ADDITIONAL permutations per pair
EXPECTED_PAIRS <- 348L
EXPECTED_GENES <- 1630L        # check against the study's Korean TWAS gene set
N_CORES <- 1L                 # parallelize across traits on Unix if > 1
MIN_SWAPPABLE_FRACTION <- 0.70
MIN_VALID_FRACTION <- 0.90
AGE_BINS <- 5L
BMI_BINS <- 5L
CHECKPOINT_EVERY <- 100L
BONFERRONI_ALPHA <- 0.05
SAVE_NULL_PVALUES <- TRUE    # retain null P-values for exact lambda_GC / QQ

# The analysis-ready dataset may retain original 1/2/3 phenotype coding.
# If coded as 0/1 already, binary variables are left unchanged.

# 2. Phenotype-specific TWAS models ------------------------------------------
CONTINUOUS <- c(
  "BMI", "HEIGHT", "HIP", "WAIST", "DBP", "SBP", "MAP", "HCT", "HB", "MCH",
  "MCHC", "MCV", "PLAT", "RBC", "WBC", "CREATININE", "eGFR", "BUN",
  "ALT", "AST", "R_GTP", "T_BIL", "ALBUMIN", "GLU", "HbA1C", "HDL",
  "LDL", "TCHL", "TG"
)
BINARY <- c(
  "SMOKE", "BRCA", "UTCA", "COLCA", "GCA", "HCCCA", "LCA", "PROCA",
  "THYCA", "ANG", "ARRHY", "CVA", "CAD", "CHF", "HTN", "MI", "PV", "STR",
  "TIA", "VD", "LIV", "GASTRO", "HEPATITISB", "HEPATITISC", "FLIV", "GB",
  "STOMUL", "CIRRHOSIS", "ULCER", "THY", "CATA", "GLAU", "KD", "LIP",
  "DM", "ARTH", "FRAC", "OSTE", "NOI", "PARK", "ALLER", "GT", "POL",
  "PER", "ASTH", "BRON", "COPD", "CLD", "TB", "UB", "BPH", "UT"
)

config_for_trait <- function(trait) {
  if (trait %in% CONTINUOUS) {
    covariates <- c("SEX", "AGE", "SMOKE", if (!trait %in% c("BMI", "HEIGHT", "HIP", "WAIST")) "BMI")
    group <- if (trait %in% c("BMI", "HEIGHT", "HIP", "WAIST")) "Continuous_All_AnotherCov" else "Continuous_All"
    return(list(type = "continuous", group = group, covariates = covariates, sex = NA_integer_))
  }
  if (!trait %in% BINARY) stop("Unknown trait: ", trait)
  if (trait == "SMOKE") return(list(type = "binary", group = "Binary_SMOKE", covariates = c("SEX", "AGE"), sex = NA_integer_))
  if (trait %in% c("PROCA", "BPH")) return(list(type = "binary", group = "Binary_Male", covariates = c("AGE", "SMOKE", "BMI"), sex = 0L))
  if (trait %in% c("BRCA", "UTCA")) return(list(type = "binary", group = "Binary_Female", covariates = c("AGE", "SMOKE", "BMI"), sex = 1L))
  list(type = "binary", group = "Binary_All", covariates = c("SEX", "AGE", "SMOKE", "BMI"), sex = NA_integer_)
}

# Original selected manifest used B_-prefixed coding. Only map the known
# renamed variables; all phenotype columns in DATA_FILE remain abbreviations.
normalize_manifest_trait <- function(x) {
  x <- sub("^B_", "", trimws(as.character(x)))
  mapping <- c(RBC_B = "RBC", WBC_B = "WBC", CREATINE = "CREATININE", EGFR = "eGFR", GLU0 = "GLU", HBA1C = "HbA1C")
  changed <- x %in% names(mapping)
  x[changed] <- unname(mapping[x[changed]])
  x
}
original_trait_id <- function(x) {
  original <- c(RBC = "B_RBC_B", WBC = "B_WBC_B", CREATININE = "B_CREATINE", eGFR = "B_EGFR", GLU = "B_GLU0", HbA1C = "B_HBA1C")
  if (x %in% names(original)) unname(original[[x]]) else paste0("B_", x)
}

# 3. Utilities and stratified permutation ------------------------------------
need_cols <- function(dt, cols) {
  missing <- setdiff(cols, names(dt))
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))
}
recode_binary <- function(v, name) {
  v <- suppressWarnings(as.numeric(v))
  x <- sort(unique(v[!is.na(v)]))
  if (all(x %in% c(0, 1))) return(v)
  if (name == "SMOKE" && all(x %in% c(1, 2, 3)) && any(x > 1)) return(ifelse(is.na(v), NA, as.numeric(v != 1)))
  if (all(x %in% c(1, 2)) && any(x == 2)) return(ifelse(is.na(v), NA, as.numeric(v == 2)))
  stop("Unsupported coding for binary variable ", name, ": ", paste(x, collapse = ","))
}
quantile_bin <- function(x, n) {
  br <- unique(as.numeric(quantile(x, probs = seq(0, 1, length.out = n + 1L), na.rm = TRUE)))
  if (length(br) < 2L) return(rep("Q1", length(x)))
  as.character(cut(x, breaks = br, include.lowest = TRUE, labels = paste0("Q", seq_len(length(br) - 1L))))
}
group_indices <- function(dt, cols) {
  cols <- cols[cols %in% names(dt)]
  cols <- cols[vapply(cols, function(x) length(unique(na.omit(dt[[x]]))) > 1L, logical(1))]
  if (!length(cols)) return(list(ALL = seq_len(nrow(dt))))
  keys <- do.call(paste, c(lapply(dt[, ..cols], function(x) { x <- as.character(x); x[is.na(x)] <- "NA"; x }), sep = "|"))
  split(seq_len(nrow(dt)), keys)
}
swappable <- function(y, groups) {
  sum(vapply(groups, function(ii) length(ii) >= 2L && length(unique(y[ii])) > 1L, logical(1)) * lengths(groups)) / length(y)
}
binary_strata <- function(dt, trait, group) {
  dt[, AGE_BIN := quantile_bin(AGE, AGE_BINS)]
  dt[, BMI_BIN := quantile_bin(BMI, BMI_BINS)]
  candidates <- switch(group,
    Binary_All = list(c("SEX", "SMOKE", "AGE_BIN", "BMI_BIN"), c("SEX", "SMOKE", "AGE_BIN"), c("SEX", "SMOKE"), "SEX", character()),
    Binary_SMOKE = list(c("SEX", "AGE_BIN"), "SEX", "AGE_BIN", character()),
    Binary_Male = list(c("SMOKE", "AGE_BIN", "BMI_BIN"), c("SMOKE", "AGE_BIN"), "SMOKE", "AGE_BIN", character()),
    Binary_Female = list(c("SMOKE", "AGE_BIN", "BMI_BIN"), c("SMOKE", "AGE_BIN"), "SMOKE", "AGE_BIN", character())
  )
  groups <- lapply(candidates, function(cols) group_indices(dt, cols))
  fractions <- vapply(groups, function(g) swappable(dt[[trait]], g), numeric(1))
  chosen <- which(fractions >= MIN_SWAPPABLE_FRACTION)[1]
  if (is.na(chosen)) chosen <- which.max(fractions)
  list(groups = groups[[chosen]], cols = paste(candidates[[chosen]], collapse = ";"), swappable_fraction = fractions[chosen])
}
seed_hash <- function(key) {
  h <- as.numeric(SEED) %% 2147483647
  for (x in utf8ToInt(enc2utf8(key))) h <- (h * 131 + x) %% 2147483647
  as.integer(max(1, h))
}
permutation_seed <- function(trait_seed, index) {
  as.integer(max(1, (as.numeric(trait_seed) + as.numeric(index) * 1009) %% 2147483647))
}
shuffle_within <- function(y, groups, seed) {
  set.seed(seed)
  out <- y
  for (ii in groups) if (length(ii) > 1L) out[ii] <- y[ii][sample.int(length(ii))]
  out
}

# 4. Fit same tests as primary TWAS ------------------------------------------
fit_statistic <- function(dt, response, gene, covariates, trait_type) {
  cols <- unique(c(response, gene, covariates))
  d <- dt[complete.cases(dt[, ..cols]), ..cols]
  if (nrow(d) < length(covariates) + 4L || length(unique(d[[gene]])) < 2L || length(unique(d[[response]])) < 2L) return(c(stat = NA_real_, p = NA_real_))
  formula <- reformulate(c(gene, covariates), response = response)
  if (trait_type == "continuous") {
    fit <- tryCatch(lm(formula, data = d), error = function(e) NULL)
    if (is.null(fit)) return(c(stat = NA_real_, p = NA_real_))
    tab <- summary(fit)$coefficients
    ii <- match(gene, gsub("`", "", rownames(tab), fixed = TRUE))
    if (is.na(ii)) return(c(stat = NA_real_, p = NA_real_))
    return(c(stat = abs(unname(tab[ii, "t value"])), p = unname(tab[ii, "Pr(>|t|)"])))
  }
  if (!setequal(unique(d[[response]]), c(0, 1))) return(c(stat = NA_real_, p = NA_real_))
  fit <- tryCatch(logistf::logistf(formula, data = d, pl = FALSE), error = function(e) NULL)
  if (is.null(fit)) return(c(stat = NA_real_, p = NA_real_))
  drop <- tryCatch(as.data.frame(drop1(fit, scope = gene, data = d, test = "PLR")), error = function(e) NULL)
  if (is.null(drop) || !nrow(drop)) return(c(stat = NA_real_, p = NA_real_))
  ii <- match(gene, gsub("`", "", rownames(drop), fixed = TRUE))
  if (is.na(ii) && nrow(drop) == 1L) ii <- 1L
  if (is.na(ii)) return(c(stat = NA_real_, p = NA_real_))
  chi_col <- grep("chisq|chi.sq|chi-square", names(drop), ignore.case = TRUE)[1]
  p_col <- grep("p.value|p-value|Pr\\(|prob", names(drop), ignore.case = TRUE)[1]
  if (is.na(chi_col)) chi_col <- 1L
  if (is.na(p_col)) p_col <- 3L
  if (ncol(drop) < max(chi_col, p_col)) return(c(stat = NA_real_, p = NA_real_))
  c(stat = as.numeric(drop[ii, chi_col]), p = as.numeric(drop[ii, p_col]))
}

# 5. Task setup and checkpoint ------------------------------------------------
prepare_task <- function(data, trait, genes, cfg) {
  cols <- unique(c(trait, cfg$covariates, genes, "SEX", "SMOKE", "AGE", "BMI"))
  need_cols(data, cols)
  d <- copy(data[, ..cols])
  if (!is.na(cfg$sex)) d <- d[SEX == cfg$sex]
  required_for_model <- unique(c(trait, cfg$covariates))
  d <- d[complete.cases(d[, ..required_for_model])]
  if (nrow(d) <= length(cfg$covariates) + 3L) stop("Insufficient data for ", trait)
  if (cfg$type == "continuous") {
    null <- lm(reformulate(cfg$covariates, response = trait), data = d)
    if (nobs(null) != nrow(d)) stop("Continuous null model has missing rows for ", trait)
    fitted_y <- unname(fitted(null))
    resid_y <- unname(residuals(null))
    groups <- group_indices(d, c("SEX", "SMOKE"))
    make_y <- function(seed) fitted_y + shuffle_within(resid_y, groups, seed)
    strata_info <- "SEX;SMOKE"
    swappable_frac <- NA_real_
  } else {
    if (!setequal(unique(d[[trait]]), c(0, 1))) stop("Binary outcome is not 0/1 with both classes: ", trait)
    s <- binary_strata(d, trait, cfg$group)
    original_y <- d[[trait]]
    make_y <- function(seed) shuffle_within(original_y, s$groups, seed)
    strata_info <- s$cols
    swappable_frac <- s$swappable_fraction
  }
  list(data = d, make_y = make_y, strata = strata_info, swappable = swappable_frac)
}
write_checkpoint <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".", Sys.getpid(), ".tmp")
  saveRDS(object, tmp)
  if (!file.rename(tmp, path)) { file.copy(tmp, path, overwrite = TRUE); unlink(tmp) }
}
file_stamp <- function(path) {
  info <- file.info(path)
  paste(normalizePath(path), info$size, as.numeric(info$mtime), sep = ":")
}
input_stamps <- function() {
  sources <- if (is.na(DATA_FILE)) c(PHENOTYPE_FILE, KOREAN_GREX_FILE) else DATA_FILE
  vapply(c(sources, GENE_LIST_FILE), file_stamp, character(1))
}
load_analysis_data <- function(genes) {
  if (!is.na(DATA_FILE)) return(fread(DATA_FILE))
  pheno <- fread(PHENOTYPE_FILE)
  need_cols(pheno, PHENOTYPE_ID_COL)
  headers <- names(fread(KOREAN_GREX_FILE, nrows = 0L))
  missing <- setdiff(c(GREX_ID_COL, genes), headers)
  if (length(missing)) stop("Missing predicted-expression columns: ", paste(missing, collapse = ", "))
  grex <- fread(KOREAN_GREX_FILE, select = c(GREX_ID_COL, genes))
  pheno_ids <- as.character(pheno[[PHENOTYPE_ID_COL]])
  grex_ids <- as.character(grex[[GREX_ID_COL]])
  if (anyDuplicated(pheno_ids) || anyDuplicated(grex_ids)) stop("Duplicate sample IDs in phenotype or prediction data")
  idx <- match(pheno_ids, grex_ids)
  keep <- which(!is.na(idx))
  if (!length(keep)) stop("No overlapping samples between phenotype and GReX")
  cbind(pheno[keep], grex[idx[keep], setdiff(names(grex), GREX_ID_COL), with = FALSE])
}
make_fingerprint <- function(object) {
  f <- tempfile()
  on.exit(unlink(f), add = TRUE)
  saveRDS(object, f)
  unname(tools::md5sum(f))
}

run_trait <- function(task, data, mode) {
  trait <- task$trait
  cfg <- config_for_trait(trait)
  genes <- task$genes
  if (!length(genes)) stop("Empty gene set for ", trait)
  prep <- prepare_task(data, trait, genes, cfg)
  dt <- prep$data
  dt[, OBS_PHENO := get(trait)]
  n <- if (mode == "null") if (cfg$type == "continuous") N_NULL_CONTINUOUS else N_NULL_BINARY else SELECTED_END - SELECTED_START + 1L
  start <- if (mode == "null") 1L else SELECTED_START
  end <- start + n - 1L
  # Preserve the original 100k runner's trait-seed labels despite abbreviated
  # input-column names; permutation indices likewise remain 1001:101000.
  sex_label <- if (is.na(cfg$sex)) "ALL" else if (cfg$sex == 0L) "MALE_ONLY" else "FEMALE_ONLY"
  key <- paste(cfg$group, original_trait_id(trait), cfg$type, sex_label, sep = "|")
  tseed <- seed_hash(key)
  target_dir <- file.path(OUTPUT_DIR, mode, cfg$type, trait)
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  checkpoint <- file.path(target_dir, "checkpoint.rds")
  summary_file <- file.path(target_dir, "summary.tsv")
  fingerprint <- make_fingerprint(list(mode, trait, cfg, genes, start, end, SEED, nrow(dt), key,
    prep$strata, MIN_SWAPPABLE_FRACTION, MIN_VALID_FRACTION, AGE_BINS, BMI_BINS,
    input_stamps()))
  if (file.exists(summary_file)) {
    x <- fread(summary_file)
    if ("fingerprint" %in% names(x) && all(x$fingerprint == fingerprint)) return(as.data.table(x))
    stop("Existing output has different settings: ", summary_file)
  }
  if (mode == "selected348") {
    observed <- vapply(genes, function(g) fit_statistic(dt, "OBS_PHENO", g, cfg$covariates, cfg$type)["stat"], numeric(1))
    if (any(!is.finite(observed))) stop("Observed fit failed: ", trait, " / ", paste(genes[!is.finite(observed)], collapse = ","))
  }
  k <- valid <- failed <- integer(length(genes))
  perm_summaries <- if (mode == "null") vector("list", n) else NULL
  null_p <- if (mode == "null") matrix(NA_real_, nrow = n, ncol = length(genes),
    dimnames = list(as.character(seq.int(start, end)), genes)) else NULL
  last <- start - 1L
  if (file.exists(checkpoint)) {
    ck <- readRDS(checkpoint)
    if (!identical(ck$fingerprint, fingerprint)) stop("Incompatible checkpoint: ", checkpoint)
    k <- ck$k; valid <- ck$valid; failed <- ck$failed
    perm_summaries <- ck$perm_summaries; null_p <- ck$null_p; last <- ck$last
  }
  message("[", mode, "] ", trait, ": ", length(genes), " genes; permutations ", start, "-", end, "; resume from ", last + 1L)
  if (last < end) for (b in seq.int(last + 1L, end)) {
    seed <- permutation_seed(tseed, b)
    dt[, PERM_PHENO := prep$make_y(seed)]
    fits <- lapply(genes, function(g) fit_statistic(dt, "PERM_PHENO", g, cfg$covariates, cfg$type))
    stats <- vapply(fits, function(z) unname(z["stat"]), numeric(1))
    pvals <- vapply(fits, function(z) unname(z["p"]), numeric(1))
    ok <- is.finite(stats) & is.finite(pvals) & pvals >= 0 & pvals <= 1
    valid <- valid + as.integer(ok)
    failed <- failed + as.integer(!ok)
    if (mode == "selected348") k <- k + as.integer(ok & stats >= observed)
    if (mode == "null") {
      pv <- pvals[ok]
      null_p[b - start + 1L, ] <- ifelse(ok, pvals, NA_real_)
      chisq <- qchisq(pmax(pv, .Machine$double.xmin), df = 1, lower.tail = FALSE)
      lambda <- if (length(pv)) median(chisq) / qchisq(0.5, 1, lower.tail = FALSE) else NA_real_
      perm_summaries[[b - start + 1L]] <- data.table(index = b, n_valid = length(pv), n_failed = sum(!ok),
        lambda_GC = lambda, prop_p_below_005 = if (length(pv)) mean(pv < 0.05) else NA_real_)
    }
    if ((b - start + 1L) %% CHECKPOINT_EVERY == 0L || b == end) {
      write_checkpoint(list(fingerprint = fingerprint, k = k, valid = valid, failed = failed,
        null_p = null_p, perm_summaries = perm_summaries, last = b), checkpoint)
      message("  ", trait, ": ", b - start + 1L, "/", n, " permutations")
    }
  }
  result <- data.table(trait = trait, trait_type = cfg$type, analysis_group = cfg$group,
    n_samples = nrow(dt), gene = genes, permutations_requested = n,
    permutations_valid = valid, permutations_failed = failed,
    chosen_strata = prep$strata, swappable_fraction = prep$swappable, fingerprint = fingerprint)
  if (mode == "selected348") {
    result[, `:=`(observed_statistic = observed, exceedances_K = k,
      monte_carlo_p = fifelse(valid / n >= MIN_VALID_FRACTION, (k + 1) / (valid + 1), NA_real_),
      min_attainable_p = 1 / (valid + 1), bonferroni_threshold = BONFERRONI_ALPHA / EXPECTED_GENES)]
    result[, bonferroni_pass := !is.na(monte_carlo_p) & monte_carlo_p <= bonferroni_threshold]
  } else {
    perm_dt <- rbindlist(perm_summaries, fill = TRUE)
    fwrite(perm_dt, file.path(target_dir, "null_permutation_summary.tsv"), sep = "\t")
    all_p <- as.vector(null_p)
    all_p <- all_p[is.finite(all_p) & all_p >= 0 & all_p <= 1]
    lambda_pool <- if (length(all_p)) median(qchisq(pmax(all_p, .Machine$double.xmin), df = 1, lower.tail = FALSE)) /
      qchisq(0.5, 1, lower.tail = FALSE) else NA_real_
    if (SAVE_NULL_PVALUES) saveRDS(null_p, file.path(target_dir, "null_pvalues.rds"), compress = TRUE)
    result[, `:=`(lambda_GC_pooled = lambda_pool,
      lambda_GC_permutation_median = median(perm_dt$lambda_GC, na.rm = TRUE),
      null_p_values = length(all_p), null_p_below_005 = sum(all_p < 0.05))]
  }
  fwrite(result, summary_file, sep = "\t")
  unlink(checkpoint)
  result
}

# 6. Input, tasks, and execution ---------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) tolower(args[1]) else "null"
trait_filter <- if (length(args) >= 2L) args[2] else NULL
if (!mode %in% c("null", "selected348", "all")) stop("Usage: Rscript 08_run_permutation.R [null|selected348|all] [optional_trait]")
input_files <- if (is.na(DATA_FILE)) c(PHENOTYPE_FILE, KOREAN_GREX_FILE) else DATA_FILE
for (f in c(input_files, GENE_LIST_FILE)) if (!file.exists(f)) stop("Missing input: ", f)
model_genes <- fread(GENE_LIST_FILE)$gene
if (length(unique(model_genes)) != EXPECTED_GENES) stop("Expected ", EXPECTED_GENES, " R2-filtered genes; found ", length(unique(model_genes)))
if (anyDuplicated(model_genes)) stop("Duplicated genes in GENE_LIST_FILE")
all_data <- load_analysis_data(model_genes)
need_cols(all_data, c("SEX", "AGE", "SMOKE", "BMI", model_genes))
if (anyDuplicated(names(all_data))) stop("Duplicate column names in analysis-ready data")
all_data[, SEX := recode_binary(SEX, "SEX")]
all_data[, SMOKE := recode_binary(SMOKE, "SMOKE")]
for (trait in setdiff(intersect(BINARY, names(all_data)), "SMOKE")) all_data[, (trait) := recode_binary(get(trait), trait)]

selected <- NULL
if (mode %in% c("selected348", "all")) {
  if (!file.exists(SELECTED_PAIRS_FILE)) stop("Missing 348-pair manifest: ", SELECTED_PAIRS_FILE)
  selected <- fread(SELECTED_PAIRS_FILE)
  need_cols(selected, c("trait_id", "Gene"))
  if (nrow(selected) != EXPECTED_PAIRS) stop("Expected exactly ", EXPECTED_PAIRS, " rows in selected manifest; found ", nrow(selected))
  selected[, trait := normalize_manifest_trait(trait_id)]
  if (anyDuplicated(selected[, .(trait, Gene)])) stop("Duplicate selected gene-trait pairs")
  if (any(!selected$Gene %in% model_genes)) stop("Selected genes absent from the Korean R2-filtered model list")
  if (any(!selected$trait %in% c(CONTINUOUS, BINARY))) stop("Selected manifest has unrecognized traits: ", paste(setdiff(selected$trait, c(CONTINUOUS, BINARY)), collapse = ", "))
  if ("analysis_group" %in% names(selected)) {
    expected_groups <- vapply(selected$trait, function(x) config_for_trait(x)$group, character(1))
    if (any(selected$analysis_group != expected_groups)) stop("Selected manifest analysis_group inconsistent with primary TWAS settings")
  }
  need_cols(all_data, unique(selected$trait))
}
run_mode <- function(which_mode) {
  if (which_mode == "null") {
    traits <- c(CONTINUOUS, BINARY)
    need_cols(all_data, traits)
    jobs <- lapply(traits, function(t) list(trait = t, genes = model_genes))
  } else {
    traits <- unique(selected$trait)
    jobs <- lapply(traits, function(t) list(trait = t, genes = selected[trait == t, Gene]))
  }
  if (!is.null(trait_filter)) jobs <- Filter(function(x) x$trait == trait_filter, jobs)
  if (!length(jobs)) stop("No tasks for requested trait: ", trait_filter)
  worker <- function(job) tryCatch(run_trait(job, all_data, which_mode), error = function(e) {
    message("ERROR [", which_mode, "] ", job$trait, ": ", conditionMessage(e))
    data.table(trait = job$trait, gene = job$genes, status = "ERROR", error = conditionMessage(e))
  })
  if (N_CORES > 1L && .Platform$OS.type != "windows") {
    rows <- parallel::mclapply(jobs, worker, mc.cores = N_CORES, mc.preschedule = FALSE)
  } else rows <- lapply(jobs, worker)
  summary <- rbindlist(rows, fill = TRUE)
  if (which_mode == "selected348") {
    manifest <- copy(selected)
    manifest[, manifest_order := .I]
    summary <- merge(manifest, summary, by.x = c("trait", "Gene"), by.y = c("trait", "gene"), all.x = FALSE, sort = FALSE)
    setorder(summary, manifest_order)
    if (is.null(trait_filter) && nrow(summary) != EXPECTED_PAIRS) stop("Selected-pair summary must contain exactly 348 rows")
  }
  dir.create(file.path(OUTPUT_DIR, which_mode), recursive = TRUE, showWarnings = FALSE)
  fwrite(summary, file.path(OUTPUT_DIR, which_mode, "ALL_trait_results.tsv"), sep = "\t")
  if ("status" %in% names(summary) && any(summary$status == "ERROR", na.rm = TRUE)) stop("One or more permutation tasks failed. See ALL_trait_results.tsv")
  message("Completed mode: ", which_mode, "; traits = ", length(jobs), "; rows = ", nrow(summary))
}
if (mode %in% c("null", "all")) run_mode("null")
if (mode %in% c("selected348", "all")) run_mode("selected348")
