# KoreanTWAS

Code repository for the study:

**Transcriptome-Wide Association Studies of 81 Traits in 79,294 Korean Individuals**

This repository contains the study-specific analysis scripts used to train and apply a Korean whole-blood genetically regulated gene expression (GReX) model, perform transcriptome-wide association studies (TWAS), conduct conditional fine-mapping, evaluate external replication, prepare gene set enrichment analyses, perform computational drug-repurposing analyses, and assess TWAS null calibration and association-level robustness using phenotype permutations.

## Overview

The analysis workflow consists of eight main steps:

1. **Train a Korean whole-blood GReX model** using the PredictDB-Tutorial framework.
2. **Apply Korean- and GTEx-based GReX models** and evaluate prediction performance.
3. **Perform individual-level TWAS** for 81 traits in KoGES+GENIE.
4. **Perform GIFT conditional fine-mapping** of Bonferroni-significant Korean-based TWAS associations.
5. **Perform external S-PrediXcan analyses** using BioBank Japan (BBJ) and China Kadoorie Biobank (CKB) GWAS summary statistics.
6. **Prepare ranked gene lists for gene set enrichment analysis (GSEA)** using WebGestalt.
7. **Prepare and summarize Connectivity Map (CMap) drug-repurposing analyses** using the CLUE platform.
8. **Perform phenotype-permutation analyses** to assess TWAS null calibration and the robustness of 348 Bonferroni-significant associations.

The repository contains study-specific analysis code only. Individual-level cohort data and third-party software are not redistributed.

---

## Repository structure

```text
KoreanTWAS/
├── 01_train_GReX.R
├── 02_predict_GReX.R
├── 03_run_TWAS.R
├── 04_run_GIFT.R
├── 05_external_replication.R
├── 06_prepare_GSEA.R
├── 07_DrugRepurposing.R
├── 08_run_permutation.R
└── README.md
```

All file paths in the scripts are placeholders and must be replaced with local paths before execution.

---

## Analysis workflow

### 1. Korean GReX model training

Script:

```text
01_train_GReX.R
```

The Korean whole-blood GReX model was trained using the PredictDB-Tutorial implementation:

https://github.com/hakyimlab/PredictDB-Tutorial

Required preprocessed inputs include:

```text
<input_dir>/
├── gene_annot.parsed.txt
├── genotype/
│   ├── snp_annot.chr1.txt
│   ├── ...
│   ├── snp_annot.chr22.txt
│   ├── genotype.chr1.txt
│   ├── ...
│   └── genotype.chr22.txt
├── expression/
│   └── transformed_expression.txt
└── PEER/
    └── covariates.txt
```

Model-training settings used in the study:

- Minor allele frequency threshold: `MAF >= 0.01`
- cis-window: `±1 Mb`
- Elastic-Net mixing parameter: `alpha = 0.5`
- Inner cross-validation: `10 folds`
- Outer nested cross-validation: `5 folds`

Example:

```bash
Rscript 01_train_GReX.R all \
  /path/to/input \
  /path/to/model_output \
  /path/to/PredictDB-Tutorial
```

PredictDB database construction and filtering were subsequently performed following the **Make a database** and **Filter the database** procedures in PredictDB-Tutorial. The study retained models satisfying:

```text
zscore_pval < 0.05
rho_avg > 0.10
```

The sample population label was set to `KOR`.

---

### 2. GReX prediction and prediction-performance evaluation

Script:

```text
02_predict_GReX.R
```

Predicted expression was generated using the original PrediXcan implementation:

https://github.com/hakyimlab/PrediXcan

The GTEx whole-blood Elastic-Net model was obtained from PredictDB:

https://predictdb.org/post/2021/07/21/gtex-v8-models-on-eqtl-and-sqtl/elastic_net_eqtl.tar

The following analyses are implemented:

| Analysis | GReX model | Dataset |
|---|---|---|
| Internal prediction performance | Korean | Asan+Chosun |
| Internal prediction performance | GTEx | GTEx |
| External prediction performance | GTEx | Asan+Chosun |
| External prediction performance | Korean | CODA |
| External prediction performance | GTEx | CODA |
| GReX prediction for TWAS | Korean | KoGES+GENIE |
| GReX prediction for TWAS | GTEx | KoGES+GENIE |

For external prediction-performance analyses, predicted and observed expression are compared using:

- gene-wise Pearson correlation;
- gene-wise \(R^2\); and
- sample-wise Pearson correlation.

PrediXcan-ready dosage files are assumed to have been generated beforehand using `convert_plink_to_dosage.py`.

Example:

```bash
Rscript 02_predict_GReX.R all
```

Available modes are:

```text
all
cv
external
twas
```

---

### 3. Individual-level TWAS

Script:

```text
03_run_TWAS.R
```

TWAS was performed in **79,294 KoGES+GENIE participants** across **81 traits** using both Korean- and GTEx-based whole-blood GReX models.

The analysis included:

- 29 continuous traits;
- 52 binary traits.

Continuous traits were analyzed using linear regression, and binary traits were analyzed using Firth bias-reduced logistic regression.

GReX models were restricted to genes satisfying:

```text
rho_avg_squared > 0.05
```

Genes for which more than 90% of predicted-expression values were zero were excluded from association testing.

The number of GReX models satisfying `rho_avg_squared > 0.05` was used as the Bonferroni multiple-testing denominator for each prediction model.

#### Covariate adjustment

The primary TWAS models were adjusted for:

- age;
- sex;
- body mass index (BMI); and
- smoking status.

No genotype principal components were included as covariates in the individual-level TWAS.

Covariate exceptions reflected the phenotype being analyzed:

- BMI was omitted as a covariate for analyses of **BMI, height, hip circumference, and waist circumference**.
- For the **smoking-status** analysis, BMI and smoking status were not included as covariates.
- Sex was not included as a covariate in sex-specific analyses.

Sex-specific analyses were performed for:

- male-only: prostate cancer (`PROCA`) and benign prostatic hyperplasia (`BPH`);
- female-only: breast cancer (`BRCA`) and uterine cancer (`UTCA`).

Example:

```bash
Rscript 03_run_TWAS.R all
```

Available modes are:

```text
all
Korean
GTEx
```

---

### 4. GIFT conditional fine-mapping

Script:

```text
04_run_GIFT.R
```

Bonferroni-significant Korean-based TWAS associations were further evaluated using GIFT conditional fine-mapping.

Candidate genes were defined as GReX models with:

```text
rho_avg_squared > 0.05
```

that overlapped the East Asian linkage-disequilibrium block containing the focal TWAS gene.

GRCh38 East Asian LD blocks were obtained from:

https://github.com/jmacdon/LDblocks_GRCh38

For each candidate gene, cis-genotypes were extracted within:

```text
±100 kb
```

A minimum of 10 common cis-SNPs between the Asan+Chosun and KoGES+GENIE datasets was required.

GIFT settings used in the study were:

```text
maxiter = 100
tol     = 1e-3
pleio   = 0
filter  = TRUE
```

Asan+Chosun expression used for GIFT was covariate-adjusted before running this script. Expression adjustment included:

- sex;
- age;
- three genotype principal components; and
- 60 PEER factors.

The GIFT workflow can be run in separate stages:

```bash
Rscript 04_run_GIFT.R prepare
Rscript 04_run_GIFT.R run
Rscript 04_run_GIFT.R summarize
```

or together:

```bash
Rscript 04_run_GIFT.R all
```

A single trait can optionally be specified:

```bash
Rscript 04_run_GIFT.R all MCHC
```

The summary output contains both:

- regional Bonferroni-adjusted GIFT P values; and
- transcriptome-wide Bonferroni-adjusted GIFT P values.

---

### 5. External replication using BBJ and CKB

Script:

```text
05_external_replication.R
```

External association analyses were performed using S-PrediXcan from MetaXcan:

https://github.com/hakyimlab/MetaXcan

Specifically:

```text
MetaXcan/software/SPrediXcan.py
```

The Korean whole-blood GReX model and its corresponding covariance matrix were applied to publicly available GWAS summary statistics from:

- **BioBank Japan (BBJ)**  
  https://pheweb.jp/downloads

- **China Kadoorie Biobank (CKB)**  
  https://pheweb.ckbiobank.org/

The script contains separate GWAS column mappings for:

- BBJ continuous traits;
- BBJ binary traits; and
- CKB traits.

Example:

```bash
Rscript 05_external_replication.R all
```

Individual modes are:

```text
BBJ_continuous
BBJ_binary
CKB
collect
```

After the cohort-specific analyses have completed, results can be combined using:

```bash
Rscript 05_external_replication.R collect
```

This generates:

```text
S_PrediXcan_external_replication_ALL.tsv
```

Phenotype harmonization and interpretation of replication results should be performed using the phenotype definitions reported in the study and the corresponding BBJ/CKB phenotype documentation.

---

### 6. Gene set enrichment analysis

Script:

```text
06_prepare_GSEA.R
```

The script prepares trait-specific ranked gene lists for GSEA using the signed Korean-based TWAS Z-score.

GSEA was performed separately for each trait using WebGestalt:

https://www.webgestalt.org/

The following gene-set collections were evaluated:

1. MSigDB Hallmark
2. KEGG
3. Reactome
4. Gene Ontology Biological Process, noRedundant

Study settings were:

- analysis method: GSEA;
- ranked input: gene symbol and signed TWAS Z-score;
- minimum gene-set size: 10;
- maximum gene-set size: 500;
- permutations: 1,000;
- significance threshold: `FDR < 0.05`.

The R script prepares `.rnk` files only. The enrichment analysis itself was performed through the WebGestalt web interface.

Example:

```bash
Rscript 06_prepare_GSEA.R
```

---

### 7. Computational drug repurposing

Script:

```text
07_DrugRepurposing.R
```

Computational drug repurposing was performed using the Connectivity Map approach through the CLUE platform:

https://clue.io/

For each trait, the query signature consisted of:

- the 50 genes with the strongest positive TWAS associations; and
- the 50 genes with the strongest negative TWAS associations.

Genes were ranked using the signed Korean-based TWAS Z-score.

Prepare CLUE input files using:

```bash
Rscript 07_DrugRepurposing.R prepare
```

The generated signatures are submitted manually to CLUE. After downloading the CLUE `ps_pert_summary.gctx` results, post-processing can be performed using:

```bash
Rscript 07_DrugRepurposing.R summarize
```

Only small-molecule perturbagens (`pert_type == "trt_cp"`) are retained.

Connectivity is summarized using the CMap connectivity score \(\tau\). Negative values indicate inverse connectivity between the TWAS-derived trait signature and the perturbagen-induced expression signature.

Compounds satisfying:

```text
tau <= -90
```

were considered candidate repurposing drugs.

For clinical prioritization, compounds were restricted to FDA-approved, non-withdrawn drugs based on:

- **DrugBank**  
  https://go.drugbank.com/

Mechanisms of action were annotated using the CLUE Drug Repurposing Hub.

---

### 8. Phenotype-permutation analyses

Script:

```text
08_run_permutation.R
```

Phenotype permutations were used for two distinct purposes: **genome-wide null calibration** of the Korean-based TWAS and **association-level robustness assessment** of the 348 Bonferroni-significant Korean TWAS gene-trait associations. The analyses use the same phenotype-specific covariates and regression families as the primary individual-level TWAS.

| Mode | Analysis target | Phenotype permutations | Output purpose |
|---|---|---:|---|
| `null` (continuous) | All 1,630 Korean TWAS genes across 29 continuous traits | 1,000 per trait | Null P-value distributions and genomic inflation factors |
| `null` (binary) | All 1,630 Korean TWAS genes across 52 binary traits | 100 per trait | Null P-value distributions and genomic inflation factors |
| `selected348` | The original 348 Bonferroni-significant Korean TWAS gene-trait pairs | 100,000 additional per pair | Association-level Monte Carlo permutation P values |

#### Permutation schemes

- **Continuous traits:** Freedman-Lane residual permutations under the covariate-only linear model, with residuals shuffled within sex/smoking strata. The test statistic is the absolute gene-coefficient *t* statistic.
- **Binary traits:** Phenotype labels are shuffled within strata defined adaptively from sex, smoking, age-quantile and BMI-quantile groupings, as applicable to each trait. The test statistic is the profile-likelihood-ratio chi-square from Firth logistic regression (`logistf`).
- **Trait-specific covariates:** The same exceptions as in the primary TWAS are applied for BMI, height, hip, waist, smoking status, and sex-specific diseases. Genotype principal components are not added as TWAS covariates.

For `null`, the script evaluates the regression P values obtained under each permuted phenotype, summarizes their distributions, and calculates `lambda_GC`. Trait-specific null P-value matrices are saved for further calibration checks and QQ plots; this mode does **not** estimate extreme-tail empirical P values for observed associations.

For `selected348`, the script uses **the original 348-pair manifest**, rather than reselecting genes from a new TWAS run. It performs 100,000 **additional** phenotype permutations per pair, indexed `1001:101000`, preserving the original script's seed scheme. For each association, the empirical Monte Carlo P value is:

```text
P_perm = (K + 1) / (B_valid + 1)
```

Here, `K` is the number of valid permutation test statistics at least as extreme as the observed statistic, and `B_valid` is the number of valid permutations. Runs with a valid-permutation fraction below 0.90 are not assigned an empirical P value. No generalized Pareto distribution (GPD) tail extrapolation is used. The 100,000-permutation results are reported separately; the script does not automatically combine them with an earlier 1,000-permutation analysis.

#### Required inputs

Before running the script, configure the file paths near its beginning:

- `PHENOTYPE_FILE`: KoGES+GENIE analysis-ready phenotype data with abbreviation-based trait columns (such as `BMI`, `HDL`, `SMOKE`, `LIP`), `IID`, `SEX`, and `AGE`.
- `KOREAN_GREX_FILE`: Korean-based predicted-expression data for the same individuals, with ENSG gene identifiers as columns.
- `GENE_LIST_FILE`: `TWAS_results/Korean/R2_filtered_GReX_models.tsv` generated by `03_run_TWAS.R`, containing the 1,630 tested genes (`gene` column).
- `SELECTED_PAIRS_FILE`: original `selected_348_gene_trait_pairs.tsv`, containing **exactly 348 rows**, with columns `trait_id` and `Gene` (and optionally `analysis_group`). This manifest must be supplied separately; the script does not infer or recreate the selected 348 associations.
- `OUTPUT_DIR`: writable directory for permutation outputs.

Alternatively, `DATA_FILE` may point to a combined analysis-ready phenotype/predicted-expression table. Adjust `GREX_ID_COL` if the sample identifier in the predicted-expression file is not `FID`. Other relevant settings include `SEED`, `N_CORES`, `CHECKPOINT_EVERY`, and the permutation counts.

#### Execution

Run the full null-calibration analysis:

```bash
Rscript 08_run_permutation.R null
```

Run the 100,000 additional permutations for the selected 348 associations:

```bash
Rscript 08_run_permutation.R selected348
```

Or run both modes:

```bash
Rscript 08_run_permutation.R all
```

For a limited trial run, supply a trait abbreviation:

```bash
Rscript 08_run_permutation.R null HDL
Rscript 08_run_permutation.R selected348 HDL
```

The code runs permutations on the fly and uses trait-level checkpoints to resume interrupted analyses. Individual-level permutation phenotypes and the complete sequence of permutation test statistics are not exported.

#### Outputs

```text
<OUTPUT_DIR>/
|-- null/
|   |-- ALL_trait_results.tsv
|   |-- continuous/<TRAIT>/
|   |   |-- summary.tsv
|   |   |-- null_permutation_summary.tsv
|   |   `-- null_pvalues.rds
|   `-- binary/<TRAIT>/
|       |-- summary.tsv
|       |-- null_permutation_summary.tsv
|       `-- null_pvalues.rds
`-- selected348/
    |-- ALL_trait_results.tsv
    |-- continuous/<TRAIT>/summary.tsv
    `-- binary/<TRAIT>/summary.tsv
```

The `selected348/ALL_trait_results.tsv` table preserves the input manifest columns and reports the observed test statistic, exceedance count (`K`), number of valid and failed permutations, empirical P value, and Bonferroni threshold (`0.05 / 1,630`). Full-scale runs, especially the binary Firth regressions, can be computationally intensive; available memory and compute resources should be checked before execution.

---

## Software and R packages

The scripts use the following R packages (including `logistf` and `data.table` for permutation analyses):

```text
data.table
dplyr
DBI
RSQLite
logistf
stringr
genio
GIFT
cmapR
```

External software required for different parts of the workflow includes:

- R
- Python 3
- PLINK
- PredictDB-Tutorial
- PrediXcan
- MetaXcan / S-PrediXcan

Because versions can affect reproducibility, users should record the software versions or Git commit hashes used in their local environment.

---

## Input data and access

Individual-level genotype, phenotype, and transcriptomic data are not redistributed through this repository.

Access to restricted cohort data is subject to the corresponding data-access procedures.

Public resources used in the study include:

| Resource | Purpose | Link |
|---|---|---|
| PredictDB | GTEx v8 Elastic-Net GReX model | https://predictdb.org/ |
| PredictDB-Tutorial | Korean GReX model training | https://github.com/hakyimlab/PredictDB-Tutorial |
| PrediXcan | Individual-level predicted expression | https://github.com/hakyimlab/PrediXcan |
| MetaXcan | S-PrediXcan external analyses | https://github.com/hakyimlab/MetaXcan |
| BBJ PheWeb | External GWAS summary statistics | https://pheweb.jp/downloads |
| CKB PheWeb | External GWAS summary statistics | https://pheweb.ckbiobank.org/ |
| EAS LD blocks | GIFT regional definition | https://github.com/jmacdon/LDblocks_GRCh38 |
| WebGestalt | Gene set enrichment analysis | https://www.webgestalt.org/ |
| CLUE | Connectivity Map analysis | https://clue.io/ |
| DrugBank | FDA-approved drug annotation | https://go.drugbank.com/ |

---

## Reproducibility notes

Several steps depend on resources that cannot be redistributed directly in this repository:

1. Individual-level cohort genotype, phenotype, and expression data require authorized access.
2. PrediXcan, MetaXcan, PredictDB-Tutorial, and PLINK must be installed separately.
3. WebGestalt GSEA was performed through the WebGestalt web interface.
4. Connectivity Map analysis was performed through the CLUE web platform.
5. DrugBank-derived annotations are subject to DrugBank's applicable access and licensing terms.
6. The permutation analyses require the original 348-pair manifest, the 1,630-gene model list, and authorized access to individual-level phenotype and predicted-expression data; the 100,000-permutation results are not automatically pooled with previous runs.

Accordingly, this repository should be interpreted as the collection of **study-specific analysis scripts and analysis settings** rather than a redistribution of all underlying datasets and third-party software.

---

## Citation

If you use this repository, please cite the associated manuscript.

The final manuscript citation can be added here after publication.
