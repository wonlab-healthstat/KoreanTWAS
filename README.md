# KoreanTWAS

Code repository for the study:

**Transcriptome-Wide Association Studies of 81 Traits in 79,294 Korean Individuals**

This repository contains the study-specific analysis scripts used to train and apply a Korean whole-blood genetically regulated gene expression (GReX) model, perform transcriptome-wide association studies (TWAS), conduct conditional fine-mapping, evaluate external replication, prepare gene set enrichment analyses, and perform computational drug-repurposing analyses.

## Overview

The analysis workflow consists of seven main steps:

1. **Train a Korean whole-blood GReX model** using the PredictDB-Tutorial framework.
2. **Apply Korean- and GTEx-based GReX models** and evaluate prediction performance.
3. **Perform individual-level TWAS** for 81 traits in KoGES+GENIE.
4. **Perform GIFT conditional fine-mapping** of Bonferroni-significant Korean-based TWAS associations.
5. **Perform external S-PrediXcan analyses** using BioBank Japan (BBJ) and China Kadoorie Biobank (CKB) GWAS summary statistics.
6. **Prepare ranked gene lists for gene set enrichment analysis (GSEA)** using WebGestalt.
7. **Prepare and summarize Connectivity Map (CMap) drug-repurposing analyses** using the CLUE platform.

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

- **DrugBank version 5.1.12**  
  https://go.drugbank.com/

Mechanisms of action were annotated using the CLUE Drug Repurposing Hub.

---

## Software and R packages

The scripts use the following R packages:

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

Accordingly, this repository should be interpreted as the collection of **study-specific analysis scripts and analysis settings** rather than a redistribution of all underlying datasets and third-party software.

---

## Suggested execution order

A typical analysis sequence is:

```bash
# 1. Train Korean GReX models
Rscript 01_train_GReX.R all <input_dir> <output_dir> <PredictDB-Tutorial_dir>

# 2. Generate predicted expression and evaluate prediction performance
Rscript 02_predict_GReX.R all

# 3. Run Korean- and GTEx-based individual-level TWAS
Rscript 03_run_TWAS.R all

# 4. Run GIFT fine-mapping
Rscript 04_run_GIFT.R all

# 5. Run external S-PrediXcan analyses
Rscript 05_external_replication.R all
Rscript 05_external_replication.R collect

# 6. Prepare WebGestalt GSEA ranked files
Rscript 06_prepare_GSEA.R

# 7. Prepare CLUE signatures
Rscript 07_DrugRepurposing.R prepare

# After manual CLUE analysis and result download:
Rscript 07_DrugRepurposing.R summarize
```

---

## Code availability statement

A manuscript-level code availability statement can be written as:

> Study-specific scripts for GReX model training and prediction, individual-level TWAS, GIFT conditional fine-mapping, external S-PrediXcan analyses, and the preparation and post-processing of WebGestalt and Connectivity Map analyses are publicly available in this GitHub repository. Third-party software is available from the respective developers, and individual-level data are subject to the access policies of the corresponding cohorts.

---

## Citation

If you use this repository, please cite the associated manuscript.

The final manuscript citation can be added here after publication.

---

## License

Please specify the repository license before public release (for example, MIT, BSD-3-Clause, or another license compatible with the study and institutional requirements).

Third-party software and external datasets remain subject to their respective licenses and terms of use.
