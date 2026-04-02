# ProteinBatcher

**ProteinBatcher** is a Bioconductor-style R package for **standardized, design downstream analysis of label-free quantitative proteomics** tables (with a focus on **DIA-NN-style DIA outputs**).
It provides utilities to **read, validate and harmonize quantitative proteomics tables**, convert them to **SummarizedExperiment**, run **imputation-aware differential abundance testing with limma**,
and optionally generate **publication-ready volcano plots and dereg ul ograms**.

---

## Why ProteinBatcher?

Bioconductor already provides excellent tools for statistical modeling and container infrastructure (notably **limma** and **SummarizedExperiment**), but **end-to-end proteomics workflows** that combine:

- missingness-based filtering,
- mean/LDV imputation tailored to label-free proteomics,
- flexible multi-variable designs (batch/sex/other covariates),
- correlation-aware models via `duplicateCorrelation`,
- and standardized output objects for downstream reporting

are less common in a single, reproducible pipeline.

ProteinBatcher focuses on **orchestrating these steps** in a consistent way, so repeated analyses across experiments become less error-prone and more comparable.
---

## Core features

### 1) Standardized pipeline wrapper
The main entry point is:

- `run_proteomics_pipeline()`

It performs: import -> missingness filtering -> imputation -> limma testing -> output organization -> (optional) plots and tables.

### 2) Imputation-aware differential testing
`test_limma_customized()` fits limma models and stores results in `rowData(se)`.

Key feature: **BH-adjusted p-values are recomputed after excluding features imputed as LDV in both groups of a contrast (LDV/LDV)**. This avoids artificially significant results driven by left-censored imputations.

### 3) Flexible experimental designs (including interactions)
You provide a `formula` (e.g. `~ 0 + condition + batch + condition:batch`) plus contrasts (`tests`) and optional interaction tests (`tests_interaction`).

### 4) Paired and correlation-aware designs
- `paired = TRUE`: ensures `replicate` is included in the design for paired/technical repeated measures.
- `block_effect = TRUE`: uses `limma::duplicateCorrelation()` with `colData(se)$block` and fits correlation-aware models.

### 5) Plotting for interpretation
- `plot_volcano_customized()` produces volcano plots from `SummarizedExperiment` results stored in `rowData()`.
- `plot_deregulogram()` produces deregulograms for **2-level interaction factors** (e.g. female/male, batch 1/batch 2), comparing **full effects** between levels.

**Important interpretation note:** the deregulogram is not expected to label the exact same proteins as the interaction volcano plot. The interaction volcano highlights significant interaction coefficients;
the deregulogram emphasizes a more interpretable subset where interaction occurs in the context of statistically significant condition/full effects. 

---

## Installation

### Bioconductor (recommended once available)
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ProteinBatcher")
``
