test_that("run_proteomics_pipeline runs end-to-end without plots", {

    # ---- Locate extdata shipped with the package ----
    path_annotation <- system.file(
        "extdata", "annotation_HaCaT.tsv",
        package = "ProteinBatcher"
    )

    path_pgmatrix <- system.file(
        "extdata",
        "2024MK017_HaCaT_Stimulation_1to24_Astral_report.pg_matrix.tsv",
        package = "ProteinBatcher"
    )

    # Defensive checks (skip on CRAN/BioC if files missing)
    if (path_annotation == "" || path_pgmatrix == "") {
        testthat::skip("Required extdata files not found.")
    }

    # Temporary output directory (nothing written when plots = FALSE)
    path_output <- tempdir()

    # ---- Workflow parameters (as in vignette example) ----
    experiment <- "HaCaT"
    percent_missing <- 50

    formula <- ~ 0 + condition + batch + condition:batch

    tests <- c(
        "IL13_vs_NoTreated",
        "IL22_vs_NoTreated",
        "COMBO_vs_NoTreated"
    )

    tests_interaction <- c(
        "IL13.batchday2",
        "IL22.batchday2",
        "COMBO.batchday2"
    )
    reference_condition <- "NoTreated"
    # ---- Run pipeline (plots = FALSE) ----
    res <- ProteinBatcher::run_proteomics_pipeline(
        path_pgmatrix       = path_pgmatrix,
        path_annotation     = path_annotation,
        path_output         = path_output,
        tests               = tests,
        tests_interaction   = tests_interaction,
        formula             = formula,
        reference_condition = reference_condition,
        percent_missing     = percent_missing,
        ldv_source          = "per-condition",
        experiment          = experiment,
        plots               = FALSE
    )

    # ---- Basic structure checks ----
    expect_true(is.list(res))
    expect_true(all(
        c("se_raw", "se_filt", "removed", "se_imp") %in% names(res)
    ))

    # ---- SummarizedExperiment outputs ----
    expect_s4_class(res$se_raw,  "SummarizedExperiment")
    expect_s4_class(res$se_filt, "SummarizedExperiment")
    expect_s4_class(res$se_imp,  "SummarizedExperiment")

    # ---- Filtering step sanity ----
    expect_true(nrow(res$se_filt) <= nrow(res$se_raw))
    expect_true(is.data.frame(res$removed))

    # ---- Imputation sanity ----
    X_imp <- SummarizedExperiment::assay(res$se_imp)
    expect_true(is.numeric(X_imp))
    expect_false(anyNA(X_imp))

    # Metadata added during imputation
    md <- S4Vectors::metadata(res$se_imp)
    expect_true("imputation_map" %in% names(md))

    rd <- SummarizedExperiment::rowData(res$se_imp)
    expect_true("imputed" %in% colnames(rd))
})

test_that("run_proteomics_pipeline fails clearly when input files are missing", {

    expect_error(
        ProteinBatcher::run_proteomics_pipeline(
            path_pgmatrix        = "non_existing_pg.tsv",
            path_annotation      = "non_existing_ann.tsv",
            path_output          = tempdir(),
            tests                = "A_vs_B",
            tests_interaction    = "NA",
            formula              = ~ 0 + condition,
            reference_condition  = "A",
            percent_missing      = 50,
            experiment           = "FAIL",
            plots                = FALSE
        )
    )
})
