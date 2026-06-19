# End-to-end test of run_proteomics_pipeline() with plots = TRUE.
# This covers the plotting branches of R/Step6_Wrapper_Workflow.R (volcano
# export for the main/common/interaction effects and the deregulogram) and
# asserts the corrected `effects` return structure: a contrast-keyed list that
# can be subset with `$`, with statistics kept in rowData() under the
# "<stat>_<contrast>" naming convention.

test_that("run_proteomics_pipeline(plots = TRUE) writes plots and returns a contrast-keyed `effects`", {

    path_annotation <- system.file(
        "extdata", "annotation_HaCaT.tsv",
        package = "ProteinBatcher"
    )
    path_pgmatrix <- system.file(
        "extdata",
        "2024MK017_HaCaT_Stimulation_1to24_Astral_report.pg_matrix.tsv",
        package = "ProteinBatcher"
    )
    if (path_annotation == "" || path_pgmatrix == "") {
        testthat::skip("Required extdata files not found.")
    }

    # Dedicated output directory so we can check the files that were written.
    path_output <- file.path(tempdir(), "pb_plots_test")
    dir.create(path_output, showWarnings = FALSE)
    on.exit(unlink(path_output, recursive = TRUE), add = TRUE)

    experiment <- "HaCaT"
    formula <- ~ 0 + condition + batch + condition:batch
    tests <- c("IL13_vs_NoTreated", "IL22_vs_NoTreated", "COMBO_vs_NoTreated")
    tests_interaction <- c("IL13.batchday2", "IL22.batchday2",
                           "COMBO.batchday2")

    set.seed(123)
    res <- ProteinBatcher::run_proteomics_pipeline(
        path_pgmatrix       = path_pgmatrix,
        path_annotation     = path_annotation,
        path_output         = path_output,
        tests               = tests,
        tests_interaction   = tests_interaction,
        formula             = formula,
        reference_condition = "NoTreated",
        percent_missing     = 50,
        ldv_source          = "per-condition",
        experiment          = experiment,
        plots               = TRUE
    )

    # ---- `effects` is always present and keyed by contrast ----
    expect_true("effects" %in% names(res))
    expect_type(res$effects, "list")
    expect_true(all(tests %in% names(res$effects)))

    # A contrast can be accessed with `$` (the reported bug) ...
    one <- res$effects[["IL13_vs_NoTreated"]]
    expect_type(one, "list")
    expect_true(all(c("all_common_effect", "common_effect",
                      "interaction_effect") %in% names(one)))

    # ... and each effect slot is a SummarizedExperiment whose per-contrast
    # statistics live in rowData() under "<stat>_<contrast>".
    se_eff <- one$all_common_effect
    expect_s4_class(se_eff, "SummarizedExperiment")
    rd_cols <- colnames(SummarizedExperiment::rowData(se_eff))
    expect_true("log2FC_IL13_vs_NoTreated" %in% rd_cols)
    expect_true("p.adj_IL13_vs_NoTreated" %in% rd_cols)

    # ---- Plot files were written to disk (side effect) ----
    pdfs <- list.files(path_output, pattern = "\\.pdf$")
    expect_true(length(pdfs) > 0)
    expect_true(any(grepl("MainEffect_AllProteins", pdfs)))
    expect_true(any(grepl("CommonEffect", pdfs)))
    expect_true(any(grepl("InteractionEffect", pdfs)))
    expect_true(any(grepl("Deregulogram", pdfs)))

    # ---- Imputation CSVs were written (base-R writers) ----
    csvs <- list.files(path_output, pattern = "\\.csv$")
    expect_true(any(grepl("Imputation_before", csvs)))
    expect_true(any(grepl("Imputation_after", csvs)))
    expect_true(any(grepl("filtered_out_proteins", csvs)))
})

test_that("run_proteomics_pipeline(plots = TRUE) also works without interaction", {

    path_annotation <- system.file(
        "extdata", "annotation_HaCaT.tsv",
        package = "ProteinBatcher"
    )
    path_pgmatrix <- system.file(
        "extdata",
        "2024MK017_HaCaT_Stimulation_1to24_Astral_report.pg_matrix.tsv",
        package = "ProteinBatcher"
    )
    if (path_annotation == "" || path_pgmatrix == "") {
        testthat::skip("Required extdata files not found.")
    }

    path_output <- file.path(tempdir(), "pb_plots_test_noint")
    dir.create(path_output, showWarnings = FALSE)
    on.exit(unlink(path_output, recursive = TRUE), add = TRUE)

    # No interaction term: tests_interaction = "NA" disables the interaction /
    # deregulogram branch, but main- and common-effect volcanoes still run.
    set.seed(123)
    res <- ProteinBatcher::run_proteomics_pipeline(
        path_pgmatrix       = path_pgmatrix,
        path_annotation     = path_annotation,
        path_output         = path_output,
        tests               = c("IL13_vs_NoTreated"),
        tests_interaction   = c("NA"),
        formula             = ~ 0 + condition + batch,
        reference_condition = "NoTreated",
        percent_missing     = 50,
        ldv_source          = "per-condition",
        experiment          = "HaCaT",
        plots               = TRUE
    )

    expect_true("effects" %in% names(res))
    expect_true("IL13_vs_NoTreated" %in% names(res$effects))
    expect_s4_class(res$effects$IL13_vs_NoTreated$all_common_effect,
                    "SummarizedExperiment")

    pdfs <- list.files(path_output, pattern = "\\.pdf$")
    expect_true(any(grepl("MainEffect_AllProteins", pdfs)))
    # No interaction -> no deregulogram files
    expect_false(any(grepl("Deregulogram", pdfs)))
})
