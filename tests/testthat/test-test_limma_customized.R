test_that("test_limma_customized runs on extdata (main + interaction) and writes results to rowData", {

    # Skip if limma is not available for some reason
    if (!requireNamespace("limma", quietly = TRUE)) {
        testthat::skip("limma not installed.")
    }

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

    if (path_annotation == "" || path_pgmatrix == "") {
        testthat::skip("Required extdata files not found.")
    }

    # ---- Build SE via the same import path used by the pipeline ----
    se_raw <- ProteinBatcher:::make_se_from_files(
        quant_table_pat = path_pgmatrix,
        exp_anno_path   = path_annotation,
        type            = "DIA",
        level           = "protein",
        log2transform   = TRUE
    )

    expect_s4_class(se_raw, "SummarizedExperiment")
    expect_true("condition" %in% colnames(SummarizedExperiment::colData(se_raw)))

    # ---- Filter and impute (deterministic) ----
    filt <- ProteinBatcher::filter_se_missing(se_raw, percentage = 50)

    se_imp <- ProteinBatcher::impute_se(
        filt$se_filt,
        threshold  = 0.5,
        ldv_source = "per-condition"
    )

    # Must have imputation_map for the imputation-aware FDR logic
    md <- S4Vectors::metadata(se_imp)
    expect_true("imputation_map" %in% names(md))
    expect_true(is.matrix(md$imputation_map))

    # ---- Limma parameters (as in vignette) ----
    tests <- c("IL13_vs_NoTreated", "IL22_vs_NoTreated", "COMBO_vs_NoTreated")
    tests_interaction <- c("IL13.batchday2", "IL22.batchday2", "COMBO.batchday2")
    design_formula <- stats::formula(~ 0 + condition + batch + condition:batch)
    ref_condition <- "NoTreated"

    # ---- Run: main-only first (no interaction) ----
    se_main <- ProteinBatcher::test_limma_customized(
        se = se_imp,
        type = "manual",
        test = tests,
        test_interaction = "NA",
        design_formula = design_formula,
        ref_condition  = ref_condition,
        paired = FALSE,
        block_effect = FALSE
    )

    expect_s4_class(se_main, "SummarizedExperiment")

    rd_main <- SummarizedExperiment::rowData(se_main)
    expect_true(ncol(rd_main) >= ncol(SummarizedExperiment::rowData(se_imp)))

    # We don't assume exact column naming, but we do require:
    # - new numeric columns exist
    new_cols <- setdiff(colnames(rd_main), colnames(SummarizedExperiment::rowData(se_imp)))
    expect_true(length(new_cols) > 0)

    is_num <- vapply(rd_main[, new_cols, drop = FALSE], is.numeric, logical(1))
    expect_true(any(is_num))

    # Heuristic check: for at least one test, columns mention the contrast name
    # (robust to different separator conventions)
    any_il13 <- any(grepl("IL13", colnames(rd_main), fixed = FALSE))
    expect_true(any_il13)

    # ---- Run: with interaction ----
    se_int <- ProteinBatcher::test_limma_customized(
        se = se_imp,
        type = "manual",
        test = tests,
        test_interaction = tests_interaction,
        design_formula = design_formula,
        ref_condition  = ref_condition,
        paired = FALSE,
        block_effect = FALSE
    )

    rd_int <- SummarizedExperiment::rowData(se_int)

    # Should have at least as many columns as main-only run
    expect_true(ncol(rd_int) >= ncol(rd_main))

    # Heuristic: interaction-related naming typically includes "batchday2" somewhere
    # (again robust to exact formatting)
    expect_true(any(grepl("batchday2", colnames(rd_int), ignore.case = TRUE)))

    # Basic sanity on adjusted p-values: find any p.adj-like columns and ensure in [0,1] where finite
    padj_cols <- grep("p\\.?adj|adj\\.p|FDR", colnames(rd_int), value = TRUE, ignore.case = TRUE)
    expect_true(length(padj_cols) > 0)

    vals <- unlist(as.data.frame(rd_int[, padj_cols, drop = FALSE]), use.names = FALSE)
    vals <- vals[is.finite(vals)]
    # If there are finite values, they must lie between 0 and 1
    if (length(vals) > 0) {
        expect_true(all(vals >= 0 & vals <= 1))
    }
})

test_that("block_effect=TRUE requires a block column (or errors cleanly)", {

    if (!requireNamespace("limma", quietly = TRUE)) {
        testthat::skip("limma not installed.")
    }

    X <- matrix(rnorm(2 * 4), nrow = 2)
    colnames(X) <- paste0("s", 1:4)
    rownames(X) <- paste0("p", 1:2)

    cd <- S4Vectors::DataFrame(
        condition  = c("A", "A", "B", "B"),
        batch      = c("d1", "d1", "d1", "d1"),
        replicate  = c("r1", "r2", "r1", "r2")
        # block intentionally absent here
    )

    se <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(intensity = X),
        colData = cd
    )

    imp_map <- matrix("none", nrow = nrow(X), ncol = ncol(X),
                      dimnames = list(rownames(X), colnames(X)))
    S4Vectors::metadata(se)$imputation_map <- imp_map

    expect_error(
        ProteinBatcher::test_limma_customized(
            se = se,
            test = "B_vs_A",
            test_interaction = "NA",
            design_formula = ~ 0 + condition + batch,
            ref_condition = "A",
            paired = FALSE,
            block_effect = TRUE
        )
    )
})

test_that("LDV/LDV features are not spuriously significant after imputation-aware FDR", {

    if (!requireNamespace("limma", quietly = TRUE)) {
        testthat::skip("limma not installed.")
    }

    X <- matrix(
        c(
            10, 10, 20, 20,   # p1: true signal
            1,  1,  1,  1     # p2: LDV/LDV (should not be significant)
        ),
        nrow = 2, byrow = TRUE,
        dimnames = list(c("p1", "p2"), paste0("s", 1:4))
    )

    cd <- S4Vectors::DataFrame(
        condition = c("A", "A", "B", "B"),
        batch     = c("d1", "d2", "d1", "d2"),
        replicate = c("r1", "r2", "r1", "r2"),
        donor_id  = c("1", "1", "1", "1"),
        label     = paste0("s", 1:4)
    )

    rd <- S4Vectors::DataFrame(
        name  = rownames(X),
        ID    = rownames(X),
        Genes = c("G1", "G2")
    )
    rownames(rd) <- rownames(X)

    se <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(intensity = X),
        colData = cd,
        rowData = rd
    )

    imp_map <- matrix("none", nrow = 2, ncol = 4,
                      dimnames = list(rownames(X), colnames(X)))
    imp_map["p2", ] <- "ldv"
    S4Vectors::metadata(se)$imputation_map <- imp_map

    se_out <- ProteinBatcher::test_limma_customized(
        se = se,
        test = "B_vs_A",
        test_interaction = "NA",
        design_formula = ~ 0 + condition,
        ref_condition = "A",
        paired = FALSE,
        block_effect = FALSE
    )

    rd_out <- SummarizedExperiment::rowData(se_out)
    padj_cols <- grep("p\\.?adj|adj\\.p|FDR", colnames(rd_out),
                      value = TRUE, ignore.case = TRUE)
    testthat::expect_true(length(padj_cols) >= 1)

    p2_padj <- as.numeric(rd_out["p2", padj_cols[1]])
    testthat::expect_true(is.na(p2_padj) || p2_padj >= 0.05)
})
