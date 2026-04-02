test_that(".pp_export_volcano_tables exports PDFs and rowData tables correctly", {

    ## --- synthetic assay ---
    X <- matrix(
        rnorm(3 * 4),
        nrow = 3,
        dimnames = list(paste0("p", 1:3), paste0("s", 1:4))
    )

    ## --- colData (required but not used here) ---
    cd <- S4Vectors::DataFrame(
        condition = c("A","A","B","B"),
        batch     = c("d1","d2","d1","d2"),
        replicate = c("r1","r2","r1","r2"),
        donor_id  = c("1","1","1","1"),
        label     = colnames(X),
        row.names = colnames(X)
    )

    ## --- rowData with CURRENT VALID volcano scheme ---
    contrast <- "B_vs_A"

    rd <- S4Vectors::DataFrame(
        Genes   = c("G1","G2","G3"),
        ID      = c("P1","P2","P3"),
        imputed = c(FALSE, TRUE, FALSE),
        row.names = rownames(X)
    )

    rd[[paste0("log2FC_", contrast)]] <- c(2.0, 0.4, -1.8)
    rd[[paste0("p.adj_",  contrast)]] <- c(0.01, 0.8, 0.03)

    se <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(intensity = X),
        colData = cd,
        rowData = rd
    )

    ## --- effects object as produced by .pp_organize_effects() ---
    effects <- list(
        B_vs_A = list(
            all_common_effect = se,
            common_effect     = se,
            interaction_effect = se[0, ]  # empty
        )
    )

    tmp <- tempfile()
    dir.create(tmp)

    ## --- run export ---
    res <- ProteinBatcher:::.pp_export_volcano_tables(
        effects = effects,
        tests = "B_vs_A",
        path_output = tmp,
        experiment = "TESTEXP",
        effect_slot = "all_common_effect",
        file_tag = "MainEffect_AllProteins"
    )

    ## --- PDF created ---
    pdf_file <- file.path(
        tmp,
        "TESTEXP_MainEffect_AllProteins_B_vs_A.pdf"
    )
    testthat::expect_true(file.exists(pdf_file))

    ## --- returned structure ---
    testthat::expect_type(res, "list")
    testthat::expect_length(res, 1)

    tab <- res[[1]][["B_vs_A"]]
    testthat::expect_s3_class(tab, "data.frame")
    testthat::expect_equal(nrow(tab), 3)
    testthat::expect_true(all(colnames(rd) %in% colnames(tab)))
})


test_that(".pp_export_volcano_tables skips empty effect objects silently", {

    se_empty <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(intensity = matrix(numeric(), nrow = 0, ncol = 0)),
        colData = S4Vectors::DataFrame(),
        rowData = S4Vectors::DataFrame()
    )

    effects <- list(
        B_vs_A = list(
            all_common_effect = se_empty,
            common_effect     = se_empty,
            interaction_effect = se_empty
        )
    )

    tmp <- tempfile()
    dir.create(tmp)

    res <- ProteinBatcher:::.pp_export_volcano_tables(
        effects = effects,
        tests = "B_vs_A",
        path_output = tmp,
        experiment = "TESTEXP",
        effect_slot = "common_effect"
    )

    testthat::expect_null(res[[1]])
})


test_that(".pp_export_volcano_tables errors on invalid plot_contrasts length", {

    effects <- list(
        A = list(
            common_effect = NULL,
            all_common_effect = NULL,
            interaction_effect = NULL
        )
    )

    testthat::expect_error(
        ProteinBatcher:::.pp_export_volcano_tables(
            effects = effects,
            tests = c("A","B"),
            plot_contrasts = c("X","Y","Z"),
            path_output = tempdir(),
            experiment = "TEST"
        ),
        "plot_contrasts must be length 1 or length\\(tests\\)"
    )
})
