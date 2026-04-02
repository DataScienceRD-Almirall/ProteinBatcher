test_that("plot_deregulogram writes PDFs for valid contrasts and interactions", {

    ## --- synthetic assay ---
    X <- matrix(
        rnorm(3 * 4),
        nrow = 3,
        dimnames = list(paste0("p", 1:3), paste0("s", 1:4))
    )

    ## --- colData (not used directly but realistic) ---
    cd <- S4Vectors::DataFrame(
        condition = c("A","A","B","B"),
        sex       = c("female","male","female","male"),
        batch     = c("d1","d1","d2","d2"),
        replicate = c("r1","r2","r1","r2"),
        donor_id  = c("1","1","1","1"),
        label     = colnames(X),
        row.names = colnames(X)
    )

    ## --- contrasts ---
    main_contrast  <- "B_vs_A"
    inter_contrast <- "B_vs_A.sexmale"

    ## --- rowData with CURRENT VALID scheme ---
    rd <- S4Vectors::DataFrame(
        Genes = c("G1","G2","G3"),
        ID    = c("P1","P2","P3"),
        row.names = rownames(X)
    )

    ## main effects
    rd[[paste0("log2FC_", main_contrast)]] <- c( 1.5, 0.2, -1.2)
    rd[[paste0("p.adj_",  main_contrast)]] <- c( 0.01, 0.8,  0.03)

    ## interaction effects
    rd[[paste0("log2FC_", inter_contrast)]] <- c( 0.8, -0.1,  1.0)
    rd[[paste0("p.adj_",  inter_contrast)]] <- c( 0.02, 0.9,  0.01)

    se <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(intensity = X),
        colData = cd,
        rowData = rd
    )

    tmp <- tempfile()
    dir.create(tmp)

    ## --- run ---
    res <- plot_deregulogram(
        se_limma        = se,
        tests           = main_contrast,
        tests_interaction = inter_contrast,
        factor_name     = "sex",
        factor_levels   = c("female","male"),
        path_output     = tmp,
        experiment      = "TESTEXP",
        alpha           = 0.05,
        lfc             = 1,
        label_col       = "Genes"
    )

    ## --- expected PDF ---
    pdf_file <- file.path(
        tmp,
        paste0("TESTEXP_Deregulogram_", main_contrast, "_sex.pdf")
    )

    testthat::expect_true(file.exists(pdf_file))
    testthat::expect_null(res)
})
