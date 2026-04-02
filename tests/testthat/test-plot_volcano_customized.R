test_that("plot_volcano_customized detects valid contrasts and computes significance", {
    ## --- synthetic assay ---
    X <- matrix(
        rnorm(3 * 4),
        nrow = 3,
        dimnames = list(c("p1","p2","p3"), paste0("s", 1:4))
    )
    ## --- required colData ---
    cd <- S4Vectors::DataFrame(
        condition = c("A","A","B","B"),
        batch     = c("d1","d2","d1","d2"),
        replicate = c("r1","r2","r1","r2"),
        donor_id  = c("1","1","1","1"),
        label     = colnames(X),
        row.names = colnames(X)
    )
    ## --- rowData with REQUIRED naming (log2FC_*) ---
    rd <- S4Vectors::DataFrame(
        ID      = c("Prot1","Prot2","Prot3"),
        imputed = c(FALSE, TRUE, FALSE),
        row.names = rownames(X)
    )
    contrast <- "B_vs_A"
    rd[[paste0("log2FC_", contrast)]] <- c(2.0, 0.5, -1.5)
    rd[[paste0("p.adj_",  contrast)]] <- c(0.01, 0.9, 0.03)

    se <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(intensity = X),
        colData = cd,
        rowData = rd
    )

    ## --- plot = FALSE -> data.frame ---
    df <- ProteinBatcher::plot_volcano_customized(
        dep = se,
        contrast = contrast,
        plot = FALSE,
        alpha = 0.05,
        lfc   = 1
    )

    testthat::expect_s3_class(df, "data.frame")
    ## --- plot = TRUE -> ggplot ---
    p <- ProteinBatcher::plot_volcano_customized(
        dep = se,
        contrast = contrast,
        plot = TRUE,
        add_names = FALSE
    )

    testthat::expect_s3_class(p, "ggplot")
})

