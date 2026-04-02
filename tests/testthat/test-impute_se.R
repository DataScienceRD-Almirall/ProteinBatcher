test_that("impute_se validates inputs", {
    X <- matrix(c(1, NA, 2, 3), nrow = 2, byrow = TRUE)
    se <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(intensity = X),
        colData = S4Vectors::DataFrame(condition = c("A", "A"))
    )

    expect_error(ProteinBatcher::impute_se(list(), threshold = 0.5))
    expect_error(ProteinBatcher::impute_se(se, threshold = 0))      # must be in (0,1]
    expect_error(ProteinBatcher::impute_se(se, threshold = 1.1))
    expect_error(ProteinBatcher::impute_se(se, ldv_source = "x"))   # invalid choice

    se_no_cond <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(intensity = X),
        colData = S4Vectors::DataFrame(group = c("A", "A"))
    )
    expect_error(ProteinBatcher::impute_se(se_no_cond, threshold = 0.5))

    se_chr <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(intensity = matrix("x", nrow = 1, ncol = 2)),
        colData = S4Vectors::DataFrame(condition = c("A", "A"))
    )
    expect_error(ProteinBatcher::impute_se(se_chr, threshold = 0.5))
})

test_that("impute_se imputes NAs and is reproducible with seed", {
    # Two conditions, 2 samples each; design encourages both mean and LDV paths
    X <- matrix(
        c(
            10, NA, 12, 13,   # row1: one NA in A (1/2 missing = 0.5) -> at threshold=0.75 mean-impute in A
            NA, NA,  2,  3,   # row2: A all missing -> LDV for A; B ok
            5,  6,  NA, NA    # row3: B all missing -> LDV for B
        ),
        nrow = 3, byrow = TRUE,
        dimnames = list(paste0("p", 1:3), paste0("s", 1:4))
    )
    se <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(intensity = X),
        colData = S4Vectors::DataFrame(condition = c("A", "A", "B", "B")),
        rowData = S4Vectors::DataFrame(Protein.Names = c("P1","P2","P3"))
    )

    set.seed(123)
    se1 <- ProteinBatcher::impute_se(se, threshold = 0.75, ldv_source = "global")
    se2 <- ProteinBatcher::impute_se(se, threshold = 0.75, ldv_source = "global")

    X1 <- SummarizedExperiment::assay(se1)
    X2 <- SummarizedExperiment::assay(se2)

    # No missing values after imputation for valid samples
    expect_false(anyNA(X1))

    # Metadata: imputation_map recorded; rowData has "imputed" flag (as per docs)
    md <- S4Vectors::metadata(se1)
    expect_true("imputation_map" %in% names(md))
    expect_true(is.matrix(md$imputation_map))
    expect_true(all(md$imputation_map %in% c("none", "mean", "ldv")))

    rd <- SummarizedExperiment::rowData(se1)
    expect_true("imputed" %in% colnames(rd))
})
