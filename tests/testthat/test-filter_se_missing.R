test_that("filter_se_missing validates inputs", {
    # Not a SummarizedExperiment
    expect_error(
        ProteinBatcher::filter_se_missing(se = list(), percentage = 50),
        "`se` must be a SummarizedExperiment"
    )

    # percentage must be numeric scalar
    se0 <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(1, nrow = 1, ncol = 2)),
        colData = S4Vectors::DataFrame(condition = c("A", "B"))
    )
    expect_error(ProteinBatcher::filter_se_missing(se0, percentage = NA_real_))
    expect_error(ProteinBatcher::filter_se_missing(se0, percentage = c(1, 2)))
    expect_error(ProteinBatcher::filter_se_missing(se0, percentage = -1))
    expect_error(ProteinBatcher::filter_se_missing(se0, percentage = 101))

    # missing condition in colData
    se_no_cond <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(1, nrow = 1, ncol = 2)),
        colData = S4Vectors::DataFrame(group = c("A", "B"))
    )
    expect_error(
        ProteinBatcher::filter_se_missing(se_no_cond, percentage = 50),
        "colData\\(se\\)\\$condition"
    )

    # assay must be numeric
    se_chr <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix("x", nrow = 1, ncol = 2)),
        colData = S4Vectors::DataFrame(condition = c("A", "B"))
    )
})

test_that("filter_se_missing filters by per-condition non-missing proportion", {
    # Construct a minimal SE:
    # 2 conditions (A,B), each with 2 samples
    X <- matrix(
        c(
            1, NA,  5,  6,   # row1: A has 1/2 observed (=0.5), B has 2/2 (=1.0) -> keep at 50+
            NA, NA, 7, NA,   # row2: A 0/2 (=0.0), B 1/2 (=0.5) -> keep at 50, drop at 60
            NA, NA, NA, NA   # row3: all missing -> drop always
        ),
        nrow = 3, byrow = TRUE,
        dimnames = list(paste0("p", 1:3), paste0("s", 1:4))
    )

    se <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(intensity = X),
        colData = S4Vectors::DataFrame(condition = c("A", "A", "B", "B")),
        rowData = S4Vectors::DataFrame(
            Protein.Names = c("P1", "P2", "P3"),
            Genes         = c("G1", "G2", "G3")
        )
    )

    out50 <- ProteinBatcher::filter_se_missing(se, percentage = 50)
    expect_true(is.list(out50))
    expect_true(all(c("se_filt", "removed") %in% names(out50)))
    expect_s4_class(out50$se_filt, "SummarizedExperiment")
    expect_true(is.data.frame(out50$removed))
    expect_true(all(c("protein_id", "gene_name", "reason") %in% colnames(out50$removed)))

    # At 50%: keep row1 & row2, drop row3
    expect_equal(nrow(out50$se_filt), 2)
    expect_true("P3" %in% out50$removed$protein_id)

    out60 <- ProteinBatcher::filter_se_missing(se, percentage = 60)
    # At 60%: row1 keep (B fully observed), row2 drop (best is 0.5), row3 drop
    expect_equal(nrow(out60$se_filt), 1)
    expect_true(all(c("P2", "P3") %in% out60$removed$protein_id))
})

test_that("filter_se_missing fills removed annotations safely when rowData lacks fields", {
    X <- matrix(c(1, NA, NA, NA), nrow = 2, byrow = TRUE)
    se <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(intensity = X),
        colData = S4Vectors::DataFrame(condition = c("A", "A"))
        # rowData intentionally missing Protein.Names / Genes
    )

    out <- ProteinBatcher::filter_se_missing(se, percentage = 100)
    expect_true(is.data.frame(out$removed))
    # These columns exist and are character (may contain NA)
    expect_true(is.character(out$removed$protein_id))
    expect_true(is.character(out$removed$gene_name))
})
