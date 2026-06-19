# Unit tests for the internal helpers in R/Wrapper_Imputation_helpers.R.
# These use small synthetic inputs (no shipped extdata required) so they run
# everywhere and exercise the error branches and base-R I/O helpers that are
# not reached by the end-to-end pipeline tests.

# ---------------------------------------------------------------------------
# .pp_validate_inputs
# ---------------------------------------------------------------------------
test_that(".pp_validate_inputs validates files, params and annotation header", {
    tmp <- tempfile(fileext = ".tsv")
    ann <- tempfile(fileext = ".tsv")
    writeLines("dummy", tmp)

    good_header <- paste(
        c("file", "sample", "sample_name", "condition",
          "replicate", "batch", "donor_id"),
        collapse = "\t"
    )
    writeLines(c(good_header, "f1\ts1\tA_1\tA\t1\tb1\td1"), ann)

    out <- ProteinBatcher:::.pp_validate_inputs(
        path_pgmatrix = tmp, path_annotation = ann, level = "protein",
        type = "DIA", percent_missing = 50, path_output = tempdir(),
        experiment = "EXP"
    )
    expect_type(out, "list")
    expect_identical(out$experiment, "EXP")
    expect_identical(out$percent_missing, 50)

    # Missing pgmatrix file
    expect_error(
        ProteinBatcher:::.pp_validate_inputs(
            "does_not_exist.tsv", ann, "protein", "DIA", 50,
            tempdir(), "EXP"),
        "path_pgmatrix"
    )
    # Missing annotation file
    expect_error(
        ProteinBatcher:::.pp_validate_inputs(
            tmp, "no_ann.tsv", "protein", "DIA", 50, tempdir(), "EXP"),
        "path_annotation"
    )
    # Bad percent_missing
    expect_error(
        ProteinBatcher:::.pp_validate_inputs(
            tmp, ann, "protein", "DIA", 0, tempdir(), "EXP"),
        "percent_missing"
    )
    expect_error(
        ProteinBatcher:::.pp_validate_inputs(
            tmp, ann, "protein", "DIA", 150, tempdir(), "EXP"),
        "percent_missing"
    )
    # Output directory does not exist
    expect_error(
        ProteinBatcher:::.pp_validate_inputs(
            tmp, ann, "protein", "DIA", 50,
            file.path(tempdir(), "no_such_dir_xyz"), "EXP"),
        "path_output"
    )
    # Empty experiment
    expect_error(
        ProteinBatcher:::.pp_validate_inputs(
            tmp, ann, "protein", "DIA", 50, tempdir(), ""),
        "experiment"
    )
})

test_that(".pp_validate_inputs detects missing and duplicated columns", {
    tmp <- tempfile(fileext = ".tsv"); writeLines("dummy", tmp)

    # Missing required columns
    ann_missing <- tempfile(fileext = ".tsv")
    writeLines(c("file\tsample\tcondition", "f1\ts1\tA"), ann_missing)
    expect_error(
        ProteinBatcher:::.pp_validate_inputs(
            tmp, ann_missing, "protein", "DIA", 50, tempdir(), "EXP"),
        "missing required columns"
    )

    # Duplicated column names
    ann_dup <- tempfile(fileext = ".tsv")
    dup_header <- paste(
        c("file", "sample", "sample_name", "condition",
          "replicate", "batch", "donor_id", "file"),
        collapse = "\t"
    )
    writeLines(c(dup_header, "f1\ts1\tA_1\tA\t1\tb1\td1\tf1"), ann_dup)
    expect_error(
        ProteinBatcher:::.pp_validate_inputs(
            tmp, ann_dup, "protein", "DIA", 50, tempdir(), "EXP"),
        "duplicated column names"
    )
})

# ---------------------------------------------------------------------------
# .pp_write_filtered / .pp_write_before_after  (base-R CSV writers)
# ---------------------------------------------------------------------------
test_that(".pp_write_filtered writes a single CSV", {
    removed <- data.frame(Protein.Group = c("P1", "P2"),
                          reason = c("missing", "missing"),
                          stringsAsFactors = FALSE)
    out <- ProteinBatcher:::.pp_write_filtered(removed, tempdir(), "EXP")
    expect_true(file.exists(out))
    expect_match(out, "EXP_filtered_out_proteins\\.csv$")
    back <- utils::read.csv(out, stringsAsFactors = FALSE)
    expect_equal(nrow(back), 2L)
    expect_true("Protein.Group" %in% colnames(back))
})

test_that(".pp_write_before_after writes two CSVs (before and after)", {
    X <- matrix(c(1, 2, 3, 4), nrow = 2,
                dimnames = list(c("P1", "P2"), c("s1", "s2")))
    se <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(intensity = X),
        rowData = S4Vectors::DataFrame(Genes = c("g1", "g2"))
    )
    paths <- ProteinBatcher:::.pp_write_before_after(se, se, tempdir(), "EXP")
    expect_length(paths, 2L)
    expect_true(all(file.exists(paths)))
    expect_match(paths[["before"]], "EXP_Imputation_before\\.csv$")
    expect_match(paths[["after"]],  "EXP_Imputation_after\\.csv$")

    before <- utils::read.csv(paths[["before"]], stringsAsFactors = FALSE)
    expect_true("gene" %in% colnames(before))
    expect_equal(nrow(before), 2L)
})

# ---------------------------------------------------------------------------
# readExpDesign
# ---------------------------------------------------------------------------
test_that("readExpDesign reads a DIA annotation and applies make.names", {
    ann <- tempfile(fileext = ".tsv")
    df <- data.frame(
        file = c("f1", "f2"),
        sample = c("s1", "s2"),
        sample_name = c("A_1", "B_1"),
        condition = c("No Treated", "IL13"),  # space -> make.names
        replicate = c(1L, 1L),
        stringsAsFactors = FALSE
    )
    utils::write.table(df, ann, sep = "\t", row.names = FALSE, quote = FALSE)

    out <- ProteinBatcher:::readExpDesign(ann, type = "DIA")
    expect_s3_class(out, "data.frame")
    expect_true("label" %in% colnames(out))            # replicate not all NA
    expect_identical(out$condition[1], make.names("No Treated"))
})

test_that("readExpDesign errors on missing columns, duplicates and bad type", {
    base_df <- data.frame(
        file = c("f1", "f2"), sample = c("s1", "s2"),
        sample_name = c("A_1", "B_1"), condition = c("A", "B"),
        replicate = c(1L, 1L), stringsAsFactors = FALSE
    )
    write_tsv <- function(d) {
        p <- tempfile(fileext = ".tsv")
        utils::write.table(d, p, sep = "\t", row.names = FALSE, quote = FALSE)
        p
    }

    # Missing 'file'
    expect_error(
        ProteinBatcher:::readExpDesign(write_tsv(base_df[, -1]), type = "DIA"),
        "'file' column"
    )
    # Duplicated sample_name
    dup_name <- base_df; dup_name$sample_name <- c("A_1", "A_1")
    expect_error(
        ProteinBatcher:::readExpDesign(write_tsv(dup_name), type = "DIA"),
        "Duplicated 'sample_name'"
    )
    # Duplicated sample
    dup_samp <- base_df; dup_samp$sample <- c("s1", "s1")
    expect_error(
        ProteinBatcher:::readExpDesign(write_tsv(dup_samp), type = "DIA"),
        "Duplicated 'sample'"
    )
    # Unsupported type
    expect_error(
        ProteinBatcher:::readExpDesign(write_tsv(base_df), type = "TMT"),
        "not supported"
    )
})

# ---------------------------------------------------------------------------
# readQuantTable
# ---------------------------------------------------------------------------
test_that("readQuantTable reshapes a DIA pg_matrix to wide format", {
    quant <- tempfile(fileext = ".tsv")
    df <- data.frame(
        Protein.Group = c("P1", "P2"),
        Protein.Names = c("n1", "n2"),
        Genes = c("g1", "g2"),
        Run1 = c(10, 20),
        Run2 = c(11, 21),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    utils::write.table(df, quant, sep = "\t", row.names = FALSE, quote = FALSE)

    out <- ProteinBatcher:::readQuantTable(quant, type = "DIA")
    expect_s3_class(out, "data.frame")
    expect_true(all(c("Protein.Group", "Protein.Names", "Genes") %in%
                        colnames(out)))
    expect_true(all(c("Run1", "Run2") %in% colnames(out)))
    expect_equal(nrow(out), 2L)
})

test_that("readQuantTable rejects unsupported types", {
    quant <- tempfile(fileext = ".tsv")
    writeLines(c("Protein.Group\tRun1", "P1\t10"), quant)
    expect_error(
        ProteinBatcher:::readQuantTable(quant, type = "TMT"),
        "not supported"
    )
})

# ---------------------------------------------------------------------------
# make_unique
# ---------------------------------------------------------------------------
test_that("make_unique builds unique name/ID columns from base data.frame", {
    prot <- data.frame(
        Protein.Group = c("P1", "P2", "P3"),
        Genes = c("g1", NA, ""),         # NA and "" must fall back to ID
        stringsAsFactors = FALSE
    )
    out <- ProteinBatcher:::make_unique(prot, names = "Genes",
                                        ids = "Protein.Group")
    expect_true(all(c("name", "ID") %in% colnames(out)))
    expect_identical(out$ID, c("P1", "P2", "P3"))
    # 'g1' kept; NA -> "P2"; "" -> "P3"
    expect_identical(out$name, c("g1", "P2", "P3"))

    # make.unique should disambiguate repeated fallbacks
    prot2 <- data.frame(
        Protein.Group = c("P1", "P1"),
        Genes = c(NA, NA),
        stringsAsFactors = FALSE
    )
    # both NA in Genes but IDs present -> not "double NA"
    out2 <- ProteinBatcher:::make_unique(prot2, "Genes", "Protein.Group")
    expect_identical(out2$name, c("P1", "P1.1"))
})

test_that("make_unique errors on missing columns and double NAs", {
    prot <- data.frame(Protein.Group = "P1", Genes = "g1",
                       stringsAsFactors = FALSE)
    expect_error(
        ProteinBatcher:::make_unique(prot, "NotThere", "Protein.Group"),
        "is not a column"
    )
    expect_error(
        ProteinBatcher:::make_unique(prot, "Genes", "NotThere"),
        "is not a column"
    )
    prot_na <- data.frame(Protein.Group = NA, Genes = NA,
                          stringsAsFactors = FALSE)
    expect_error(
        ProteinBatcher:::make_unique(prot_na, "Genes", "Protein.Group"),
        "NAs in both"
    )
})

test_that("make_unique accepts tbl_df-classed input", {
    # Simulate a tibble without depending on the tibble package: the function
    # branch only checks inherits(x, "tbl_df").
    prot <- data.frame(Protein.Group = c("P1", "P2"),
                       Genes = c("g1", "g2"), stringsAsFactors = FALSE)
    class(prot) <- c("tbl_df", "tbl", "data.frame")
    out <- ProteinBatcher:::make_unique(prot, "Genes", "Protein.Group")
    expect_s3_class(out, "data.frame")
    expect_identical(out$name, c("g1", "g2"))
})

# ---------------------------------------------------------------------------
# make_se_customized
# ---------------------------------------------------------------------------
make_demo_inputs <- function() {
    prot <- data.frame(
        Protein.Group = c("P1", "P2"),
        Protein.Names = c("n1", "n2"),
        Genes = c("g1", "g2"),
        s1 = c(10, 20),
        s2 = c(11, 21),
        name = c("g1", "g2"),
        ID = c("P1", "P2"),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    cols <- which(colnames(prot) %in% c("s1", "s2"))
    exp <- data.frame(
        label = c("s1", "s2"),
        sample_name = c("A_1", "B_1"),
        sample = c("A", "B"),
        condition = c("A", "B"),
        replicate = c(1L, 1L),
        stringsAsFactors = FALSE
    )
    list(prot = prot, cols = cols, exp = exp)
}

test_that("make_se_customized builds a SummarizedExperiment with named assay", {
    d <- make_demo_inputs()
    se <- ProteinBatcher:::make_se_customized(
        d$prot, d$cols, d$exp, exp = "LFQ"
    )
    expect_s4_class(se, "SummarizedExperiment")
    # Assay must be named 'intensity' (the bug fix for assays(1): '')
    expect_identical(SummarizedExperiment::assayNames(se), "intensity")
    expect_equal(nrow(se), 2L)
    expect_equal(ncol(se), 2L)
    expect_true(all(c("name", "ID") %in%
                        colnames(SummarizedExperiment::rowData(se))))
})

test_that("make_se_customized validates its inputs", {
    d <- make_demo_inputs()

    # Missing name/ID in proteins
    bad_prot <- d$prot[, setdiff(colnames(d$prot), c("name", "ID"))]
    cols_bad <- which(colnames(bad_prot) %in% c("s1", "s2"))
    expect_error(
        ProteinBatcher:::make_se_customized(bad_prot, cols_bad, d$exp),
        "'name' and/or 'ID'"
    )

    # Missing condition in expdesign
    bad_exp <- d$exp[, setdiff(colnames(d$exp), "condition")]
    expect_error(
        ProteinBatcher:::make_se_customized(d$prot, d$cols, bad_exp),
        "experimental design"
    )

    # Non-numeric 'columns'
    prot_chr <- d$prot
    prot_chr$s1 <- as.character(prot_chr$s1)
    expect_error(
        ProteinBatcher:::make_se_customized(prot_chr, d$cols, d$exp),
        "should be numeric"
    )
})

# ---------------------------------------------------------------------------
# Lower-level imputation helpers (edge branches not hit by impute_se defaults)
# ---------------------------------------------------------------------------
test_that(".impute_exclude_na_condition drops NA-condition samples", {
    X <- matrix(1:6, nrow = 2,
                dimnames = list(c("p1", "p2"), c("s1", "s2", "s3")))
    cd <- S4Vectors::DataFrame(condition = c("A", NA, "B"))
    expect_message(
        res <- ProteinBatcher:::.impute_exclude_na_condition(X, cd),
        "Excluding 1 sample"
    )
    expect_true(res$dropped)
    expect_equal(ncol(res$X), 2L)
    expect_identical(res$cond, c("A", "B"))
})

test_that(".impute_sd_protein uses a robust fallback for degenerate SDs", {
    # Row 2 has no variance -> sd 0 -> must be replaced by fallback
    X <- matrix(c(1, 5, 3, 5, 7, 5), nrow = 2,
                dimnames = list(c("p1", "p2"), c("s1", "s2", "s3")))
    out <- ProteinBatcher:::.impute_sd_protein(X)
    expect_length(out$sd, 2L)
    expect_true(all(is.finite(out$sd)))
    expect_true(all(out$sd > 0))
    # Explicit fallback is honoured
    out2 <- ProteinBatcher:::.impute_sd_protein(X, sd_fallback = 0.25)
    expect_equal(unname(out2$sd[2]), 0.25)
})

test_that(".impute_ldv_fill respects the lower bound and records the map", {
    X1 <- matrix(c(10, NA, 12, 13), nrow = 2, byrow = TRUE,
                 dimnames = list(c("p1", "p2"), c("s1", "s2")))
    cond <- c("A", "A")
    ldv <- list(type = "global", mu = -100)        # very low mean
    sd_prot <- c(0.1, 0.1)
    imp_map <- matrix("none", nrow = 2, ncol = 2,
                      dimnames = dimnames(X1))
    set.seed(1)
    out <- ProteinBatcher:::.impute_ldv_fill(
        X1, cond, ldv, sd_prot, "global",
        lower_bound = 0, imp_map = imp_map
    )
    expect_false(anyNA(out$X))
    expect_true(all(out$X >= 0))                    # lower_bound enforced
    expect_identical(out$map["p1", "s2"], "ldv")
})

test_that(".impute_write_back expands the map when samples were dropped", {
    X0 <- matrix(1:6, nrow = 2,
                 dimnames = list(c("p1", "p2"), c("s1", "s2", "s3")))
    se <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(intensity = X0),
        colData = S4Vectors::DataFrame(condition = c("A", "A", "B"))
    )
    # Imputed sub-matrix covers only s1 and s3 (s2 was dropped)
    X_imp <- X0[, c("s1", "s3"), drop = FALSE]
    imp_map <- matrix("mean", nrow = 2, ncol = 2,
                      dimnames = list(c("p1", "p2"), c("s1", "s3")))
    out <- ProteinBatcher:::.impute_write_back(se, X0, X_imp, imp_map,
                                               dropped = TRUE)
    md <- S4Vectors::metadata(out)$imputation_map
    expect_equal(ncol(md), 3L)                      # expanded back to 3 cols
    expect_identical(md["p1", "s2"], "none")        # dropped col stays 'none'
    expect_true("imputed" %in%
                    colnames(SummarizedExperiment::rowData(out)))
})
