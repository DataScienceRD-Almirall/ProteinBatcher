test_that(".pp_organize_effects splits common vs interaction effects correctly", {

    # --- Minimal SummarizedExperiment resembling output of test_limma_customized() ---
    X <- matrix(rnorm(3 * 4), nrow = 3,
                dimnames = list(paste0("p", 1:3), paste0("s", 1:4)))

    cd <- S4Vectors::DataFrame(
        condition = c("A", "A", "B", "B"),
        batch     = c("d1", "d2", "d1", "d2"),
        replicate = c("r1", "r2", "r1", "r2"),
        donor_id  = c("1", "1", "1", "1"),
        label     = colnames(X)
    )
    rownames(cd) <- colnames(X)

    # Base annotation columns (some of base_cols)
    rd <- S4Vectors::DataFrame(
        Protein.Group = c("PG1", "PG2", "PG3"),
        Protein.Names = c("P1", "P2", "P3"),
        Genes         = c("G1", "G2", "G3"),
        name          = rownames(X),
        ID            = rownames(X),
        imputed       = c(FALSE, TRUE, FALSE)
    )
    rownames(rd) <- rownames(X)

    # Add limma-like stats columns using the NEW naming style required by .pp_organize_effects
    # main test prefix:
    main <- "B_vs_A"
    rd[[paste0("diff_",  main)]] <- c( 1.0,  0.2, -0.5)
    rd[[paste0("p.val_", main)]] <- c( 0.01, 0.50, 0.20)
    rd[[paste0("p.adj_", main)]] <- c( 0.02, 0.80, 0.30)

    # interaction test prefix:
    inter <- "B.batchd2"
    rd[[paste0("diff_",  inter)]] <- c(0.1,  0.2,  0.3)
    rd[[paste0("p.val_", inter)]] <- c(0.2,  0.01, 0.9)

    # Make only protein2 significant for interaction at alpha=0.05
    rd[[paste0("p.adj_", inter)]] <- c(0.20, 0.01, 0.60)

    se_limma <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(intensity = X),
        colData = cd,
        rowData = rd
    )

    out <- ProteinBatcher:::.pp_organize_effects(
        se_limma = se_limma,
        tests = main,
        tests_interaction = inter,
        alpha = 0.05
    )

    # Structure
    testthat::expect_true(is.list(out))
    testthat::expect_true(main %in% names(out))
    testthat::expect_true(all(c("common_effect","interaction_effect","all_common_effect") %in% names(out[[main]])))

    common_se <- out[[main]]$common_effect
    inter_se  <- out[[main]]$interaction_effect
    all_se    <- out[[main]]$all_common_effect

    testthat::expect_s4_class(common_se, "SummarizedExperiment")
    testthat::expect_s4_class(inter_se,  "SummarizedExperiment")
    testthat::expect_s4_class(all_se,    "SummarizedExperiment")

    # all_common_effect keeps all proteins
    testthat::expect_equal(nrow(all_se), 3)

    # interaction significant only for p2 -> interaction_effect has 1 row (p2)
    testthat::expect_equal(rownames(inter_se), "p2")
    testthat::expect_equal(nrow(inter_se), 1)

    # common_effect are the non-interacting proteins -> p1 and p3
    testthat::expect_equal(rownames(common_se), c("p1","p3"))
    testthat::expect_equal(nrow(common_se), 2)

    # Column selection: should keep base cols (subset) + stats for main/inter
    rd_all <- SummarizedExperiment::rowData(all_se)
    testthat::expect_true(all(c("Genes","name","ID","imputed") %in% colnames(rd_all)))
    testthat::expect_true(paste0("p.adj_", main) %in% colnames(rd_all))

    rd_int <- SummarizedExperiment::rowData(inter_se)
    testthat::expect_true(paste0("p.adj_", inter) %in% colnames(rd_int))
})

test_that(".pp_organize_effects handles tests_interaction = 'NA' (no interaction)", {

    X <- matrix(rnorm(2 * 2), nrow = 2,
                dimnames = list(c("p1","p2"), c("s1","s2")))

    cd <- S4Vectors::DataFrame(
        condition = c("A","B"),
        batch     = c("d1","d2"),
        replicate = c("r1","r1"),
        donor_id  = c("1","1"),
        label     = colnames(X)
    )
    rownames(cd) <- colnames(X)

    rd <- S4Vectors::DataFrame(
        Genes = c("G1","G2"),
        name  = rownames(X),
        ID    = rownames(X),
        imputed = c(FALSE, FALSE)
    )
    rownames(rd) <- rownames(X)

    main <- "B_vs_A"
    rd[[paste0("p.adj_", main)]] <- c(0.1, 0.2)
    rd[[paste0("diff_",  main)]] <- c(1, -1)

    se_limma <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(intensity = X),
        colData = cd,
        rowData = rd
    )

    out <- ProteinBatcher:::.pp_organize_effects(
        se_limma = se_limma,
        tests = main,
        tests_interaction = "NA",
        alpha = 0.05
    )

    all_se <- out[[main]]$all_common_effect
    common <- out[[main]]$common_effect
    inter  <- out[[main]]$interaction_effect

    # all_common == common when no interaction is tested
    testthat::expect_equal(nrow(all_se), nrow(common))
    testthat::expect_equal(rownames(all_se), rownames(common))

    # interaction_effect should be empty
    testthat::expect_equal(nrow(inter), 0)
})

test_that(".pp_organize_effects recycles tests_interaction of length 1", {

    X <- matrix(rnorm(3 * 2), nrow = 3,
                dimnames = list(paste0("p", 1:3), paste0("s", 1:2)))

    cd <- S4Vectors::DataFrame(
        condition = c("A","B"),
        batch     = c("d1","d2"),
        replicate = c("r1","r1"),
        donor_id  = c("1","1"),
        label     = colnames(X)
    )
    rownames(cd) <- colnames(X)

    rd <- S4Vectors::DataFrame(
        Genes = c("G1","G2","G3"),
        name  = rownames(X),
        ID    = rownames(X),
        imputed = c(FALSE, TRUE, FALSE)
    )
    rownames(rd) <- rownames(X)

    tests <- c("T1_vs_T0", "T2_vs_T0")
    inter <- "T1.batchd2"

    # required main p.adj columns
    rd[[paste0("p.adj_", tests[1])]] <- c(0.1, 0.2, 0.3)
    rd[[paste0("p.adj_", tests[2])]] <- c(0.2, 0.3, 0.4)

    # required interaction p.adj column (only for the single inter pref)
    rd[[paste0("p.adj_", inter)]] <- c(0.01, 0.2, 0.6)

    se_limma <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(intensity = X),
        colData = cd,
        rowData = rd
    )

    out <- ProteinBatcher:::.pp_organize_effects(
        se_limma = se_limma,
        tests = tests,
        tests_interaction = inter,  # length 1 -> should recycle
        alpha = 0.05
    )

    testthat::expect_equal(names(out), tests)
    testthat::expect_true(all(vapply(out, is.list, logical(1))))
})

test_that(".pp_organize_effects errors on missing required columns or wrong lengths", {

    X <- matrix(rnorm(2 * 2), nrow = 2,
                dimnames = list(c("p1","p2"), c("s1","s2")))
    cd <- S4Vectors::DataFrame(
        condition = c("A","B"),
        batch     = c("d1","d2"),
        replicate = c("r1","r1"),
        donor_id  = c("1","1"),
        label     = colnames(X)
    )
    rownames(cd) <- colnames(X)

    rd <- S4Vectors::DataFrame(
        Genes = c("G1","G2"),
        name  = rownames(X),
        ID    = rownames(X)
    )
    rownames(rd) <- rownames(X)

    se_limma <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(intensity = X),
        colData = cd,
        rowData = rd
    )

    # Missing p.adj_<main> should error
    testthat::expect_error(
        ProteinBatcher:::.pp_organize_effects(
            se_limma = se_limma,
            tests = "B_vs_A",
            tests_interaction = "NA",
            alpha = 0.05
        ),
        "Missing rowData column: p\\.adj_"
    )

    # Wrong length for tests_interaction (not 1 and not length(tests))
    rd[["p.adj_B_vs_A"]] <- c(0.1, 0.2)
    SummarizedExperiment::rowData(se_limma) <- rd

    testthat::expect_error(
        ProteinBatcher:::.pp_organize_effects(
            se_limma = se_limma,
            tests = c("B_vs_A", "C_vs_A"),
            tests_interaction = c("X", "Y", "Z"),
            alpha = 0.05
        ),
        "tests_interaction must be length 1 or length\\(tests\\)"
    )
})

