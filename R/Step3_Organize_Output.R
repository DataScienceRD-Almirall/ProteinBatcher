#' Organize limma outputs into common vs interaction effects (plot-ready)
#'
#' @param se_limma SummarizedExperiment returned by test_limma_customized().
#' @param tests Character vector of main contrasts (prefixes in rowData), e.g.
#' "IL13_vs_NoTreated".
#' @param tests_interaction Character vector of interaction tests (prefixes in
#' rowData) or "NA". If length 1, it will be recycled to match length(tests).
#' @param alpha FDR threshold for interaction significance (default 0.05).
#' @param base_cols Annotation columns to keep (if present).
#' @return Named list (by `tests`). Each element is a list with:
#'   common_effect, interaction_effect, all_common_effect
#'   (SummarizedExperiment objects).
#' @export
#' @examples
#' if (requireNamespace("limma", quietly = TRUE)) {
#'   library(SummarizedExperiment)
#'   library(S4Vectors)
#'
#'   ## Minimal imputed SE-like object
#'   mat <- matrix(
#'     c(10,11, 9,10,
#'       20,21,19,20,
#'       5, NA, 6, 7),
#'     nrow = 3, byrow = TRUE,
#'     dimnames = list(paste0("p",1:3), paste0("s",1:4))
#'   )
#'
#'   cd <- DataFrame(
#'     condition = c("A","A","B","B"),
#'     batch     = c("d1","d2","d1","d2"),
#'     replicate = c("r1","r2","r1","r2"),
#'     donor_id  = c("1","1","1","1"),
#'     label     = colnames(mat),
#'     row.names = colnames(mat)
#'   )
#'
#'   rd <- DataFrame(
#'     name  = rownames(mat),
#'     ID    = rownames(mat),
#'     Genes = paste0("G",1:3),
#'     row.names = rownames(mat)
#'   )
#'
#'   se_imp <- SummarizedExperiment(
#'     assays  = list(intensity = mat),
#'     colData = cd,
#'     rowData = rd
#'   )
#'
#'   ## minimal imputation map required by your limma wrapper
#'   metadata(se_imp)$imputation_map <- matrix(
#'     "none", nrow = nrow(mat), ncol = ncol(mat),
#'     dimnames = dimnames(mat)
#'   )
#'
#'   ## Run limma wrapper -> this is se_limma
#'   se_limma <- test_limma_customized(
#'     se = se_imp,
#'     test = "B_vs_A",
#'     test_interaction = "NA",
#'     design_formula = ~ 0 + condition + batch,
#'     ref_condition = "A"
#'   )
#'
#'   ## Organize effects (no interaction)
#'   eff <- ProteinBatcher:::.pp_organize_effects(
#'     se_limma = se_limma,
#'     tests = "B_vs_A",
#'     tests_interaction = "NA"
#'   )
#'
#'   eff[["B_vs_A"]]$all_common_effect
#' }
#' @keywords internal
.pp_organize_effects <- function(
        se_limma, tests, tests_interaction = "NA", alpha = 0.05,
        base_cols = c("Protein.Group","Protein.Names","Genes",
                      "First.Protein.Description","name","ID","Index","imputed")
){
    rd <- SummarizedExperiment::rowData(se_limma)
    rn <- colnames(rd)
    if (length(tests_interaction) == 1L)
        tests_interaction <- rep(tests_interaction, length(tests))
    if (length(tests_interaction) != length(tests))
        stop("tests_interaction must be length 1 or length(tests).",call.=FALSE)

    keep_for <- function(pref){
        base <- intersect(base_cols, rn)
        # new style: <stat>_<pref>
        stats <- c("diff","p.val","p.adj","CI.L","CI.R")
        new_style <- paste0(stats, "_", pref)
        new_style <- new_style[new_style %in% rn]
        # old style: <pref>_<stat>
        old_style <- paste0(pref, "_", stats)
        old_style <- old_style[old_style %in% rn]
        # fallback (covers any future extras): either starts with
        # "<pref>_" or ends with "_<pref>"
        loose <- rn[startsWith(rn, paste0(pref, "_")) |
                        endsWith(rn,paste0("_", pref))]
        unique(c(base, new_style, old_style, loose))
    }

    slice_se <- function(idx, cols){
        x <- se_limma[idx, , drop = FALSE]
        SummarizedExperiment::rowData(x) <-
            SummarizedExperiment::rowData(x)[, cols, drop = FALSE]
        x
    }

    out <- vector("list", length(tests))
    for (i in seq_along(tests)) {
        main <- tests[i]
        main_p <- paste0("p.adj_", main)
        if (!main_p %in% rn)
            stop("Missing rowData column: ", main_p, call.=FALSE)

        cols_main <- keep_for(main)
        all_common <- slice_se(TRUE, cols_main)

        inter <- tests_interaction[i]
        if ("NA" %in% inter) {
            out[[i]] <- list(common_effect = all_common,
                             interaction_effect = slice_se(FALSE, cols_main),
                             all_common_effect = all_common)
            next
        }
        inter_p <- paste0("p.adj_", inter)
        if (!inter_p %in% rn)
            stop("Missing rowData column: ", inter_p, call.=FALSE)
        sig_inter <- rd[[inter_p]] <= alpha
        sig_inter[is.na(sig_inter)] <- FALSE
        common <- slice_se(!sig_inter, cols_main)
        cols_int <- keep_for(inter)
        inter_se <- slice_se(sig_inter, cols_int)
        out[[i]] <- list(common_effect = common,
                         interaction_effect = inter_se,
                         all_common_effect = all_common)
    }
    names(out) <- tests
    out
}
