#' Deregulogram plot for interaction models (2-level factors)
#'
#' Creates a deregulogram: a scatter of full effects in two levels of an
#' interaction factor (e.g. sex female vs male, batch day1 vs day2), colored by
#' whether the interaction is significant. This visualization complements
#' volcano plots by showing effect concordance and sex-/batch-bias in magnitude
#' and direction, as described in the proteomics workflow documentation.
#'
#' The function assumes that differential testing has already been performed and
#' stored in \code{rowData(se_limma)} using the \code{"<stat>_<contrast>"}
#' naming convention (e.g., \code{"log2FC_IL13_vs_NoTreated"},
#' \code{"p.adj_IL13_vs_NoTreated"}, \code{"log2FC_IL13.batchday2"},
#' \code{"p.adj_IL13.batchday2"}).
#'
#' Deregulogram axes:
#' \itemize{
#'   \item y-axis: full effect in the reference level (level 1)
#'   \item x-axis: full effect in the non-reference level (level 2), computed as
#'   \code{log2FC_main + log2FC_interaction}
#' }
#'
#' This plot is meaningful only when the interaction factor has exactly two
#' levels. Use \code{.pp_infer_interaction_factor()} (internal) to infer the
#' factor from \code{colData()} and \code{tests_interaction}, or pass factor
#' name and levels explicitly.
#'
#' @details
#' \strong{Important: gene sets do not have to match the interaction volcano
#' plot.}
#' The interaction volcano plot (see \code{plot_volcano_customized} on the
#' \code{interaction_effect} object) identifies proteins with a statistically
#' significant interaction term (e.g. \code{p.adj_interaction < alpha} and
#' \code{|log2FC_interaction| >= lfc}). In contrast, the deregulogram is
#' designed to highlight a more interpretable subset where the interaction is
#' not only statistically significant, but also occurs in the context of a
#' relevant condition effect (i.e., the main/full effects reach the chosen
#' \code{alpha} and \code{lfc} thresholds). Therefore, the deregulogram
#' typically shows a subset of the interaction-volcano hits and is not expected
#' to reproduce the same labeled proteins one-to-one.
#'
#' Deregulogram axes compare full effects between two levels of the interaction
#' factor. For a main contrast \code{C} and a 2-level factor (ref vs level2),
#' the full effect in the reference level is \code{log2FC_C}, while the full
#' effect in the second level is \code{log2FC_C + log2FC_{C:factor}}. Points are
#' colored by whether the interaction is significant and optionally labeled to
#' emphasize proteins whose condition effect is factor-dependent.
#'
#' @param se_limma A SummarizedExperiment containing limma results
#'   in \code{rowData()} (e.g., output of \code{test_limma_customized()}).
#' @param tests Character vector of main contrasts (e.g.
#' \code{"IL13_vs_NoTreated"}).
#' @param tests_interaction Character vector of interaction
#' contrasts/coefficients aligned with \code{tests} (e.g. \code{"IL13.sexmale"}
#' or \code{"IL13.batchday2"}).
#' @param factor_name Character scalar, name of the interaction factor
#'   (e.g. \code{"sex"}, \code{"batch"}). Used only for titles/labels.
#' @param factor_levels Character vector of length 2 giving the factor levels
#' in order \code{c(ref, other)} (e.g. \code{c("female","male")}).
#' @param path_output Output directory where PDFs are written.
#' @param experiment Experiment identifier used in output filenames.
#' @param alpha Numeric. Significance threshold (default 0.05).
#' @param lfc Numeric. Absolute log2FC threshold used for highlighting
#' (default 1).
#' @param label_col Column in \code{rowData(se_limma)} used for point labels
#'   (default \code{"Genes"}).
#'
#' @return Invisibly returns \code{NULL}. Side effect: writes one PDF per
#' contrast.
#'
#' @export
#' @examples
#' library(SummarizedExperiment)
#' library(S4Vectors)
#'
#' ## ---------------------------------------------------------------
#' ## Minimal SummarizedExperiment with main + interaction statistics
#' ## ---------------------------------------------------------------
#'
#' mat <- matrix(
#'   rnorm(12),
#'   nrow = 3,
#'   dimnames = list(
#'     paste0("p", 1:3),
#'     paste0("s", 1:4)
#'   )
#' )
#'
#' ## Sample metadata (interaction factor has 2 levels)
#' cd <- DataFrame(
#'   condition = c("A","A","B","B"),
#'   sex       = c("female","male","female","male"),
#'   row.names = colnames(mat)
#' )
#'
#' ## Main and interaction contrasts
#' main_contrast  <- "B_vs_A"
#' inter_contrast <- "B_vs_A.sexmale"
#'
#' ## rowData with columns expected by plot_deregulogram()
#' rd <- DataFrame(
#'   Genes = c("G1","G2","G3"),
#'   row.names = rownames(mat)
#' )
#'
#' ## Main effect
#' rd[[paste0("log2FC_", main_contrast)]] <- c( 1.5,  0.2, -1.3)
#' rd[[paste0("p.adj_",  main_contrast)]] <- c( 0.01, 0.8,  0.03)
#'
#' ## Interaction effect
#' rd[[paste0("log2FC_", inter_contrast)]] <- c( 0.7, -0.1,  1.1)
#' rd[[paste0("p.adj_",  inter_contrast)]] <- c( 0.02, 0.9,  0.01)
#'
#' se_limma <- SummarizedExperiment(
#'   assays  = list(intensity = mat),
#'   colData = cd,
#'   rowData = rd
#' )
#'
#' ## ---------------------------------------------------------------
#' ## Create deregulogram (PDF written to tempdir())
#' ## ---------------------------------------------------------------
#'
#' plot_deregulogram(
#'   se_limma = se_limma,
#'   tests = main_contrast,
#'   tests_interaction = inter_contrast,
#'   factor_name = "sex",
#'   factor_levels = c("female", "male"),
#'   path_output = tempdir(),
#'   experiment = "EXAMPLE",
#'   alpha = 0.05,
#'   lfc = 1,
#'   label_col = "Genes"
#' )
plot_deregulogram <- function(
        se_limma, tests, tests_interaction, factor_name, factor_levels,
        path_output, experiment, alpha = 0.05, lfc = 1, label_col = "Genes"
){
    stopifnot(length(factor_levels) == 2)
    rd <- SummarizedExperiment::rowData(se_limma)
    rn <- colnames(rd)
    for (i in seq_along(tests)) {
        main <- tests[i]
        inter <- tests_interaction[i]
        if ("NA" %in% inter) next
        cols_needed <- .pd_cols_needed(main, inter, label_col)
        if (!all(cols_needed %in% rn)) next
        df <- as.data.frame(rd[, cols_needed])
        df <- df[stats::complete.cases(df), , drop = FALSE]
        if (nrow(df) == 0) next
        df$fc_ref  <- df[[paste0("log2FC_", main)]]
        df$fc_lvl2 <- df$fc_ref + df[[paste0("log2FC_", inter)]]
        df$p_ref   <- df[[paste0("p.adj_", main)]]
        df$p_lvl2  <- df$p_ref
        df$p_inter <- df[[paste0("p.adj_", inter)]]
        df$fc_inter <- df[[paste0("log2FC_", inter)]]
        df$class <- .pd_classify(df$p_ref, df$p_lvl2, df$fc_ref, df$fc_lvl2,
                                 df$p_inter, df$fc_inter, alpha, lfc)
        out_file <- file.path(path_output,
                              paste0(experiment, "_Deregulogram_", main, "_",
                                    factor_name, ".pdf"))
        grDevices::pdf(out_file, width = 14, height = 10)
        print(.pd_plot(df, main, factor_name, factor_levels, label_col))
        grDevices::dev.off()
    }
    invisible(NULL)
}
