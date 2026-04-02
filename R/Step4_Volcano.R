#' Volcano plot for differential abundance results
#'
#' Generates a volcano plot from differential testing results stored in a
#' SummarizedExperiment. The function expects contrast-specific
#' columns in \code{rowData(dep)} following the naming convention
#' \code{"<contrast>_diff"} and either \code{"<contrast>_p.adj"} or
#' \code{"<contrast>_p.val"}.
#'
#' Significance is determined internally using a log2 fold-change threshold
#' (\code{lfc}) and a p-value or adjusted p-value threshold (\code{alpha}).
#' Optionally, significant features can be labeled using a column from
#' \code{rowData(dep)} (e.g. gene symbols).
#'
#' @param dep A SummarizedExperiment containing differential
#'   testing results in \code{rowData()}.
#' @param contrast Character scalar specifying the contrast to plot. Must match
#'   the prefix used in \code{rowData(dep)} (e.g. \code{"IL13_vs_NoTreated"}).
#' @param label_size Numeric. Text size for feature labels in the plot.
#' @param name_col Character. Column in \code{rowData(dep)} used for labeling
#'   significant points (default: \code{"ID"}).
#' @param name_imputed Character or \code{NULL}. Column in \code{rowData(dep)}
#'   indicating imputed features (logical). Used to color points. If \code{NULL}
#'   or not present, all features are treated as non-imputed.
#' @param add_names Logical. If \code{TRUE}, labels significant features using
#'   \pkg{ggrepel}.
#' @param adjusted Logical. If \code{TRUE}, uses adjusted p-values
#'   (\code{"<contrast>_p.adj"}); otherwise uses raw p-values
#'   (\code{"<contrast>_p.val"}).
#' @param plot Logical. If \code{TRUE}, returns a \link[ggplot2]{ggplot} object;
#'   otherwise returns a data.frame with volcano statistics.
#' @param alpha Numeric. Significance threshold for p-values.
#' @param lfc Numeric. Absolute log2 fold-change threshold.
#'
#' @details
#' Points are colored according to significance and imputation status:
#' imputed & significant, non-imputed & significant, or non-significant.
#' The plot title and axis annotations are derived from the contrast name.
#'
#' This function does not rely on precomputed \code{*_significant} columns;
#' significance is evaluated on the fly from \code{diff}, \code{p.adj/p.val},
#' \code{alpha}, and \code{lfc}.
#'
#' @return
#' If \code{plot = TRUE}, a \link[ggplot2]{ggplot} volcano plot.
#' If \code{plot = FALSE}, a data.frame with log2 fold change,
#' -log10(p-value), and significance flag.
#'
#' @export
#' @examples
#' library(SummarizedExperiment)
#' library(S4Vectors)
#'
#' ## Minimal SummarizedExperiment with differential results in rowData
#' mat <- matrix(
#'   rnorm(12),
#'   nrow = 3,
#'   dimnames = list(
#'     paste0("p", 1:3),
#'     paste0("s", 1:4)
#'   )
#' )
#'
#' cd <- DataFrame(
#'   condition = c("A","A","B","B"),
#'   row.names = colnames(mat)
#' )
#'
#' rd <- DataFrame(
#'   ID = c("Prot1","Prot2","Prot3"),
#'   imputed = c(FALSE, TRUE, FALSE),
#'   row.names = rownames(mat)
#' )
#'
#' contrast <- "B_vs_A"
#'
#' ## Columns expected by plot_volcano_customized()
#' rd[[paste0("log2FC_", contrast)]] <- c(2.0, 0.3, -1.5)
#' rd[[paste0("p.adj_",  contrast)]] <- c(0.01, 0.9, 0.03)
#'
#' se <- SummarizedExperiment(
#'   assays  = list(intensity = mat),
#'   colData = cd,
#'   rowData = rd
#' )
#'
#' ## Return volcano statistics as a data.frame (no plotting)
#' df <- plot_volcano_customized(
#'   dep = se,
#'   contrast = contrast,
#'   plot = FALSE
#' )
#'
#' df
plot_volcano_customized <- function(
        dep, contrast, label_size = 3, name_col = "ID",
        name_imputed = "imputed", add_names = TRUE, adjusted = TRUE,
        plot = TRUE, alpha = 0.05, lfc = 1
){
    if (is.integer(label_size)) label_size <- as.numeric(label_size)
    assertthat::assert_that(
        inherits(dep, "SummarizedExperiment"),
        is.character(contrast), length(contrast) == 1,
        is.numeric(label_size), length(label_size) == 1,
        is.logical(add_names), length(add_names) == 1,
        is.logical(adjusted), length(adjusted) == 1,
        is.logical(plot), length(plot) == 1
    )
    df <- .pv_build_volcano_df(dep, contrast, name_col, name_imputed, adjusted,
                                alpha, lfc)
    if (!plot) return(.pv_export_df(df, adjusted))
    .pv_plot_volcano(df, contrast, label_size, add_names, adjusted)
}


#' Export volcano plots and result tables from organized limma effects
#'
#' Generates volcano plots (PDF) and corresponding result tables for a given
#' effect type (main, common, or interaction) produced by
#' \code{.pp_organize_effects()}.
#'
#' For each contrast, the function extracts the requested effect slot,
#' generates a volcano plot using \code{\link{plot_volcano_customized}},
#' and returns the corresponding \code{rowData()} as a named list.
#'
#' @param effects Named list as returned by \code{.pp_organize_effects()}.
#'   Each element corresponds to a main contrast and contains
#'   SummarizedExperiment objects for different effect types.
#' @param tests Character vector of main contrast names (e.g.
#'   \code{"IL13_vs_NoTreated"}).
#' @param path_output Character. Output directory where PDF files are written.
#' @param experiment Character. Experiment identifier used in output filenames.
#' @param effect_slot Character. Which effect to export. One of
#'   \code{"all_common_effect"}, \code{"common_effect"}, or
#'   \code{"interaction_effect"}.
#' @param file_tag Character. Tag used in output filenames describing the
#'   exported effect (e.g. \code{"MainEffect_AllProteins"},
#'   \code{"CommonEffect"}, \code{"InteractionEffect"}).
#' @param plot_contrasts Character vector of contrasts to be plotted.
#'   By default, this is set to \code{tests}. For
#'   \code{effect_slot = "interaction_effect"}, this should typically be
#'   the corresponding interaction contrasts (e.g. \code{tests_interaction}),
#'   since interaction result objects do not contain statistics for the
#'   main contrasts.
#' @param name_col Character. Column in \code{rowData()} used for labeling
#'   significant points in volcano plots (default: \code{"Genes"}).
#' @param name_imputed Character. Column in \code{rowData()} indicating
#'   imputed features.
#' @param alpha Numeric. Significance threshold used in volcano plots.
#' @param lfc Numeric. Absolute log2 fold-change threshold used in volcano
#' plots.
#'
#' @details
#' If the selected effect object is \code{NULL} or contains zero rows,
#' the corresponding contrast is skipped silently.
#'
#' Output PDF filenames follow the pattern:
#' \preformatted{
#' <experiment>_<file_tag>_<plot_contrast>.pdf
#' }
#'
#' This function does not perform additional filtering beyond what is already
#' encoded in the provided \code{effects} object.
#'
#' @return
#' A list (one element per main contrast) containing named lists of
#' \code{rowData()} tables. Contrasts that are skipped return \code{NULL}.
#'
#' @seealso
#' \code{\link{.pp_organize_effects}},
#' \code{\link{plot_volcano_customized}}
#'
#' @importFrom stats setNames
#' @keywords internal
.pp_export_volcano_tables <- function(
        effects, tests, path_output, experiment,
        effect_slot = c("all_common_effect","common_effect",
                        "interaction_effect"),
        file_tag = c("MainEffect_AllProteins","CommonEffect",
                    "InteractionEffect"),
        plot_contrasts = NULL,
        name_col = "Genes", name_imputed = "imputed", alpha = 0.05, lfc = 1
){
    effect_slot <- match.arg(effect_slot)
    file_tag <- match.arg(file_tag)
    if (is.null(plot_contrasts)) plot_contrasts <- tests
    if (length(plot_contrasts) == 1L) plot_contrasts <- rep(plot_contrasts,
                                                            length(tests))
    if (length(plot_contrasts) != length(tests))
        stop("plot_contrasts must be length 1 or length(tests).", call.=FALSE)
    lapply(seq_along(tests), function(i){
        obj <- effects[[tests[i]]][[effect_slot]]
        if (is.null(obj) || nrow(obj) == 0) return(NULL)

        pdf_file <- file.path(path_output, paste0(experiment, "_", file_tag,
                                                "_", plot_contrasts[i], ".pdf"))
        grDevices::pdf(pdf_file)
        print(plot_volcano_customized(obj, plot_contrasts[i],
                                        name_col = name_col,
                                        name_imputed = name_imputed,
                                        alpha = alpha, lfc = lfc,
                                        adjusted = TRUE, plot = TRUE))
        grDevices::dev.off()
        rd <- SummarizedExperiment::rowData(obj)
        setNames(list(as.data.frame(rd)), tests[i])
    })
}
