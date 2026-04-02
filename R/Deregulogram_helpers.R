#' @keywords internal
.pd_cols_needed <- function(main, inter, label_col){
    c(paste0("log2FC_", main),
      paste0("p.adj_", main),
      paste0("log2FC_", inter),
      paste0("p.adj_", inter),
      label_col)
}

#' @keywords internal
.pd_classify <- function(p_ref, p_lvl2, fc_ref, fc_lvl2, p_inter, fc_inter,
                        alpha, lfc){
    ifelse(
        p_ref <= alpha & p_lvl2 <= alpha &
            abs(fc_ref) >= lfc & abs(fc_lvl2) >= lfc &
            p_inter <= alpha & abs(fc_inter) >= lfc,
        "sig_interaction",
        ifelse(
            p_ref <= alpha & p_lvl2 <= alpha &
                (abs(fc_ref) >= lfc | abs(fc_lvl2) >= lfc),
            "sig_no_interaction",
            "non_sig"
        )
    )
}

#' @keywords internal
.pd_plot <- function(df, main, factor_name, factor_levels, label_col){
    df_subset <- df[df$class == "sig_interaction", , drop = FALSE]
    ggplot2::ggplot(df, ggplot2::aes(x = fc_lvl2, y = fc_ref, color = class)) +
        ggplot2::geom_point() +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
        ggrepel::geom_label_repel(
            data = df_subset,
            ggplot2::aes(label = .data[[label_col]]),
            na.rm = TRUE,
            max.overlaps = Inf
        ) +
        ggplot2::scale_color_manual(
            values = c(sig_interaction = "red",
                       sig_no_interaction = "blue",
                       non_sig = "grey70"),
            name = "Significance",
            labels = c(sig_interaction = paste0("Significant condition x",
                                                factor_name),
                       sig_no_interaction = "Significant condition only",
                       non_sig = "Not significant")
        ) +
        ggplot2::labs(
            title = paste0(main, " - Deregulogram (", factor_name, ": ",
                           factor_levels[1], " vs ", factor_levels[2], ")"),
            x = paste0("Full effect (", factor_levels[2], ")"),
            y = paste0("Full effect (", factor_levels[1], ")")
        ) +
        ggplot2::theme_minimal()
}

#' Infer a 2-level interaction factor from colData() and tests_interaction
#'
#' Internal helper to infer the interaction factor name and its two levels by
#' scanning \code{colData(se)} for variables with exactly two levels and
#' matching the non-reference level against \code{tests_interaction} (e.g.
#' "sexmale", "batchday2"). Returns \code{NULL} if no suitable factor is found.
#'
#' @param se SummarizedExperiment.
#' @param tests_interaction Character vector of interaction
#' contrasts/coefficients.
#' @return A list with \code{name} and \code{levels} (length 2) or \code{NULL}.
#' @keywords internal
.pp_infer_interaction_factor <- function(se, tests_interaction){
    cd <- as.data.frame(SummarizedExperiment::colData(se),
                        stringsAsFactors = FALSE)
    cand <- setdiff(names(cd), c("condition","label","file","sample",
                                "replicate"))
    ti <- tests_interaction[!("NA" %in% tests_interaction)]
    if (!length(ti)) return(NULL)

    for (v in cand) {
        lev <- levels(factor(cd[[v]]))
        if (length(lev) != 2) next
        key1 <- paste0(v, make.names(lev[2]))
        key2 <- make.names(lev[2])
        if (any(grepl(key1, ti, fixed = TRUE)) || any(grepl(key2, ti,
                                                            fixed = TRUE))) {
            return(list(name = v, levels = lev))
        }
    }
    NULL
}
