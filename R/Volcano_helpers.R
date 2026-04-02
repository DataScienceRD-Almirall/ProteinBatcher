# Build volcano data frame
#' @keywords internal
#' @noRd
.pv_build_volcano_df <- function(dep, contrast, name_col, name_imputed,
                                    adjusted, alpha, lfc){
    rd <- SummarizedExperiment::rowData(dep)
    rn <- colnames(rd)

    # Accept both naming conventions:
    # 1) <contrast>_log2FC   (old)
    # 2) diff_<contrast>   (new)
    diff_candidates <- c(paste0(contrast, "_log2FC"), paste0("log2FC_", contrast))
    padj_candidates <- c(paste0(contrast, if (adjusted) "_p.adj" else "_p.val"),
                         paste0(if (adjusted) "p.adj_" else "p.val_", contrast))

    diff_col <- diff_candidates[diff_candidates %in% rn][1]
    pcol     <- padj_candidates[padj_candidates %in% rn][1]

    if (is.na(diff_col) || is.na(pcol)) {
        # build "valid contrasts" from whatever exists
        valid <- unique(c(
            sub("^log2FC_", "", rn[startsWith(rn, "_log2FC")]),
            sub("_log2FC$", "", rn[endsWith(rn, "_log2FC")])
        ))
        stop("Not a valid contrast: ", contrast,
             "\nValid contrasts: ", paste(valid, collapse = ", "),
             call. = FALSE)
    }

    if (!name_col %in% rn)
        stop("Column '", name_col, "' not found in rowData(dep).", call.=FALSE)

    imputed <- if (!is.null(name_imputed) && name_imputed %in% rn)
        rd[[name_imputed]] else FALSE
    pv <- rd[[pcol]]
    fc <- rd[[diff_col]]
    keep <- !is.na(pv) & !is.na(fc)
    pv <- pv[keep]
    fc <- fc[keep]

    data.frame(
        log2FC = fc,
        p_values = -log10(pv),
        signif = (abs(fc) >= lfc) & (pv <= alpha),
        label = rd[[name_col]][keep],
        imputed = imputed[keep],
        stringsAsFactors = FALSE
    )
}

# Plot volcano
#' @keywords internal
#' @noRd
.pv_plot_volcano <- function(df, contrast, label_size, add_names, adjusted){
    name1 <- sub("_vs_.*", "", contrast); name2 <- sub(".*_vs_", "", contrast)
    df$color <-
        ifelse(df$signif,ifelse(df$imputed, "imputed_significant",
                                "non_imputed_significant"), "non_significant")

    p <- ggplot2::ggplot(df, ggplot2::aes(log2FC, p_values)) +
        ggplot2::geom_vline(xintercept = 0) +
        ggplot2::geom_point(ggplot2::aes(col = color)) +
        ggplot2::geom_text(data = data.frame(), ggplot2::aes(
            x = c(Inf, -Inf), y = c(-Inf, -Inf), hjust = c(1,0),
            vjust = c(-1,-1), label = c(name1, name2)
        ), size = 5, fontface = "bold") +
        ggplot2::labs(title = contrast,
                        x = expression(log[2] ~ "Fold change")) +
        ggplot2::theme_bw() +
        ggplot2::theme(panel.border = ggplot2::element_blank(),
                       panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(),
                       axis.line = ggplot2::element_line(colour = "black"),
                       legend.position = "bottom") +
        ggplot2::scale_color_manual(
            values = c(imputed_significant="black",
                        non_imputed_significant="red",
                        non_significant="lightgrey"),
            name = "Imputation & Significance",
            breaks = c("imputed_significant","non_imputed_significant",
                        "non_significant"),
            labels = c("Imputed & Significant","Non-Imputed & Significant",
                        "Non-Significant")
        ) +
        ggplot2::labs(y = if (adjusted) expression(-log[10]~ "Adjusted p-value")
                        else expression(-log[10] ~ "P-value"))

    if (add_names) {
        p <- p + ggrepel::geom_text_repel(
            data = df[df$signif, , drop = FALSE],
            ggplot2::aes(label = label),
            size = label_size,
            box.padding = grid::unit(0.1, "lines"),
            point.padding = grid::unit(0.1, "lines"),
            segment.size = 0.5
        )
    }
    p
}

# Export volcano data frame
#' @keywords internal
#' @noRd
.pv_export_df <- function(df, adjusted){
    out <- df[, c("log2FC","p_values","signif"), drop = FALSE]
    out <-
        stats::setNames(out, c("log2_fold_change",
                        if (adjusted) "adjusted_p_value_-log10"
                        else "p_value_-log10", "significant"))
    out
}
