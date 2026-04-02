#' Validate input to limma
#' @keywords internal
#' @noRd
.tl_check_inputs <- function(se, design_formula, ref_condition, test, test_interaction){
    assertthat::assert_that(
        inherits(se, "SummarizedExperiment"),
        is(design_formula, "formula"),
        is.character(ref_condition), length(ref_condition) == 1,
        is.character(test), length(test) >= 1,
        is.character(test_interaction)
    )
    cd <- as.data.frame(SummarizedExperiment::colData(se))
    need_cd <- c("label", "condition", "replicate", "batch", "donor_id")
    if (!all(need_cd %in% colnames(cd))) {
        stop("Missing in colData(se): ", paste(setdiff(need_cd, colnames(cd)), collapse=", "), call. = FALSE)
    }
    rd <- as.data.frame(SummarizedExperiment::rowData(se))
    if (!("name" %in% colnames(rd))) stop("rowData(se) must contain column 'name'.", call. = FALSE)
    imp <- S4Vectors::metadata(se)$imputation_map
    if (is.null(imp)) stop("metadata(se)$imputation_map is missing.", call. = FALSE)
}

#' Aggregation of imputation map by condition
#' @keywords internal
#' @noRd
.tl_condition_imputation_map <- function(se, cd){
    imp_map <- S4Vectors::metadata(se)$imputation_map
    conds <- levels(factor(cd$condition))
    out <- lapply(conds, function(cond){
        sub <- imp_map[, cd$condition == cond, drop = FALSE]
        apply(sub, 1, function(x){
            if (any(x == "ldv")) "ldv" else if (any(x == "mean")) "mean" else "none"
        })
    })
    names(out) <- conds
    out
}

#' Normalize the formuala for model matrix
#' @keywords internal
#' @noRd
.tl_normalize_formula_for_model_matrix <- function(f) {
    f_chr <- as.character(f)

    # case update(): . ~ RHS
    if (length(f_chr) == 3 && f_chr[2] == ".") {
        return(stats::as.formula(paste("~", f_chr[3])))
    }

    # case normal: ~ RHS
    if (length(f_chr) == 2 && f_chr[1] == "~") {
        return(f)
    }

    stop("Unsupported design_formula structure: ", deparse(f), call. = FALSE)
}

#' Robust design matrix
#' @keywords internal
#' @noRd
.tl_make_design <- function(design_formula, cd){
    design_formula <- .tl_normalize_formula_for_model_matrix(design_formula)
    design <- stats::model.matrix(design_formula, data = cd)
    colnames(design) <- make.names(gsub("^condition", "", colnames(design)))
    design
}

#' Fit limma with or without "block" (duplicateCorrelation)
#' @keywords internal
#' @noRd
.tl_fit_limma <- function(raw, design, cd, block_effect){
    if (!block_effect) {
        message("Fitting limma model without block")
        return(limma::lmFit(raw, design = design))
    }
    if (!("block" %in% colnames(cd))) stop("Block variable is missing in colData(se).", call. = FALSE)
    if(length(unique(cd$block)) <= 1) stop("Block only applicable with 2 or more factor levels")
    message("Fitting limma model with block")
    corfit <- limma::duplicateCorrelation(raw, design, block = cd$block)
    limma::lmFit(raw, design, block = cd$block, correlation = corfit$consensus)
}

#' Execute contrasts + BH "recomputed" after filtering genes ldv/ldv
#' @keywords internal
#' @noRd
.tl_run_contrasts <- function(fit0, design, tests, cond_imp_map, map_mn, ref_condition){
    cn_expr <- gsub("_vs_", " - ", tests)
    cn_name <- gsub(" - ", "_vs_", cn_expr)

    cm <- limma::makeContrasts(contrasts = stats::setNames(cn_expr, cn_name),
                                levels = design)
    fit <- limma::eBayes(limma::contrasts.fit(fit0, cm))

    cond_cols <- intersect(colnames(design), names(map_mn))

    one <- function(comp_name, comp_expr){
        tt <- limma::topTable(fit, coef = comp_name, number = Inf, sort.by = "none",
                              adjust.method = "none", confint = TRUE)
        gene <- rownames(tt)

        present <- cond_cols[vapply(cond_cols,function(cc) grepl(paste0("\\b", cc, "\\b"), comp_expr), logical(1))]
        c1 <- if (length(present) >= 1) map_mn[[present[1]]] else ref_condition
        c2 <- if (length(present) >= 2) map_mn[[present[2]]] else ref_condition

        m1 <- cond_imp_map[[c1]][gene]; m2 <- cond_imp_map[[c2]][gene]
        keep <- !(m1 == "ldv" & m2 == "ldv")

        adj <- rep(NA_real_, length(gene))
        adj[keep] <- stats::p.adjust(tt$P.Value[keep], method = "BH")
        tt$adj.P.Val <- adj

        dplyr::tibble(
            ID = gene, comparison = comp_name,
            log2FC = tt$logFC, CI.L = tt$CI.L, CI.R = tt$CI.R,
            p.val = tt$P.Value, p.adj = tt$adj.P.Val
        )
    }

    long <- dplyr::bind_rows(Map(one, cn_name, cn_expr))
    tidyr::pivot_wider(long, id_cols = "ID",
                       names_from = "comparison",
                       values_from = c("log2FC","CI.L","CI.R","p.val","p.adj"),
                       names_sep = "_")
}

#' Clean merge to rowData(se)
#' @keywords internal
#' @noRd
.tl_merge_rowdata <- function(se, wide_tbl){
    rd <- as.data.frame(SummarizedExperiment::rowData(se))
    rd$ID <- rownames(rd)
    out <- dplyr::left_join(rd, wide_tbl, by = "ID")
    rownames(out) <- out$ID
    out$ID <- NULL
    out
}

#' Full contrast of interaction effects and base effects
#' @keywords internal
#' @noRd
.tl_build_full_contrasts <- function(test, test_interaction){
    main <- gsub("_vs_", " - ", test)
    inter <- gsub("_vs_", " - ", test_interaction)
    out <- paste0("(", main, ") + ", inter)
    gsub(" - ", "_vs_", out)
}
