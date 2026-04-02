#' Filter proteins by missingness within conditions (SummarizedExperiment)
#' @description
#' `filter_se_missing()` removes proteins (rows) that fail a per-condition
#' missingness threshold. For each condition in `colData(se)$condition`, the
#' function computes the **proportion of non-missing intensities per protein**
#' and keeps a protein if it meets the threshold in **at least one** condition.
#'
#' This is useful in label-free proteomics pipelines prior to imputation and
#' differential analysis, ensuring that proteins are not uniformly sparse across
#' all conditions.
#'
#' @details
#' Let \eqn{X} be the assay matrix (`assay(se)`), rows = proteins and columns =
#' samples. For every distinct `condition` in `colData(se)$condition`, we
#' compute `rowMeans(!is.na(X[, condition == cn]))`. A protein is **kept** if
#' the per-condition proportion of observed values is `>= percentage / 100` for
#' **any** condition. Otherwise, it is **removed**.
#'
#' @param se A `SummarizedExperiment` (or compatible) object with:
#'   - an assay matrix accessible via `assay(se)` (numeric, may contain `NA`)
#'   - `colData(se)$condition` (character or factor) giving the experimental
#'     condition per sample.
#'   - (optional) `rowData(se)$Protein.Names` and `rowData(se)$Genes` used to
#'     annotate the returned table of removed proteins. If missing, they will be
#'     safely filled with `NA_character_`.
#' @param percentage Numeric scalar in `[0, 100]`. Minimum proportion of
#'   **non-missing** values required within a condition for a protein to be
#'   retained. Default: `50` (i.e., at least half of samples observed within
#'   any one condition).
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{`se_filt`}{`SummarizedExperiment` containing only retained proteins.}
#'   \item{`removed`}{A `data.frame` with metadata for proteins that were
#'      filtered out, columns: `protein_id`, `gene_name`, and `reason`.}
#' }
#'
#' @section Notes:
#' - Samples with `NA` in `colData(se)$condition` are **excluded** from
#'   per-condition computations and a message is emitted.
#' - If `percentage = 0`, all proteins with at least one observed value in any
#'   condition are kept. If `percentage = 100`, only proteins fully observed
#'   within at least one condition are kept.
#'
#' @importFrom SummarizedExperiment assay rowData colData
#' @importFrom methods is

#' @examples
#' library(SummarizedExperiment)
#'
#' ## Toy expression matrix: 4 proteins x 4 samples
#' mat <- matrix(
#'   c(
#'     1,  NA,  2,  3,   # p1: OK in cond B
#'     NA, NA,  NA, NA, # p2: always missing -> removed
#'     5,  6,   NA, NA, # p3: OK in cond A
#'     NA, NA,  1,  NA  # p4: too sparse everywhere
#'   ),
#'   nrow = 4, byrow = TRUE,
#'   dimnames = list(
#'     paste0("p", 1:4),
#'     paste0("s", 1:4)
#'   )
#' )
#'
#' ## Sample conditions
#' cd <- DataFrame(
#'   condition = c("A", "A", "B", "B"),
#'   row.names = colnames(mat)
#' )
#'
#' ## Protein annotations (optional)
#' rd <- DataFrame(
#'   Protein.Names = paste0("Prot", 1:4),
#'   Genes = paste0("Gene", 1:4),
#'   row.names = rownames(mat)
#' )
#'
#' se <- SummarizedExperiment(
#'   assays = list(intensity = mat),
#'   colData = cd,
#'   rowData = rd
#' )
#'
#' ## Keep proteins with >= 50% observed values in at least one condition
#' out <- filter_se_missing(se, percentage = 50)
#'
#' ## Filtered SummarizedExperiment
#' out$se_filt
#'
#' ## Table of removed proteins with reasons
#' out$removed
#'
#' @seealso
#' - `SummarizedExperiment::SummarizedExperiment`
#' - Downstream steps like `impute_se()` for imputation.
#'
#' @export
filter_se_missing <- function(se, percentage = 50) {
    if (!methods::is(se, "SummarizedExperiment")) {
        stop("`se` must be a SummarizedExperiment.")
    }
    if (!is.numeric(percentage) || length(percentage) != 1L ||
        is.na(percentage)) {
        stop("`percentage` must be a single numeric value.")
    }
    if (percentage < 0 || percentage > 100) {
        stop("`percentage` must be in [0, 100].")
    }
    if (!("condition" %in% colnames(SummarizedExperiment::colData(se)))) {
        stop("`colData(se)$condition` is required.")
    }
    min_prop <- percentage / 100
    mat  <- SummarizedExperiment::assay(se)
    meta <- SummarizedExperiment::colData(se)
    if (!is.numeric(mat)) {
        stop("`assay(se)` must be numeric (log-intensities or similar).")
    }
    # Handle missing condition labels
    cond_vec <- as.character(meta$condition)
    if (anyNA(cond_vec)) {
        bad <- which(is.na(cond_vec))
        message(sprintf(
            "Excluding %d sample(s) with NA condition from missingness
            computation.",
            length(bad)
        ))
    }
    valid <- !is.na(cond_vec)
    conds <- unique(cond_vec[valid])
    if (length(conds) == 0L) {
        stop("No valid (non-NA) conditions found in `colData(se)$condition`.")
    }
    keep <- rep(FALSE, nrow(mat))
    for (cn in conds) {
        idx <- which(valid & cond_vec == cn)
        # Guard: if a condition has zero valid samples, skip (shouldn't happen)
        if (length(idx) == 0L) next
        sub <- mat[, idx, drop = FALSE]
        prop_non_na <- rowMeans(!is.na(sub))
        keep <- keep | (prop_non_na >= min_prop)
    }
    se_filt    <- se[keep, ]
    se_removed <- se[!keep, ]
    # Removed table with robust annotations
    rd <- SummarizedExperiment::rowData(se_removed)
    protein_id <- if ("Protein.Names" %in%
                      colnames(rd)) rd$Protein.Names else NA_character_
    gene_name  <- if ("Genes"         %in%
                      colnames(rd)) rd$Genes          else NA_character_
    removed <- data.frame(
        protein_id = as.character(protein_id),
        gene_name  = as.character(gene_name),
        reason     = sprintf(
            "Failed per-condition non-missing threshold: >= %.2f within
            any condition",
            min_prop
        ),
        stringsAsFactors = FALSE
    )
    # Named return for clarity
    out <- list(se_filt = se_filt, removed = removed)
    return(out)
}

#' Impute missing proteomics intensities within conditions
#' (SummarizedExperiment)
#'
#' @description
#' `impute_se()` imputes missing values using a two‑step, condition‑aware
#' strategy: (1) within‑condition mean when missingness < `threshold`; (2)
#' left‑censored draws (LDV) for remaining NAs, using a low‑intensity reference
#' and per‑protein SD.
#'
#' @details
#' Let `assay(se)` be a numeric matrix (rows = proteins, cols = samples);
#' `colData(se)$condition` provides the condition for each sample. Step (1)
#' replaces NAs in a condition with the mean of observed values for that protein
#' in the same condition when the NA proportion is below `threshold`. Step (2)
#' replaces remaining NAs with `rnorm(1, mean = LDV, sd = sd_protein)`, where
#' LDV is the global (or per‑condition) minimum observed intensity; per‑protein
#' SDs are computed from available values with a robust fall-back when needed.
#'
#' The per‑cell imputation method is recorded in `metadata(se)$imputation_map`
#' (`"none"|"mean"|"ldv"`). A row‑level flag is added at `rowData(se)$imputed`.
#'
#' @param se A `SummarizedExperiment` with numeric assay and
#' `colData(se)$condition`.
#' @param threshold Numeric in `(0,1]`, default `0.5`.
#' @param ldv_source `c("global","per_condition")`, default `"global"`.
#' @param sd_fallback Optional numeric fallback SD; default is median non‑zero
#' SD.
#' @param lower_bound Optional minimum clip for imputed values (avoid on log
#' scale).
#'
#' @return A `SummarizedExperiment` with imputed assay, updated
#' metadata/rowData.
#'
#' @examples
#' library(SummarizedExperiment)
#'
#' ## Toy matrix with missing values
#' mat <- matrix(
#'   c(
#'     10, 11,  NA, NA,  # p1: partial missing
#'     NA, NA,  NA, NA, # p2: fully missing
#'     5,  NA,  6,  7   # p3: mixed
#'   ),
#'   nrow = 3, byrow = TRUE,
#'   dimnames = list(
#'     paste0("p", 1:3),
#'     paste0("s", 1:4)
#'   )
#' )
#'
#' ## Conditions for samples
#' cd <- DataFrame(
#'   condition = c("A", "A", "B", "B"),
#'   row.names = colnames(mat)
#' )
#'
#' ## Build SummarizedExperiment
#' se <- SummarizedExperiment(
#'   assays = list(intensity = mat),
#'   colData = cd
#' )
#'
#' ## Impute with default strategy:
#' ## - mean within condition if missingness < threshold
#' ## - LDV otherwise
#' set.seed(123)
#' se_imp <- impute_se(se, threshold = 0.5, ldv_source = "global")
#'
#' ## Imputed assay
#' assay(se_imp)
#'
#' ## Per-cell imputation method
#' metadata(se_imp)$imputation_map
#'
#' ## Row-level flag: was any value imputed?
#' rowData(se_imp)$imputed
#'
#' @seealso SummarizedExperiment::SummarizedExperiment
#' @importFrom SummarizedExperiment assay rowData colData
#' @importFrom S4Vectors metadata
#' @importFrom methods is
#' @importFrom stats rnorm sd
#' @export
impute_se <- function(se, threshold = 0.5,
                      ldv_source = c("global", "per_condition"),
                      sd_fallback = NULL, lower_bound = NULL) {
    args <- .impute_validate_inputs(se, threshold, ldv_source)
    se <- args$se; threshold <- args$threshold; ldv_source <- args$ldv_source
    X0 <- SummarizedExperiment::assay(se)
    cd <- SummarizedExperiment::colData(se)
    excl <- .impute_exclude_na_condition(X0, cd)
    X <- excl$X; cond <- excl$cond; dropped <- excl$dropped
    imp_map <- matrix("none", nrow = nrow(X), ncol = ncol(X),
                      dimnames = list(rownames(X), colnames(X)))
    step1 <- .impute_mean_within_condition(X, cond, threshold, imp_map)
    X1 <- step1$X; imp_map <- step1$map
    ldv <- .impute_ldv_mu(X1, cond, ldv_source)
    sds <- .impute_sd_protein(X, sd_fallback)
    step2 <- .impute_ldv_fill(X1, cond, ldv, sds$sd, ldv_source, lower_bound,
                              imp_map)
    out <- .impute_write_back(se, X0, step2$X, step2$map, dropped)
    return(out)
}
