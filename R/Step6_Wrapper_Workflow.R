#' Run a the downstream proteomics pipeline (filter, impute, limma, reporting)
#'
#' This workflow performs a reproducible downstream analysis for quantitative
#' proteomics: (i) imports quantification and sample annotation into a
#' SummarizedExperiment, (ii) filters features by missingness,
#' (iii) imputes remaining missing values using a mean vs LDV (left-censored)
#' strategy, (iv) runs differential abundance testing with limma (optionally
#' including an interaction term), (v) organizes results into plot-ready
#' effect objects, and (vi) optionally exports volcano plots, deregulograms and
#' results to disk.
#'
#' @section Input files:
#' \describe{
#'   \item{\code{path_pgmatrix}}{
#'     Path to a quantitative matrix file produced by the upstream
#'     quantification tool (e.g., a TSV). Rows represent quantified features
#'     (e.g. protein groups) and columns represent samples.
#'   }
#'   \item{\code{path_annotation}}{
#'     Path to a tab-delimited sample annotation file describing the
#'     experimental design. At minimum it must encode sample-to-condition
#'     mapping. Additional columns are required
#'     (\code{batch}, \code{replicate}, \code{block}).
#'   }
#' }
#'
#' @section Main steps:
#' \enumerate{
#'   \item \strong{Import}: creates a raw SummarizedExperiment
#'   (\code{se_raw}).
#'   \item \strong{Missingness filter}: removes features not sufficiently
#'   observed across samples/conditions according to \code{percent_missing}.
#'   Filtered-out features are returned and may be written to disk.
#'   \item \strong{Imputation}: imputes remaining missing values using a
#'   two-step strategy: within-condition mean imputation when missingness is
#'   below a threshold and LDV (left-censored) imputation when missingness is
#'   above the threshold. The imputation method is controlled by
#'   \code{ldv_source} and optional arguments passed via \code{...}.
#'   \item \strong{Differential testing}: calls
#'   \code{\link{test_limma_customized}} to fit a limma model and evaluate
#'   user-specified contrasts. Interaction testing is optional.
#'   \item \strong{Organize outputs}: splits the limma results into plot-ready
#'   effect objects using \code{.pp_organize_effects}:
#'         \itemize{
#'           \item \code{all_common_effect}: main-effect statistics for all
#'           proteins (no interaction filtering).
#'           \item \code{common_effect}: proteins not significant for the
#'           interaction term (main effect stable).
#'           \item \code{interaction_effect}: proteins significant for the
#'           interaction term.
#'         }
#'   \item \strong{Optional reporting} (\code{plots = TRUE}): exports volcano
#'   plots and tables for the three effect types above, and generates
#'   deregulograms when an interaction factor with exactly two levels is
#'   detected (e.g., female vs male, day1 vs day2, ...).
#' }
#'
#' @section Differential testing inputs (passed via \code{...}):
#' The limma step requires the following objects, typically supplied via
#' \code{...}:
#' \describe{
#'   \item{\code{tests}}{
#'     Character vector of main contrasts to test, using a consistent naming
#'     scheme (e.g. \code{"A_vs_B"}). These define the primary condition
#'     effects.
#'   }
#'   \item{\code{tests_interaction}}{
#'     Character vector aligned with \code{tests} defining interaction
#'     contrasts or interaction coefficients (or \code{"NA"} to disable
#'     interaction testing). These represent factor-dependent differences in the
#'     condition effect (e.g. condition-by-batch).
#'   }
#'   \item{\code{formula}}{
#'     A model formula describing the design matrix (e.g.
#'     \code{~ 0 + condition + batch + condition:batch}). Variables referenced
#'     in the formula must be present in \code{colData(se)}.
#'   }
#'   \item{\code{reference_condition}}{
#'     Character scalar indicating the reference level for \code{condition}.
#'     The function will relevel \code{colData(se)$condition} to this value
#'     before model fitting.
#'   }
#' }
#'
#' @section Paired designs and blocking (passed via \code{...}):
#' \describe{
#'   \item{\code{paired}}{
#'     Logical (forwarded to \code{test_limma_customized}). If \code{TRUE},
#'     the model is treated as a paired / repeated-measures design by ensuring
#'     \code{replicate} is included in the design formula. In the current
#'     implementation, if \code{replicate} is not already present in
#'     \code{formula}, it is added automatically. This preserves any additional
#'     covariates already present in the formula.
#'   }
#'   \item{\code{block_effect}}{
#'     Logical (forwarded to \code{test_limma_customized}). If \code{TRUE},
#'     correlation between repeated observations is modeled using
#'     \code{limma::duplicateCorrelation} and
#'     \code{lmFit(..., block=..., correlation=...)}.
#'     This requires a \code{block} column in \code{colData(se)}. Use this when
#'     there is a known blocking factor (e.g., donor/subject) inducing
#'     correlation across samples.
#'   }
#' }
#'
#' @section Volcano vs deregulogram (interpretation note):
#' Interaction volcano plots highlight proteins with significant interaction
#' coefficients (i.e., where the condition effect depends on a second factor).
#' Deregulograms are more restrictive: they visualize full effects across the
#' two levels of the interaction factor and emphasize proteins where interaction
#' occurs in the context of a relevant condition effect (effect-size and FDR
#' thresholds). Therefore, the labeled proteins in the deregulogram are not
#' expected to match one-to-one with those in the interaction volcano plot.
#'
#' @param path_pgmatrix Character scalar. Path to the quantitative matrix file.
#' @param path_annotation Character scalar. Path to the sample annotation TSV.
#' @param path_output Character scalar. Output directory where files may be
#' written.
#' @param level Character. Quantification level (e.g. \code{"protein"}).
#' @param type Character. Quantification type (e.g. \code{"DIA"}).
#' @param experiment Character scalar. Experiment identifier used in output
#' filenames.
#' @param percent_missing Numeric. Missingness threshold used for filtering
#' (percentage, e.g. 50).
#' @param ldv_source Character. One of \code{"global"} or \code{"per-condition"}
#' controlling LDV definition.
#' @param plots Logical. If \code{TRUE}, exports plots/tables to
#' \code{path_output}.
#'
#' @param tests Character vector of main contrasts to test.
#'   Each element must correspond to a valid condition contrast present
#'   in the design matrix (e.g. \code{"IL13_vs_NoTreated"}).
#'
#' @param tests_interaction Character vector defining interaction contrasts
#'   or interaction coefficients aligned with \code{tests}. Use \code{"NA"}
#'   to disable interaction testing for a given contrast. Typical values
#'   encode condition-by-factor interactions (e.g. \code{"IL13.batchday2"}).
#'
#' @param formula A model formula specifying the design matrix for limma.
#'   Common examples include \code{~ 0 + condition + batch + condition:batch}.
#'   All variables referenced in the formula must be present in
#'   \code{colData()}.
#'
#' @param reference_condition Character scalar specifying the reference
#'   level for the \code{condition} variable. This level is used to relevel
#'   the design prior to fitting the limma model.
#'
#' @param ... Additional arguments forwarded to \code{impute_se()} and
#' \code{test_limma_customized()}.
#'   In particular, \code{tests}, \code{tests_interaction}, \code{formula},
#'   \code{reference_condition}, and optionally \code{paired} /
#'   \code{block_effect} should be provided here.
#'
#' @return A named list with:
#' \describe{
#'   \item{\code{se_raw}}{Raw SummarizedExperiment
#'   imported from files.}
#'   \item{\code{se_filt}}{Filtered SummarizedExperiment.}
#'   \item{\code{removed}}{Object describing filtered-out features (typically a
#'   data.frame).}
#'   \item{\code{se_imp}}{Imputed SummarizedExperiment.}
#'   \item{\code{effects}}{
#'     A list of exported objects (or \code{NULL} entries when
#'     \code{plots = FALSE}):
#'     main-effect exports, interaction-effect exports, and common-effect
#'     exports.
#'   }
#' }
#'
#' @seealso \code{\link{test_limma_customized}},
#' \code{\link{plot_volcano_customized}},
#'   \code{\link{plot_deregulogram}}
#'
#' @examples
#' ## Load package
#' library(ProteinBatcher)
#'
#' ## ---------------------------------------------------------------
#' ## Input files shipped with the package
#' ## ---------------------------------------------------------------
#'
#' path_annotation <- system.file(
#'   "extdata", "annotation_HaCaT.tsv",
#'   package = "ProteinBatcher"
#' )
#'
#' path_pgmatrix <- system.file(
#'   "extdata",
#'   "2024MK017_HaCaT_Stimulation_1to24_Astral_report.pg_matrix.tsv",
#'   package = "ProteinBatcher"
#' )
#'
#' ## Output directory
#' path_output <- tempdir()
#'
#' ## ---------------------------------------------------------------
#' ## Workflow parameters (single contrast example)
#' ## ---------------------------------------------------------------
#'
#' experiment <- "HaCaT"
#' percent_missing <- 50
#'
#' formula <- ~ 0 + condition + batch + condition:batch
#'
#' tests <- c("IL13_vs_NoTreated")
#' tests_interaction <- c("IL13.batchday2")
#'
#' reference_condition <- "NoTreated"
#'
#' ## ---------------------------------------------------------------
#' ## Run downstream proteomics workflow
#' ## ---------------------------------------------------------------
#'
#' res <- run_proteomics_pipeline(
#'   path_pgmatrix        = path_pgmatrix,
#'   path_annotation      = path_annotation,
#'   path_output          = path_output,
#'   tests                = tests,
#'   tests_interaction    = tests_interaction,
#'   formula              = formula,
#'   reference_condition  = reference_condition,
#'   percent_missing      = percent_missing,
#'   ldv_source           = "per-condition",
#'   experiment           = experiment,
#'   plots                = FALSE
#' )
#'
#' ## Inspect outputs
#' names(res)
#' res$se_imp
#' names(res)
#' @export
run_proteomics_pipeline <- function(
        path_pgmatrix, path_annotation, path_output, level = "protein",
        type = "DIA", experiment, percent_missing,
        ldv_source = c("global", "per-condition"),
        tests, tests_interaction, formula, reference_condition,
        plots = FALSE, ...
){
    args <- .pp_validate_inputs(path_pgmatrix, path_annotation, level, type,
                                percent_missing, path_output, experiment)
    se0 <- .pp_import_se(args$path_pgmatrix, args$path_annotation,
                         args$level, args$type)
    # 1) Filter and imputation
    filt <- filter_se_missing(se0, percentage = args$percent_missing)
    .pp_write_filtered(filt$removed, args$path_output, args$experiment)
    se_imp <- impute_se(filt$se_filt, ldv_source = ldv_source, ...)
    .pp_write_before_after(filt$se_filt, se_imp,
                           args$path_output, args$experiment)
    # 2) Differential testing
    se_limma <- test_limma_customized(
        se_imp, type = "manual", test = tests,
        test_interaction = tests_interaction, design_formula = formula,
        ref_condition  = reference_condition, ...
    )
    # 3) Organize outputs (plot-ready)
    effects <- .pp_organize_effects(
        se_limma = se_limma, tests = tests,
        tests_interaction = tests_interaction, alpha = 0.05
    )
    if(plots){
        # 4) Volcano plots
        # Main effects: effects without significant interaction
        output_main_all <- .pp_export_volcano_tables(
            effects = effects, tests = tests, path_output = path_output,
            experiment = experiment, effect_slot = "all_common_effect",
            file_tag = "MainEffect_AllProteins", name_col = "Genes",
            name_imputed = "imputed", ...
        )
        # Common: effects significant for condition and NOT interacting
        output_common <- .pp_export_volcano_tables(
            effects = effects, tests = tests, path_output = path_output,
            experiment = experiment, effect_slot = "common_effect",
            file_tag="CommonEffect",name_col = "Genes",name_imputed = "imputed"
        )
        # 5) Deregulogram (only if interaction tested AND factor has 2 levels)
        if (!("NA" %in% tests_interaction)) {
            # Interaction: effects with significant interaction
            output_interaction <- .pp_export_volcano_tables(
                effects = effects, tests = tests, path_output = path_output,
                experiment = experiment, effect_slot = "interaction_effect",
                file_tag = "InteractionEffect",plot_contrasts=tests_interaction,
                name_col = "Genes", name_imputed = "imputed", ...
            )
            info <- .pp_infer_interaction_factor(se_limma, tests_interaction)
            if (!is.null(info) && length(info$levels) == 2) {
                plot_deregulogram(se_limma = se_limma, tests = tests,
                                  tests_interaction = tests_interaction,
                                  factor_name = info$name,
                                  factor_levels = info$levels,
                                  path_output = path_output,
                                  experiment = experiment, alpha = 0.05,
                                  lfc = 1, label_col = "Genes"
                )
            }
        }
        return(list(se_raw = se0, se_filt = filt$se_filt,removed = filt$removed,
                    se_imp   = se_imp,
                    effects = list(output_main_all, output_interaction,
                                   output_common)
        ))
    }else{
        return(list(se_raw = se0, se_filt = filt$se_filt,removed = filt$removed,
                    se_imp   = se_imp
        ))
    }
}
