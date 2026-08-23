#' Define a PPBR dose-finding design
#'
#' @param dose Strictly increasing numeric dose labels.
#' @param target Target toxicity probability.
#' @param cohort_size Number of patients per cohort.
#' @param max_sample Maximum trial sample size.
#' @param overdose_threshold Toxicity probability defining an overdose.
#' @param safety_cutoff Posterior probability cutoff for eliminating a dose and
#'   all higher doses.
#' @param min_evidence Minimum patients required before elimination or terminal
#'   recommendation.
#' @param start_dose Dose used for the first cohort.
#' @param state_weight Prior probability assigned to the exact-target family.
#' @param prior_alpha,prior_beta Local beta-prior parameters.
#' @return An object of class `ppbr_design`.
#' @export
ppbr_design <- function(dose, target, cohort_size = 3L, max_sample,
                        overdose_threshold = target + 0.10,
                        safety_cutoff = 0.95,
                        min_evidence = cohort_size,
                        start_dose = dose[1L], state_weight = 0.50,
                        prior_alpha = 1, prior_beta = 1) {
  if (!is.numeric(dose) || length(dose) < 2L || any(!is.finite(dose)) ||
      any(diff(dose) <= 0))
    .ppbr_stop("`dose` must contain at least two strictly increasing numeric values.")
  target <- .check_probability(target, "target")
  overdose_threshold <- .check_probability(overdose_threshold, "overdose_threshold")
  if (overdose_threshold <= target)
    .ppbr_stop("`overdose_threshold` must be greater than `target`.")
  safety_cutoff <- .check_probability(safety_cutoff, "safety_cutoff")
  state_weight <- .check_probability(state_weight, "state_weight", open = FALSE)
  ints <- c(cohort_size, max_sample, min_evidence)
  if (any(!is.finite(ints)) || any(ints < 1) || any(abs(ints - round(ints)) > 1e-8))
    .ppbr_stop("`cohort_size`, `max_sample`, and `min_evidence` must be positive integers.")
  if (max_sample < cohort_size)
    .ppbr_stop("`max_sample` must be at least one cohort.")
  if (!is.finite(prior_alpha) || !is.finite(prior_beta) ||
      prior_alpha <= 0 || prior_beta <= 0)
    .ppbr_stop("Beta-prior parameters must be positive.")
  start_index <- .dose_index(dose, start_dose, "start_dose")
  structure(list(
    dose = as.numeric(dose), J = length(dose), target = target,
    cohort_size = as.integer(cohort_size), max_sample = as.integer(max_sample),
    overdose_threshold = overdose_threshold, safety_cutoff = safety_cutoff,
    min_evidence = as.integer(min_evidence), start_index = start_index,
    state_weight = state_weight, prior_alpha = prior_alpha,
    prior_beta = prior_beta), class = "ppbr_design")
}

#' @export
print.ppbr_design <- function(x, ...) {
  cat("PPBR dose-finding design\n")
  cat("  Doses:", paste(x$dose, collapse = ", "), "\n")
  cat("  Target / overdose threshold:", x$target, "/", x$overdose_threshold, "\n")
  cat("  Cohort size / maximum sample:", x$cohort_size, "/", x$max_sample, "\n")
  cat("  Starting dose:", x$dose[x$start_index], "\n")
  invisible(x)
}
