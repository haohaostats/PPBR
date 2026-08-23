#' Fit PPBR to the current dose-level data
#'
#' @param design A `ppbr_design` object.
#' @param n Number of patients treated at each dose.
#' @param dlt Number of dose-limiting toxicities at each dose.
#' @param current_dose Current dose label. If omitted, the last administered
#'   dose in the ordered grid is used; an empty trial uses the starting dose.
#' @return An object of class `ppbr_fit`.
#' @export
ppbr <- function(design, n, dlt, current_dose = NULL) {
  if (!inherits(design, "ppbr_design"))
    .ppbr_stop("`design` must be created by `ppbr_design()`.")
  n <- .check_count_vector(n, "n", design$J)
  dlt <- .check_count_vector(dlt, "dlt", design$J)
  if (any(dlt > n)) .ppbr_stop("Each `dlt` count must be no greater than `n`.")
  if (sum(n) > design$max_sample)
    .ppbr_stop("The total number treated exceeds `design$max_sample`.")
  if (is.null(current_dose)) {
    current_index <- if (sum(n) == 0L) design$start_index else max(which(n > 0L))
  } else current_index <- .dose_index(design$dose, current_dose)

  engine <- .new_ppbr_engine(design)
  safety <- .ppbr_safety(design, n, dlt)
  posterior <- engine$posterior(n, dlt)
  final_index <- engine$terminal(n, dlt, safety$usable)

  state_table <- data.frame(
    state = names(posterior$q),
    family = ifelse(substr(names(posterior$q), 1L, 1L) == "T",
                    "Exact target", "Gap"),
    probability = unname(posterior$q))
  target_table <- data.frame(
    index = 0:design$J,
    dose = c(NA_real_, design$dose),
    probability = unname(posterior$rho),
    eligible = c(TRUE, n >= design$min_evidence & safety$usable))
  bracket_table <- data.frame(
    lower = design$dose[-design$J], upper = design$dose[-1L],
    lower_favored = posterior$h)
  safety_table <- data.frame(
    dose = design$dose, n = n, dlt = dlt,
    overdose_probability = safety$probability,
    usable = safety$usable,
    terminally_eligible = n >= design$min_evidence & safety$usable)

  structure(list(
    design = design, n = n, dlt = dlt, current_index = current_index,
    current_dose = design$dose[current_index], posterior = posterior,
    states = state_table, targets = target_table, brackets = bracket_table,
    safety = safety_table, final_index = final_index,
    final_dose = if (final_index == 0L) NA_real_ else design$dose[final_index],
    engine = engine), class = "ppbr_fit")
}

#' @export
print.ppbr_fit <- function(x, ...) {
  cat("PPBR analysis\n")
  cat("  Enrolled:", sum(x$n), "of", x$design$max_sample, "\n")
  cat("  Current dose:", x$current_dose, "\n")
  cat("  Terminal recommendation:",
      if (x$final_index == 0L) "No dose" else x$final_dose, "\n")
  cat("  Posterior selection error:", sprintf("%.3f", x$posterior$risk), "\n")
  invisible(x)
}

#' @export
summary.ppbr_fit <- function(object, ...) {
  unrestricted <- as.integer(names(object$posterior$rho)[.lower_argmax(object$posterior$rho)])
  out <- list(
    enrolled = sum(object$n), maximum_sample = object$design$max_sample,
    current_dose = object$current_dose,
    unrestricted_dose = if (unrestricted == 0L) NA_real_ else object$design$dose[unrestricted],
    recommended_dose = object$final_dose,
    posterior_error = object$posterior$risk,
    targets = object$targets, brackets = object$brackets,
    safety = object$safety)
  class(out) <- "summary.ppbr_fit"
  out
}

#' @export
print.summary.ppbr_fit <- function(x, ...) {
  cat("PPBR posterior summary\n")
  cat("  Enrolled:", x$enrolled, "of", x$maximum_sample, "\n")
  cat("  Current dose:", x$current_dose, "\n")
  cat("  Unrestricted posterior-optimal dose:",
      ifelse(is.na(x$unrestricted_dose), "No dose", x$unrestricted_dose), "\n")
  cat("  Eligible terminal recommendation:",
      ifelse(is.na(x$recommended_dose), "No dose", x$recommended_dose), "\n")
  cat("  Posterior selection error:", sprintf("%.3f", x$posterior_error), "\n\n")
  cat("Target probabilities\n")
  tab <- x$targets
  tab$dose <- ifelse(is.na(tab$dose), "No dose", as.character(tab$dose))
  tab$probability <- round(tab$probability, 4)
  print(tab, row.names = FALSE)
  invisible(x)
}
