#' Obtain the next PPBR action
#'
#' @param fit A fitted `ppbr_fit` object.
#' @return An object of class `ppbr_decision`.
#' @export
ppbr_next <- function(fit) {
  if (!inherits(fit, "ppbr_fit")) .ppbr_stop("`fit` must be created by `ppbr()`.")
  design <- fit$design
  total <- sum(fit$n)
  usable <- fit$safety$usable
  if (!usable[1L]) {
    return(structure(list(continue = FALSE, reason = "Lowest dose eliminated",
      next_index = NA_integer_, next_dose = NA_real_, final_index = 0L,
      final_dose = NA_real_, candidates = data.frame()), class = "ppbr_decision"))
  }
  if (total >= design$max_sample) {
    return(structure(list(continue = FALSE, reason = "Maximum sample reached",
      next_index = NA_integer_, next_dose = NA_real_, final_index = fit$final_index,
      final_dose = fit$final_dose, candidates = data.frame()), class = "ppbr_decision"))
  }
  if (total == 0L) {
    return(structure(list(continue = TRUE, reason = "Start trial",
      next_index = design$start_index, next_dose = design$dose[design$start_index],
      final_index = NA_integer_, final_dose = NA_real_, candidates = data.frame()),
      class = "ppbr_decision"))
  }
  current <- fit$current_index
  if (!usable[current]) {
    lower <- which(usable & seq_len(design$J) < current)
    if (!length(lower)) return(structure(list(continue = FALSE,
      reason = "No usable lower dose", next_index = NA_integer_, next_dose = NA_real_,
      final_index = 0L, final_dose = NA_real_, candidates = data.frame()),
      class = "ppbr_decision"))
    idx <- max(lower)
    return(structure(list(continue = TRUE, reason = "Current dose eliminated",
      next_index = idx, next_dose = design$dose[idx], final_index = NA_integer_,
      final_dose = NA_real_, candidates = data.frame()), class = "ppbr_decision"))
  }
  action <- fit$engine$choose_action(fit$n, fit$dlt, current, usable)
  if (is.na(action$index)) return(structure(list(continue = FALSE,
    reason = "No admissible candidate", next_index = NA_integer_, next_dose = NA_real_,
    final_index = fit$final_index, final_dose = fit$final_dose,
    candidates = data.frame()), class = "ppbr_decision"))
  cand <- data.frame(index = action$candidates,
                     dose = design$dose[action$candidates],
                     expected_risk = unname(action$risks))
  structure(list(continue = TRUE, reason = "Posterior predictive risk",
    next_index = action$index, next_dose = design$dose[action$index],
    final_index = NA_integer_, final_dose = NA_real_, candidates = cand,
    predictive = action$details), class = "ppbr_decision")
}

#' @export
print.ppbr_decision <- function(x, ...) {
  cat("PPBR decision\n")
  cat("  Continue:", if (x$continue) "Yes" else "No", "\n")
  cat("  Reason:", x$reason, "\n")
  if (x$continue) cat("  Next dose:", x$next_dose, "\n")
  if (!x$continue) cat("  Final recommendation:",
                       ifelse(is.na(x$final_dose), "No dose", x$final_dose), "\n")
  if (nrow(x$candidates)) {
    cat("\nCandidate risks\n")
    z <- x$candidates
    z$expected_risk <- round(z$expected_risk, 4)
    print(z, row.names = FALSE)
  }
  invisible(x)
}
