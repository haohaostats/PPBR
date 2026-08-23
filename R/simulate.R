.true_target_index <- function(pi, design) {
  if (pi[1L] > design$overdose_threshold) return(0L)
  .lower_argmin(abs(pi - design$target))
}

.simulate_one_ppbr <- function(design, pi, engine) {
  n <- dlt <- integer(design$J)
  current <- design$start_index
  stopped <- FALSE
  while (sum(n) < design$max_sample) {
    take <- min(design$cohort_size, design$max_sample - sum(n))
    yy <- sum(runif(take) < pi[current])
    n[current] <- n[current] + take
    dlt[current] <- dlt[current] + yy
    safety <- .ppbr_safety(design, n, dlt)
    if (!safety$usable[1L]) {
      stopped <- TRUE
      break
    }
    if (sum(n) >= design$max_sample) break
    if (!safety$usable[current]) {
      lower <- which(safety$usable & seq_len(design$J) < current)
      if (!length(lower)) {
        stopped <- TRUE
        break
      }
      current <- max(lower)
    } else {
      action <- engine$choose_action(n, dlt, current, safety$usable)
      if (is.na(action$index)) {
        stopped <- TRUE
        break
      }
      current <- action$index
    }
  }
  safety <- .ppbr_safety(design, n, dlt)
  selected <- if (stopped) 0L else engine$terminal(n, dlt, safety$usable)
  list(selected = selected, n = n, dlt = dlt, enrolled = sum(n), stopped = stopped)
}

#' Simulate operating characteristics for user-specified scenarios
#'
#' @param design A `ppbr_design` object.
#' @param scenarios A numeric vector, matrix, or list. Each scenario must contain
#'   one nondecreasing toxicity probability per dose.
#' @param nsim Number of simulated trials per scenario.
#' @param seed Optional random seed.
#' @return An object of class `ppbr_simulation`.
#' @export
ppbr_simulate <- function(design, scenarios, nsim = 1000L, seed = NULL) {
  if (!inherits(design, "ppbr_design"))
    .ppbr_stop("`design` must be created by `ppbr_design()`.")
  if (is.list(scenarios) && !is.data.frame(scenarios)) {
    scenario_names <- names(scenarios)
    scenarios <- do.call(rbind, lapply(scenarios, as.numeric))
    if (!is.null(scenario_names) && all(nzchar(scenario_names)))
      rownames(scenarios) <- scenario_names
  }
  if (is.null(dim(scenarios))) scenarios <- matrix(as.numeric(scenarios), nrow = 1L)
  scenarios <- as.matrix(scenarios)
  storage.mode(scenarios) <- "double"
  if (ncol(scenarios) != design$J || any(!is.finite(scenarios)) ||
      any(scenarios <= 0) || any(scenarios >= 1))
    .ppbr_stop("Each scenario must contain one probability strictly between 0 and 1 per dose.")
  if (any(apply(scenarios, 1L, function(z) any(diff(z) < 0))))
    .ppbr_stop("Each toxicity scenario must be nondecreasing across doses.")
  if (length(nsim) != 1L || !is.finite(nsim) || nsim < 1L || nsim != as.integer(nsim))
    .ppbr_stop("`nsim` must be a positive integer.")
  nsim <- as.integer(nsim)
  if (!is.null(seed)) set.seed(seed)
  if (is.null(rownames(scenarios))) rownames(scenarios) <- paste0("Scenario ", seq_len(nrow(scenarios)))

  engine <- .new_ppbr_engine(design)
  trials <- vector("list", nrow(scenarios) * nsim)
  k <- 0L
  for (s in seq_len(nrow(scenarios))) for (r in seq_len(nsim)) {
    z <- .simulate_one_ppbr(design, scenarios[s, ], engine)
    k <- k + 1L
    trials[[k]] <- data.frame(
      scenario = rownames(scenarios)[s], replicate = r,
      selected = z$selected, enrolled = z$enrolled,
      total_dlt = sum(z$dlt), stopped = z$stopped,
      as.list(setNames(z$n, paste0("n", seq_len(design$J)))))
  }
  trials <- do.call(rbind, trials)
  selection <- do.call(rbind, lapply(unique(trials$scenario), function(s) {
    z <- trials[trials$scenario == s, ]
    data.frame(scenario = s, index = 0:design$J,
      dose = c(NA_real_, design$dose),
      probability = tabulate(z$selected + 1L, nbins = design$J + 1L) / nrow(z))
  }))
  allocation <- do.call(rbind, lapply(unique(trials$scenario), function(s) {
    z <- trials[trials$scenario == s, , drop = FALSE]
    nn <- colMeans(z[paste0("n", seq_len(design$J))])
    data.frame(scenario = s, dose = design$dose,
               mean_patients = as.numeric(nn),
               mean_proportion = as.numeric(nn) / mean(z$enrolled))
  }))
  scenario_summary <- do.call(rbind, lapply(seq_len(nrow(scenarios)), function(s) {
    label <- rownames(scenarios)[s]
    z <- trials[trials$scenario == label, ]
    jstar <- .true_target_index(scenarios[s, ], design)
    data.frame(scenario = label, target_index = jstar,
      target_dose = if (jstar == 0L) NA_real_ else design$dose[jstar],
      correct_selection = mean(z$selected == jstar),
      no_dose = mean(z$selected == 0L), mean_enrolled = mean(z$enrolled),
      mean_dlt = mean(z$total_dlt), stopped = mean(z$stopped))
  }))
  structure(list(design = design, scenarios = scenarios, nsim = nsim,
    trials = trials, selection = selection, allocation = allocation,
    summary = scenario_summary), class = "ppbr_simulation")
}

#' @export
print.ppbr_simulation <- function(x, ...) {
  cat("PPBR operating-characteristic simulation\n")
  cat("  Scenarios:", nrow(x$scenarios), "\n")
  cat("  Replicates per scenario:", x$nsim, "\n")
  invisible(x)
}

#' @export
summary.ppbr_simulation <- function(object, ...) {
  z <- object$summary
  z$correct_selection <- round(z$correct_selection, 4)
  z$no_dose <- round(z$no_dose, 4)
  z$mean_enrolled <- round(z$mean_enrolled, 2)
  z$mean_dlt <- round(z$mean_dlt, 2)
  z$stopped <- round(z$stopped, 4)
  z
}
