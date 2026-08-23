.ppbr_stop <- function(message, call. = FALSE) stop(message, call. = call.)

.check_probability <- function(x, name, open = TRUE) {
  ok <- length(x) == 1L && is.finite(x)
  ok <- ok && if (open) x > 0 && x < 1 else x >= 0 && x <= 1
  if (!ok) .ppbr_stop(sprintf("`%s` must be a single probability %s 0 and 1.",
                              name, if (open) "strictly between" else "between"))
  as.numeric(x)
}

.check_count_vector <- function(x, name, J) {
  if (!is.numeric(x) || length(x) != J || any(!is.finite(x)) ||
      any(x < 0) || any(abs(x - round(x)) > 1e-8))
    .ppbr_stop(sprintf("`%s` must contain %d non-negative integer counts.", name, J))
  as.integer(round(x))
}

.dose_index <- function(dose, value, name = "current_dose") {
  idx <- match(value, dose)
  if (length(value) != 1L || is.na(idx))
    .ppbr_stop(sprintf("`%s` must equal one of the values in `design$dose`.", name))
  as.integer(idx)
}

.lower_argmin <- function(x) which(x == min(x, na.rm = TRUE))[1L]
.lower_argmax <- function(x) which(x == max(x, na.rm = TRUE))[1L]
.clamp_prob <- function(x, eps = 1e-12) pmin(1 - eps, pmax(eps, x))

.log_sum_exp <- function(x) {
  z <- max(x)
  z + log(sum(exp(x - z)))
}

.gauss_legendre <- function(n = 24L) {
  i <- seq_len(n - 1L)
  beta <- i / sqrt(4 * i^2 - 1)
  jacobi <- matrix(0, n, n)
  jacobi[cbind(seq_len(n - 1L), 2:n)] <- beta
  jacobi[cbind(2:n, seq_len(n - 1L))] <- beta
  ev <- eigen(jacobi, symmetric = TRUE)
  ord <- order(ev$values)
  list(x = ev$values[ord], w = 2 * ev$vectors[1, ord]^2)
}

.ppbr_safety <- function(design, n, dlt) {
  xi <- pbeta(design$overdose_threshold, dlt + design$prior_alpha,
              n - dlt + design$prior_beta, lower.tail = FALSE)
  bad <- which(n >= design$min_evidence & xi > design$safety_cutoff)
  usable <- rep(TRUE, design$J)
  if (length(bad)) usable[min(bad):design$J] <- FALSE
  list(probability = xi, usable = usable,
       eliminated_from = if (length(bad)) min(bad) else NA_integer_)
}
