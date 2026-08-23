.new_ppbr_engine <- function(design) {
  J <- design$J
  phi <- design$target
  m <- design$cohort_size
  alpha0 <- design$prior_alpha
  beta0 <- design$prior_beta
  gl <- .gauss_legendre(24L)
  states <- c(paste0("T", seq_len(J)), paste0("G", 0:J))
  prior <- c(rep(design$state_weight / J, J),
             rep((1 - design$state_weight) / (J + 1L), J + 1L))
  cache <- new.env(hash = TRUE, parent = emptyenv())

  key <- function(n, y) paste(c(rbind(n, y)), collapse = ":")

  log_below <- function(n, y) {
    A <- alpha0 + y
    B <- beta0 + n - y
    lbeta(A, B) + pbeta(phi, A, B, log.p = TRUE) -
      lbeta(alpha0, beta0) - pbeta(phi, alpha0, beta0, log.p = TRUE)
  }

  log_above <- function(n, y) {
    A <- alpha0 + y
    B <- beta0 + n - y
    lbeta(A, B) + pbeta(phi, A, B, lower.tail = FALSE, log.p = TRUE) -
      lbeta(alpha0, beta0) -
      pbeta(phi, alpha0, beta0, lower.tail = FALSE, log.p = TRUE)
  }

  log_exact <- function(n, y) y * log(phi) + (n - y) * log1p(-phi)

  gap_resolution <- function(b, n, y) {
    A1 <- alpha0 + y[b]
    B1 <- beta0 + n[b] - y[b]
    A2 <- alpha0 + y[b + 1L]
    B2 <- beta0 + n[b + 1L] - y[b + 1L]
    lo <- max(0, 2 * phi - 1)
    u <- (phi - lo) * (gl$x + 1) / 2 + lo
    den1 <- pbeta(phi, A1, B1)
    den2 <- pbeta(phi, A2, B2, lower.tail = FALSE)
    f1 <- dbeta(u, A1, B1) / pmax(den1, 1e-300)
    s2 <- pbeta(2 * phi - u, A2, B2, lower.tail = FALSE) /
      pmax(den2, 1e-300)
    .clamp_prob((phi - lo) * sum(gl$w * f1 * s2) / 2, 1e-10)
  }

  posterior <- function(n, y) {
    ky <- key(n, y)
    if (exists(ky, cache, inherits = FALSE))
      return(get(ky, cache, inherits = FALSE))
    lm <- log_below(n, y)
    l0 <- log_exact(n, y)
    lp <- log_above(n, y)
    logs <- numeric(2L * J + 1L)
    for (j in seq_len(J)) {
      below <- if (j > 1L) sum(lm[seq_len(j - 1L)]) else 0
      above <- if (j < J) sum(lp[(j + 1L):J]) else 0
      logs[j] <- below + l0[j] + above
    }
    for (b in 0:J) {
      below <- if (b > 0L) sum(lm[seq_len(b)]) else 0
      above <- if (b < J) sum(lp[(b + 1L):J]) else 0
      logs[J + 1L + b] <- below + above
    }
    lpost <- log(prior) + logs
    q <- exp(lpost - .log_sum_exp(lpost))
    names(q) <- states
    h <- if (J > 1L)
      vapply(seq_len(J - 1L), gap_resolution, numeric(1), n = n, y = y) else numeric()

    A1 <- alpha0 + y[1L]
    B1 <- beta0 + n[1L] - y[1L]
    den <- pbeta(phi, A1, B1, lower.tail = FALSE)
    num <- pbeta(design$overdose_threshold, A1, B1) - pbeta(phi, A1, B1)
    edge_safe <- min(1, max(0, num / max(den, 1e-300)))

    rho <- numeric(J + 1L)
    names(rho) <- as.character(0:J)
    rho[1L] <- q["G0"] * (1 - edge_safe)
    for (j in seq_len(J)) {
      rho[j + 1L] <- q[paste0("T", j)]
      if (j == J) rho[j + 1L] <- rho[j + 1L] + q[paste0("G", J)]
      if (j <= J - 1L)
        rho[j + 1L] <- rho[j + 1L] + q[paste0("G", j)] * h[j]
      if (j >= 2L)
        rho[j + 1L] <- rho[j + 1L] + q[paste0("G", j - 1L)] * (1 - h[j - 1L])
    }
    rho[2L] <- rho[2L] + q["G0"] * edge_safe
    rho <- rho / sum(rho)

    kappa <- numeric(J + 1L)
    names(kappa) <- as.character(0:J)
    kappa[1L] <- q["G0"]
    if (J > 1L) for (b in seq_len(J - 1L))
      kappa[b + 1L] <- q[paste0("G", b)] + q[paste0("T", b)] +
        q[paste0("T", b + 1L)]
    kappa[J + 1L] <- q[paste0("G", J)]

    out <- list(q = q, h = h, rho = rho, kappa = kappa,
                edge_safe = edge_safe, risk = 1 - max(rho))
    assign(ky, out, cache)
    out
  }

  class_at <- function(state, a) {
    if (substr(state, 1L, 1L) == "T") {
      j <- as.integer(sub("T", "", state))
      if (a < j) "-" else if (a == j) "0" else "+"
    } else {
      b <- as.integer(sub("G", "", state))
      if (a <= b) "-" else "+"
    }
  }

  predictive_component <- function(A, B, z, cls) {
    if (cls == "0") return(dbinom(z, m, phi))
    if (cls == "-") {
      logp <- lchoose(m, z) + lbeta(A + z, B + m - z) +
        pbeta(phi, A + z, B + m - z, log.p = TRUE) -
        lbeta(A, B) - pbeta(phi, A, B, log.p = TRUE)
    } else {
      logp <- lchoose(m, z) + lbeta(A + z, B + m - z) +
        pbeta(phi, A + z, B + m - z, lower.tail = FALSE, log.p = TRUE) -
        lbeta(A, B) - pbeta(phi, A, B, lower.tail = FALSE, log.p = TRUE)
    }
    exp(logp)
  }

  predictive_detail <- function(a, n, y, post = posterior(n, y)) {
    A <- alpha0 + y[a]
    B <- beta0 + n[a] - y[a]
    cls <- vapply(states, class_at, character(1), a = a)
    pz <- updated_risk <- numeric(m + 1L)
    for (z in 0:m) {
      ps <- vapply(cls, function(cc) predictive_component(A, B, z, cc), numeric(1))
      pz[z + 1L] <- sum(post$q * ps)
      nn <- n
      yy <- y
      nn[a] <- nn[a] + m
      yy[a] <- yy[a] + z
      updated_risk[z + 1L] <- posterior(nn, yy)$risk
    }
    pz <- pz / sum(pz)
    data.frame(dlt = 0:m, probability = pz, updated_risk = updated_risk,
               expected_risk = sum(pz * updated_risk))
  }

  candidates <- function(n, y, current, usable) {
    rate <- if (n[current] > 0L) y[current] / n[current] else 0
    cand <- if (rate < phi) c(current, current + 1L) else
      if (rate > phi) c(current - 1L, current) else current
    reverse <- if (rate < phi) current - 1L else
      if (rate > phi) current + 1L else integer()
    reverse <- reverse[reverse >= 1L & reverse <= J]
    reverse <- reverse[n[reverse] >= m & usable[reverse]]
    cand <- unique(c(cand, reverse))
    cand <- cand[cand >= 1L & cand <= J]
    cand[usable[cand]]
  }

  choose_action <- function(n, y, current, usable) {
    cand <- candidates(n, y, current, usable)
    if (!length(cand)) return(list(index = NA_integer_, candidates = integer(),
                                   risks = numeric()))
    details <- lapply(cand, predictive_detail, n = n, y = y)
    risks <- vapply(details, function(z) z$expected_risk[1L], numeric(1))
    names(risks) <- as.character(cand)
    list(index = cand[.lower_argmin(risks)], candidates = cand,
         risks = risks, details = setNames(details, as.character(cand)))
  }

  terminal <- function(n, y, usable) {
    post <- posterior(n, y)
    eligible <- c(TRUE, n >= design$min_evidence & usable)
    feasible <- 0L
    if (J > 1L) for (b in seq_len(J - 1L))
      if (eligible[b + 1L] || eligible[b + 2L]) feasible <- c(feasible, b)
    if (eligible[J + 1L]) feasible <- c(feasible, J)
    bhat <- feasible[.lower_argmax(post$kappa[as.character(feasible)])]
    if (bhat == 0L)
      return(if (post$edge_safe >= 0.5 && eligible[2L]) 1L else 0L)
    if (bhat == J) return(J)
    members <- c(bhat, bhat + 1L)
    members <- members[eligible[members + 1L]]
    if (!length(members)) return(0L)
    if (length(members) == 1L) return(members)
    if (post$h[bhat] >= 0.5) bhat else bhat + 1L
  }

  list(posterior = posterior, predictive_detail = predictive_detail,
       choose_action = choose_action, terminal = terminal, states = states)
}
