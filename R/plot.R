#' Plot a PPBR analysis
#'
#' @param x A `ppbr_fit` object.
#' @param type One of `"target"`, `"states"`, `"brackets"`, `"safety"`, or
#'   `"all"`.
#' @param ... Additional graphical arguments.
#' @export
plot.ppbr_fit <- function(x, type = c("target", "states", "brackets", "safety", "all"), ...) {
  type <- match.arg(type)
  accent <- "#176B75"
  muted <- "#819198"
  one <- function(which) {
    if (which == "target") {
      labs <- c("No dose", as.character(x$design$dose))
      cols <- ifelse(x$targets$eligible, accent, "white")
      bp <- barplot(x$targets$probability,
                    names.arg = rep("", length(labs)), col = cols,
                    border = accent, ylab = "Posterior probability",
                    xlab = "Dose", main = "Discrete target posterior", ...)
      axis(1, at = bp, labels = c("None", x$design$dose),
           tick = FALSE, gap.axis = -1)
      invisible(bp)
    } else if (which == "states") {
      cols <- ifelse(x$states$family == "Exact target", accent, muted)
      barplot(x$states$probability, names.arg = x$states$state, col = cols,
              border = cols, las = 2, ylab = "Posterior probability",
              main = "Target-geometry posterior", ...)
      legend("topleft", c("Exact target", "Gap"), fill = c(accent, muted), bty = "n")
    } else if (which == "brackets") {
      labs <- paste(x$brackets$lower, x$brackets$upper, sep = "-")
      plot(seq_along(labs), x$brackets$lower_favored, type = "b", pch = 19,
           col = accent, ylim = c(0, 1), xaxt = "n",
           xlab = "Adjacent dose pair", ylab = "Probability favoring lower member",
           main = "Within-gap resolution", ...)
      axis(1, at = seq_along(labs), labels = labs)
      abline(h = 0.5, lty = 2, col = muted)
    } else {
      cols <- ifelse(x$safety$usable, accent, muted)
      barplot(x$safety$overdose_probability, names.arg = x$safety$dose,
              col = cols, border = cols, ylim = c(0, 1),
              xlab = "Dose", ylab = "Posterior overdose probability",
              main = "Safety posterior", ...)
      abline(h = x$design$safety_cutoff, lty = 2, col = muted)
    }
  }
  if (type == "all") {
    old <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
    on.exit(par(old), add = TRUE)
    for (z in c("states", "target", "brackets", "safety")) one(z)
  } else one(type)
  invisible(x)
}

#' Plot PPBR simulation results
#'
#' @param x A `ppbr_simulation` object.
#' @param type `"selection"` or `"allocation"`.
#' @param ... Additional graphical arguments.
#' @export
plot.ppbr_simulation <- function(x, type = c("selection", "allocation"), ...) {
  type <- match.arg(type)
  scenarios <- unique(x$selection$scenario)
  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE)
  nc <- min(3L, length(scenarios))
  nr <- ceiling(length(scenarios) / nc)
  par(mfrow = c(nr, nc), mar = c(4.2, 4.3, 2.8, 1),
      oma = c(0, 0, 2.2, 0))
  accent <- "#168A93"
  border <- "#24444F"
  if (type == "selection") {
    mat <- sapply(scenarios, function(s)
      x$selection$probability[x$selection$scenario == s])
    rownames(mat) <- c("No dose", as.character(x$design$dose))
    keep <- rowSums(mat) > 0
    mat <- mat[keep, , drop = FALSE]
    labs <- ifelse(rownames(mat) == "No dose", "None", rownames(mat))
    ymax <- max(mat) * 1.18
    for (i in seq_along(scenarios)) {
      values <- mat[, i]
      bp <- barplot(values, names.arg = rep("", length(labs)), col = accent,
                    border = border, ylim = c(0, ymax),
                    xlab = "Final selected dose",
                    ylab = if ((i - 1L) %% nc == 0L) "Selection probability" else "",
                    ...)
      axis(1, at = bp, labels = labs, tick = FALSE, gap.axis = -1)
      text(bp, values, labels = sprintf("%.2f", values), pos = 3,
           cex = 0.8, col = border, xpd = NA)
      mtext(scenarios[i], side = 3, adj = 0, line = 0.5,
            font = 2, cex = 0.95)
    }
    mtext("Terminal selection", side = 3, outer = TRUE,
          font = 2, cex = 1.15, line = 0.5)
  } else {
    mat <- sapply(scenarios, function(s)
      x$allocation$mean_proportion[x$allocation$scenario == s])
    rownames(mat) <- as.character(x$design$dose)
    labs <- rownames(mat)
    ymax <- max(mat) * 1.18
    for (i in seq_along(scenarios)) {
      values <- mat[, i]
      bp <- barplot(values, names.arg = rep("", length(labs)), col = accent,
                    border = border, ylim = c(0, ymax), xlab = "Assigned dose",
                    ylab = if ((i - 1L) %% nc == 0L) "Mean allocation proportion" else "",
                    ...)
      axis(1, at = bp, labels = labs, tick = FALSE, gap.axis = -1)
      text(bp, values, labels = sprintf("%.2f", values), pos = 3,
           cex = 0.8, col = border, xpd = NA)
      mtext(scenarios[i], side = 3, adj = 0, line = 0.5,
            font = 2, cex = 0.95)
    }
    mtext("Patient allocation", side = 3, outer = TRUE,
          font = 2, cex = 1.15, line = 0.5)
  }
  invisible(x)
}
