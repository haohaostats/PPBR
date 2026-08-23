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
      bp <- barplot(x$targets$probability, names.arg = labs, col = cols,
                    border = accent, ylab = "Posterior probability",
                    xlab = "Dose", main = "Discrete target posterior", ...)
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
  par(mar = c(4.5, 4.7, 3.5, 7.5), xpd = NA)
  if (type == "selection") {
    mat <- sapply(scenarios, function(s)
      x$selection$probability[x$selection$scenario == s])
    rownames(mat) <- c("No dose", as.character(x$design$dose))
    cols <- c("#9AA3AA", grDevices::hcl.colors(x$design$J, "Dark 3"))
    barplot(mat, beside = TRUE, names.arg = scenarios,
            col = cols, border = "#24444F",
            ylab = "Selection probability", main = "Terminal selection", ...)
    usr <- par("usr")
    legend(x = usr[2L] + 0.025 * diff(usr[1:2]), y = usr[4L],
           xjust = 0, yjust = 1, legend = rownames(mat),
           fill = cols, border = "#24444F", title = "Dose",
           bty = "n", cex = 0.9, x.intersp = 0.8, y.intersp = 1.1)
  } else {
    mat <- sapply(scenarios, function(s)
      x$allocation$mean_proportion[x$allocation$scenario == s])
    rownames(mat) <- as.character(x$design$dose)
    cols <- grDevices::hcl.colors(nrow(mat), "Dark 3")
    barplot(mat, beside = TRUE, names.arg = scenarios,
            col = cols, border = "#24444F",
            ylab = "Mean allocation proportion", main = "Patient allocation", ...)
    usr <- par("usr")
    legend(x = usr[2L] + 0.025 * diff(usr[1:2]), y = usr[4L],
           xjust = 0, yjust = 1, legend = rownames(mat),
           fill = cols, border = "#24444F", title = "Dose", bty = "n",
           cex = 0.9, x.intersp = 0.8, y.intersp = 1.1)
  }
  invisible(x)
}
