#' Plot method for k-sample test objects
#'
#' Plots the estimated survival curves for the groups involved in
#' [ltrc.test()].
#'
#' @param x An object of class `"ksamples_test"`.
#' @param y Ignored.
#' @param col Vector of colours for the groups.
#' @param lwd Line width.
#' @param lty Line type.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param main Plot title.
#' @param legend.pos Position of the legend.
#' @param ... Additional graphical parameters passed to [graphics::plot()].
#'
#' @return Invisibly returns `x`.
#'
#' @method plot ksamples_test
#' @export
plot.ksamples_test <- function(x, y = NULL,
                               col = NULL,
                               lwd = 2,
                               lty = 1,
                               xlab = "Time",
                               ylab = "Survival probability",
                               main = "Estimated survival curves",
                               legend.pos = "topright",
                               ...) {

  if (!inherits(x, "ksamples_test")) {
    stop("Object is not of class 'ksamples_test'.", call. = FALSE)
  }

  if (is.null(x$curves) || length(x$curves) == 0L) {
    stop("No curves stored in object. Modify 'ltrc.test()' so it returns them.",
         call. = FALSE)
  }

  K <- length(x$curves)

  if (is.null(col)) {
    col <- seq_len(K)
  }

  get_est <- function(obj) {
    if (!is.null(obj$estimation)) return(obj$estimation)
    if (!is.null(obj$estimate)) return(obj$estimate)
    stop("Stored curve does not contain 'estimation'/'estimate'.", call. = FALSE)
  }

  get_time <- function(obj) {
    if (!is.null(obj$fail.time)) return(obj$fail.time)
    if (!is.null(obj$time)) return(obj$time)
    stop("Stored curve does not contain 'fail.time'/'time'.", call. = FALSE)
  }

  all_times <- unlist(lapply(x$curves, get_time), use.names = FALSE)
  if (length(all_times) == 0L) {
    stop("No event times available to plot.", call. = FALSE)
  }

  graphics::plot(NA, NA,
                 xlim = range(all_times, finite = TRUE),
                 ylim = c(0, 1),
                 xlab = xlab, ylab = ylab, main = main,
                 ...)

  for (j in seq_len(K)) {
    tt <- get_time(x$curves[[j]])
    ss <- get_est(x$curves[[j]])

    tt0 <- c(min(tt), tt)
    ss0 <- c(1, ss)

    graphics::lines(tt0, ss0, type = "s",
                    col = col[j], lwd = lwd, lty = lty)
  }

  graphics::legend(legend.pos,
                   legend = paste("Group", x$groups),
                   col = col, lwd = lwd, lty = lty, bty = "n")

  invisible(x)
}
