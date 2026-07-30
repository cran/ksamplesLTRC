#' Plot method for estimators objects
#'
#' Produces a step plot of the survival estimator.
#'
#' @param x An object of class `"estimators"`.
#' @param ... Additional graphical parameters passed to [graphics::plot()].
#' @param conf.int Logical, default to FALSE. If TRUE, pointwise confidence bands
#' are depicted for the curves.
#'
#' @return Invisibly returns `x`.
#'
#' @method plot estimators
#'
#' @export
plot.estimators <- function(x,conf.int = FALSE, ...) {

  if (!inherits(x, "estimators")) {
    stop("Object is not of class 'estimators'.", call. = FALSE)
  }

  plot(x$fail.time, x$estimation, type = "s", ylim = c(0,1), ...)

  if (!is.null(x$conf.limits)) {
    lines(x$fail.time, x$conf.limits[,1], lty = 2, type = "s",...)
    lines(x$fail.time, x$conf.limits[,2], lty = 2, type = "s",...)
  }


  invisible(x)
}
