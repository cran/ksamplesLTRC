#' Add estimators to an existing plot
#'
#' Adds the survival estimator to an existing plot.
#'
#' @param x An object of class `"estimators"`.
#' @param ... Additional graphical parameters passed to [graphics::lines()].
#' @param conf.int Logical, default to FALSE. If TRUE, pointwise confidence bands
#' are depicted for the curves.
#'
#' @return Invisibly returns `x`.
#'
#' @method lines estimators
#'
#' @export
lines.estimators <- function(x, conf.int = FALSE, ...) {

  if (!inherits(x, "estimators")) {
    stop("Object is not of class 'estimators'.", call. = FALSE)
  }

  lines(x$fail.time, x$estimation, type = "s", ...)


  if (!is.null(x$conf.limits)) {
    lines(x$fail.time, x$conf.limits[,1], lty = 2, type = "s", ...)
    lines(x$fail.time, x$conf.limits[,2], lty = 2, type = "s", ...)
  }


  invisible(x)
}
