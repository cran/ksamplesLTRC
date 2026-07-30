#' Add estimators to an existing plot
#'
#' Adds the survival estimator to an existing plot.
#'
#' @param x An object of class `"estimators"`.
#' @param ... Additional graphical parameters passed to [graphics::lines()].
#'
#' @return Invisibly returns `x`.
#'
#' @method lines t.estimators
#'
#' @export
lines.t.estimators <- function(x, ...) {

  if (!inherits(x, "t.estimators")) {
    stop("Object is not of class 't.estimators'.", call. = FALSE)
  }

  lines(x$trunc.time, x$estimation, type = "s", ...)

  invisible(x)
}
