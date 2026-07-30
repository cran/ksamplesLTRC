#' Plot method for estimators objects
#'
#' Produces a step plot of the survival estimator.
#'
#' @param x An object of class `"estimators"`.
#' @param ... Additional graphical parameters passed to [graphics::plot()].
#'
#' @return Invisibly returns `x`.
#'
#' @method plot t.estimators
#'
#' @export
plot.t.estimators <- function(x, ...) {

  if (!inherits(x, "t.estimators")) {
    stop("Object is not of class 't.estimators'.", call. = FALSE)
  }

  graphics::plot.default(
    x$trunc.time,
    x$estimation,
    type = "s",
    ylim = c(0,1),
    ...
  )

  invisible(x)
}
