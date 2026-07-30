#' Survival function estimator under left truncation and right
#' censoring
#'
#' Computes the survival function estimator accounting for censored observations.
#' It also returns pointwise confidence bands.
#'
#' @param o.sample Matrix or data frame with columns (U, Y, Delta, group).
#' @param conf.int Logical, set as default to FALSE. If set to TRUE, the function
#' returns a pointwise confidence band.
#' @param alpha Pointwise confidence bands with \eqn{1-\alpha} confidence level
#' @param conf.type Type of pointwise confidence band. Set as default to \eqn{plain}.
#'
#' @return An object with the following content:
#'
#' fail.time: ordered distinct observed event times
#'
#' estimation: estimated survival function at each observed event time
#'
#' prob: jumps of the estimator
#'
#' If conf.int is set to TRUE, it also returns the confidence band lower and upper limits
#'
#' @details
#' The estimator is defined as
#' \deqn{ S_n(t) = \prod_{y_i \leq t} \left( 1-\frac{d_i}{n_i} \right), }
#' where \eqn{y_1, \ldots, y_m} are the distinct observed event times,
#' \eqn{d_i = \# \{j : Y_j = y_i , \Delta_j=1 \}} and \eqn{n_i = \#\{j : U_j \leq y_i \leq Y_j\}}. It reduces to the
#' Kaplan--Meier estimator if data are not subject to truncation, and to the Lynden-Bell
#' estimator if no censoring is present.
#'
#' The variance for the confidence bands follows a Greenwood's like formula:
#' \deqn{ Var(S_n(t)) = S_n(t)^2 \sum_{y_i \leq t} \frac{d_i}{n_i(n_i-d_i)}.}
#' Different transformations are also possible.
#
#' @references
#' Klein, J. P. and Moeschberger, M. L. (2003). Survival Analysis: Techniques for Censored and
#' Truncated Data. Springer
#'
#' Wang, M.C. (1991). Nonparametric estimation from cross-sectional survival data. Journal of
#' the American Statistical Association, 86:130-143.
#'
#' Stute, W. and Wang, J. L. (2008). The central limit theorem under random truncation. Bernoulli,
#' 14:604–622.
#'
#' @examples
#' x <- rnorm(500,4,1)
#' c <- rnorm(500,4,1)
#' o.sample <- data.frame(y = pmin(x,c),delta = as.numeric(x <= c))
#' plot(rc.estimator(o.sample,conf.int=TRUE,conf.type='loglog'))
#' curve(1-pnorm(x,4,1),add=TRUE,col=2)
#'
#' @export

rc.estimator <- function(o.sample,
                         conf.int = FALSE,
                         alpha = 0.05,
                         conf.type = c("plain", "loglog")) {

  conf.type <- match.arg(conf.type)

  if (!is.matrix(o.sample) && !is.data.frame(o.sample)) {
    stop("'o.sample' must be a matrix or data.frame.")
  }

  if (ncol(o.sample) < 2L) {
    stop("'o.sample' must have at least 2 columns.")
  }

  o.sample <- as.data.frame(o.sample)

  Y <- as.numeric(o.sample[[1]])
  delta <- as.numeric(o.sample[[2]])

  if (any(!is.finite(Y))) {
    stop("Observed times must be finite.")
  }

  if (any(delta != 0 & delta != 1)) {
    stop("The censoring indicator must take values 0/1.")
  }

  dist.fail <- sort(unique(Y[delta == 1]))
  m <- length(dist.fail)

  if (m == 0L) {
    lista <- list(
      fail.time     = numeric(0),
      estimation    = numeric(0),
      prob          = numeric(0),
      conf.limits   = NULL,
      conf.limits.F = NULL
    )
    class(lista) <- "estimators"
    return(lista)
  }

  idx <- match(Y[delta == 1], dist.fail)
  di <- tabulate(idx, nbins = m)

  Yo <- sort(Y)
  n <- length(Y)

  eps <- (abs(dist.fail) + 1) * .Machine$double.eps
  before <- findInterval(dist.fail - eps, Yo, rightmost.closed = TRUE)
  ni <- n - before

  if (any(ni <= 0)) {
    stop("Non-positive risk set size encountered; check data consistency.")
  }

  if (any(di > ni)) {
    stop("Found di > ni at some failure time; check data consistency.")
  }

  estimation <- cumprod(1 - di / ni)

  cum.prod <- c(1, estimation)
  prob <- cum.prod[1:m] - cum.prod[2:(m + 1)]

  var.i.u <- var.i.l <- rep(NA_real_, m)
  varF.u <- varF.l <- rep(NA_real_, m)

  if (conf.int) {

    if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) {
      stop("'alpha' must be in (0, 1).")
    }

    z.alpha <- qnorm(1 - alpha / 2)

    gterm <- di / (ni * (ni - di))
    gterm[!is.finite(gterm)] <- 0
    cum.sum <- cumsum(gterm)

    if (conf.type == "plain") {

      seS <- estimation * sqrt(cum.sum)
      var.i.l <- pmax(estimation - z.alpha * seS, 0)
      var.i.u <- pmin(estimation + z.alpha * seS, 1)

    } else {

      var.i.l <- var.i.u <- numeric(m)

      ok <- (estimation > 0) & (estimation < 1)

      var.i.l[!ok & estimation >= 1] <- 1
      var.i.u[!ok & estimation >= 1] <- 1
      var.i.l[!ok & estimation <= 0] <- 0
      var.i.u[!ok & estimation <= 0] <- 0

      w <- log(-log(estimation[ok]))
      se.w <- sqrt(cum.sum[ok]) / abs(log(estimation[ok]))

      var.i.l[ok] <- exp(-exp(w + z.alpha * se.w))
      var.i.u[ok] <- exp(-exp(w - z.alpha * se.w))
    }

    varF.l <- 1 - var.i.u
    varF.u <- 1 - var.i.l
  }

  lista <- list(
    fail.time     = dist.fail,
    estimation    = estimation,
    prob          = prob,
    conf.limits   = if (conf.int) cbind(var.i.l, var.i.u) else NULL,
    conf.limits.F = if (conf.int) cbind(varF.l, varF.u) else NULL
  )
  class(lista) <- "estimators"
  return(lista)
}
