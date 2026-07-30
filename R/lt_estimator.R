#' Survival function estimator under left truncation
#'
#' Computes the survival function estimator accounting for left truncation.
#' It also returns pointwise confidence bands and includes the possibility of dealing
#' with holes in the sample.
#'
#' @param o.sample Matrix or data frame with columns (U, X), such that \eqn{U \le X}. If
#' this input contains more columns, they will not be used.
#' @param conf.int Logical, set as default to FALSE. If TRUE, the function
#' returns a pointwise confidence band.
#' @param alpha Pointwise confidence bands with \eqn{1-\alpha} confidence level
#' @param conf.type One of 'plain' (default) or 'loglog'. Type of pointwise confidence band.
#' @param holes One of 'add', 'conditional', or 'holes' (default) (see details).
#'
#' @return An object with the following content:
#'
#' * `fail.time`: ordered distinct observed event times
#'
#' * `estimation`: estimated survival function at each observed event time
#'
#' * `prob`: jumps of the estimator
#'
#' * If `conf.int = TRUE`, it also returns the confidence band lower and upper limits
#' for both the distribution and survival function.
#'
#' @details
#' The estimator of the survival function when data are subject to left truncation is
#' defined as (Lynden-Bell, 1971)
#' \deqn{ S_n(t) = \prod_{x_i \leq t} \left( 1-\frac{d_i}{r_i} \right), }
#' where \eqn{x_1, \ldots, x_m} are the distinct observed event times,
#' \eqn{d_i = \# \{j : X_j = x_i\}} and \eqn{r_i = \#\{j : U_j \leq x_i \leq X_j\}}. One
#' should note that this definition of the risk sets is not equal to that in `survival::survfit`,
#' where \eqn{r_i = \#\{j : U_j < x_i \leq X_j\}}, so these two implementations may not be
#' equal when ties are present in the data.
#'
#' The variance for the confidence bands follows a Greenwood's-like formula:
#' \deqn{ Var(S_n(t)) = S_n(t)^2 \sum_{x_i \leq t} \frac{d_i}{n_i(n_i-d_i)}.}
#' A log-log transformation is also possible. The confidence band is based on the
#' Gaussianity of the process \eqn{\sqrt{n}(S_n(t) - S(t))}.
#'
#' The function also includes a way to deal with possible holes. A hole is an event time with associated
#' unitary risk set, namely \eqn{r_i=1} for some \eqn{ i \in \{1, \ldots, m-1\} }. We exclude the largest
#' event time, since it is always a hole. If a hole is present
#' in a sample, the estimator degenerates and takes the value 1 from that time on. Different techniques
#' are available in the literature for the treatment of holes. For instance, in Stute and Wang (2008),
#' authors suggest to increase the risk set by 1 (holes='add').
#' Other references, like Klein and Moeschberger (2003), suggest conditioning on the largest hole different
#' from the largest event time (holes='conditional'). The possibility of computing the degenerate statistic
#' is also included (holes='holes').
#'
#
#' @references
#' Klein, J. P. and Moeschberger, M. L. (2003). Survival Analysis: Techniques for Censored and
#' Truncated Data. Springer.
#'
#' Lynden-Bell, D. (1971). A method of allowing for known observational selection in small
#' samples applied to 3CR quasars. Monthly Notices of the Royal Astronomical Society, 155:95-118.
#'
#' Stute, W. and Wang, J. L. (2008). The central limit theorem under random truncation. Bernoulli,
#' 14:604–622.
#'
#' @examples
#' if (requireNamespace("mvna", quietly = TRUE)) {
#' data("abortion", package = "mvna")
#' pooled.sample <- abortion[abortion$cause==3, 2:4]
#' plot( lt.estimator(pooled.sample[pooled.sample$group==0,],conf.int=TRUE,conf.type='loglog'),
#' xlab='Time (in weeks)',ylab='Estimated survival probability' )
#' lines( lt.estimator(pooled.sample[pooled.sample$group==1,],conf.int=TRUE,conf.type='loglog'),col=2 )
#' }
#' @export

lt.estimator <- function(
    o.sample,
    conf.int = FALSE,
    alpha = 0.05,
    conf.type = c("plain", "loglog"),
    holes = c("holes", "add", "conditional")
) {

  conf.type <- match.arg(conf.type)
  holes <- match.arg(holes)

  if (!is.matrix(o.sample) && !is.data.frame(o.sample)) {
    stop("'o.sample' must be a matrix or data.frame.")
  }

  if (ncol(o.sample) < 2L) {
    stop("'o.sample' must have at least 2 columns: (U, X).")
  }

  o.sample <- data.frame(
    U = as.numeric(o.sample[, 1L]),
    X = as.numeric(o.sample[, 2L])
  )

  if (anyNA(o.sample)) {
    stop("'o.sample' must not contain missing values.")
  }

  o.sample <- as.matrix(o.sample)

  empty_out <- function() {
    out <- list(
      fail.time = numeric(0),
      estimation = numeric(0),
      prob = numeric(0),
      conf.limits = NULL,
      conf.limits.F = NULL,
      hole = NULL
    )

    class(out) <- "estimators"
    out
  }

  compute_dn <- function(x) {
    dist.fail <- sort(unique(x[, 2]))
    m <- length(dist.fail)

    if (m == 0L) {
      return(list(
        dist.fail = numeric(0),
        m = 0L,
        di = numeric(0),
        ni = numeric(0)
      ))
    }

    di <- vapply(
      dist.fail,
      function(t) sum(x[, 2] == t),
      numeric(1)
    )

    ni <- vapply(
      dist.fail,
      function(t) sum(x[, 1] <= t & t <= x[, 2]),
      numeric(1)
    )

    list(
      dist.fail = dist.fail,
      m = m,
      di = di,
      ni = ni
    )
  }

  eligible_holes <- function(ni) {

    m <- length(ni)

    if (m <= 1L) {
      return(integer(0))
    }

    later_gt1 <- c(
      rev(cummax(rev(ni > 1L)))[-1L],
      FALSE
    )

    which(ni == 1L & later_gt1)
  }

  out_dn <- compute_dn(o.sample)

  if (out_dn$m == 0L) {
    return(empty_out())
  }

  hole_idx <- eligible_holes(out_dn$ni)

  largest.hole <- if (length(hole_idx) > 0L) {
    max(out_dn$dist.fail[hole_idx])
  } else {
    NULL
  }

  if (length(hole_idx) > 0L) {

    if (holes == "add") {

      out_dn$ni <- out_dn$ni + 1L

    } else if (holes == "conditional") {

      o.sample <- o.sample[
        o.sample[, 2L] > largest.hole,
        ,
        drop = FALSE
      ]

      out_dn <- compute_dn(o.sample)

      if (out_dn$m == 0L) {
        out <- empty_out()
        out$hole <- largest.hole
        return(out)
      }
    }
  }

  dist.fail <- out_dn$dist.fail
  m <- out_dn$m
  di <- out_dn$di
  ni <- out_dn$ni

  if (any(ni <= 0L)) {
    bad <- which(ni <= 0L)

    stop(
      sprintf(
        paste0(
          "Non-positive risk-set size at failure time(s): %s. ",
          "Corresponding ni: %s."
        ),
        paste(dist.fail[bad], collapse = ", "),
        paste(ni[bad], collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (any(di > ni)) {
    bad <- which(di > ni)

    stop(
      sprintf(
        paste0(
          "Found di > ni at failure time(s): %s. ",
          "Corresponding di: %s; ni: %s."
        ),
        paste(dist.fail[bad], collapse = ", "),
        paste(di[bad], collapse = ", "),
        paste(ni[bad], collapse = ", ")
      ),
      call. = FALSE
    )
  }

  estimation <- cumprod(1 - di / ni)

  cum.prod <- c(1, estimation)
  prob <- cum.prod[seq_len(m)] - cum.prod[seq_len(m) + 1L]

  var.i.u <- rep(NA_real_, m)
  var.i.l <- rep(NA_real_, m)
  varF.u <- rep(NA_real_, m)
  varF.l <- rep(NA_real_, m)

  if (conf.int) {

    if (!is.numeric(alpha) ||
        length(alpha) != 1L ||
        !is.finite(alpha) ||
        alpha <= 0 ||
        alpha >= 1) {
      stop("'alpha' must be a finite number in (0, 1).")
    }

    z.alpha <- qnorm(1 - alpha / 2)

    gterm <- di / (ni * (ni - di))

    gterm[!is.finite(gterm)] <- 0

    cum.sum <- cumsum(gterm)

    if (conf.type == "plain") {

      seS <- estimation * sqrt(cum.sum)

      var.i.l <- pmax(
        estimation - z.alpha * seS,
        0
      )

      var.i.u <- pmin(
        estimation + z.alpha * seS,
        1
      )

    } else {

      var.i.l <- numeric(m)
      var.i.u <- numeric(m)

      ok <- is.finite(estimation) &
        estimation > 0 &
        estimation < 1

      at_one <- is.finite(estimation) & estimation >= 1
      var.i.l[at_one] <- 1
      var.i.u[at_one] <- 1

      at_zero <- is.finite(estimation) & estimation <= 0
      var.i.l[at_zero] <- 0
      var.i.u[at_zero] <- 0

      if (any(ok)) {
        w <- log(-log(estimation[ok]))

        se.w <- sqrt(cum.sum[ok]) /
          abs(log(estimation[ok]))

        var.i.l[ok] <- exp(
          -exp(w + z.alpha * se.w)
        )

        var.i.u[ok] <- exp(
          -exp(w - z.alpha * se.w)
        )
      }
    }

    varF.l <- 1 - var.i.u
    varF.u <- 1 - var.i.l
  }

  out <- list(
    fail.time = dist.fail,
    estimation = estimation,
    prob = prob,
    conf.limits = if (conf.int) {
      cbind(var.i.l, var.i.u)
    } else {
      NULL
    },
    conf.limits.F = if (conf.int) {
      cbind(varF.l, varF.u)
    } else {
      NULL
    },
    hole = largest.hole
  )

  class(out) <- "estimators"
  out
}


