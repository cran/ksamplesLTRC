#' Survival function estimator under left truncation and right
#' censoring
#'
#' Computes the survival function estimator accounting for late entry
#' times and censored observations. It also returns pointwise confidence bands.
#'
#' @param o.sample Matrix or data frame with columns (U, Y, Delta, group).
#' @param conf.int Logical, set as default to FALSE. If set to TRUE, the function
#' returns a pointwise confidence band.
#' @param alpha Pointwise confidence bands with \eqn{1-\alpha} confidence level
#' @param conf.type Type of pointwise confidence band. Set as default to \eqn{plain}.
#' @param holes One of 'add','conditional', or 'holes' (default) (see details).
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
#' The confidence intervals are based on the asymptotic Gaussianity of the process
#' \eqn{\sqrt{n} (S_n(t) - S(t))}. Different transformations are also possible.
#'
#' The function also includes a way to deal with possible holes. A hole is an event time with associated
#' unitary risk set. If a hole is present in a sample, the estimator degenerates and takes the
#' value 1 from that time on. Different techniques are available in the literature for the treatment
#' of holes. For instance, in Stute and Wang (2008), authors suggest to increase the risk set by 1 (holes='add').
#' Other references, like Klein and Moeschberger (2003), suggest conditioning on the largest hole different
#' from the largest event time (holes='conditional'). The possibility of computing the degenerate statistic
#' is also included (holes='holes').
#'
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
#' observed.sample <- ltrc.sim(n=1000,rtrunc=rnorm,rtarget=rnorm,rcens=rexp,
#' trunc.args=list(4,1),target.args = list(4,1),cens.args = list(1))
#' est <- ltrc.estimator(observed.sample,conf.int=TRUE,alpha=0.05,conf.type='loglog')
#' str(est)
#' plot(est,col='firebrick',xlab='t', ylab='Estimated survival probability')
#'
#' if (requireNamespace("KMsurv", quietly = TRUE)) {
#' data("channing", package = "KMsurv")
#' o.sample <- channing[channing$gender==1,c(3,4,2)] # subgroup of men
#' head(o.sample)
#' plot( ltrc.estimator(o.sample,holes='holes'),col=2,xlab='t',
#' ylab='Estimated survival probability' ,
#' main = 'Estimated survival probability for men') # no correction
#' lines( ltrc.estimator(o.sample,holes='add'),col=3 ) # +1 to risk set
#' lines( ltrc.estimator(o.sample,holes='conditional'),col=4 )
#' legend('topright',legend=c('holes="holes" ','holes="add" ',
#' 'holes="conditional"'),col=2:4,lwd=rep(1,3))
#' }
#'
#' @export

ltrc.estimator <- function(
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

  if (ncol(o.sample) < 3L) {
    stop("'o.sample' must have at least 3 columns: (U, Y, delta).")
  }

  o.sample <- data.frame(
    U = as.numeric(o.sample[, 1L]),
    Y = as.numeric(o.sample[, 2L]),
    delta = as.numeric(o.sample[, 3L])
  )

  if (anyNA(o.sample)) {
    stop("'o.sample' must not contain missing values.")
  }


  if (any(!o.sample$delta %in% c(0, 1))) {
    stop("'delta' must contain only 0 and 1.")
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
    fail_ind <- x[, 3] == 1
    dist.fail <- sort(unique(x[fail_ind, 2]))
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
      function(t) sum(x[, 2] == t & x[, 3] == 1),
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

  prob <- cum.prod[seq_len(m)] -
    cum.prod[seq_len(m) + 1L]

  var.i.u <- rep(NA_real_, m)
  var.i.l <- rep(NA_real_, m)
  varF.u <- rep(NA_real_, m)
  varF.l <- rep(NA_real_, m)

  if (conf.int) {

    if (
      !is.numeric(alpha) ||
      length(alpha) != 1L ||
      !is.finite(alpha) ||
      alpha <= 0 ||
      alpha >= 1
    ) {
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










