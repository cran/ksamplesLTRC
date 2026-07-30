#' Kolmogorov--Smirnov-type test for the \eqn{k}-sample problem under left truncation
#' and right censoring
#'
#' Computes the Kolmogorov--Smirnov test adapted to left truncation and right
#' censoring and the corresponding \eqn{p}-value.
#'
#' @param p.sample Matrix or data frame with columns (U, Y, \eqn{\Delta}, group).
#' @param weights Vector of weights of the estimator under the null hypothesis.
#' If no vector is provided, weights proportional to each sample size are assigned automatically
#' by the function. Not used when `tests="logrank"`.
#' @param B Number of bootstrap replications to approximate the \eqn{p}-value.
#' @param plot.curves Logical, default to FALSE. If TRUE, the \eqn{k} curves to be
#' compared are depicted.
#' @param holes How to handle holes if present in a sample. One of `holes` (default),
#' `conditional` or `add`. See [ltrc.estimator()] for details.
#' @param keep.boot Logical, default to FALSE. If TRUE, the function returns the statistics
#' evaluated in the bootstrap resamples.
#'
#' @return A list with the following components:
#' * `statistic`: numerical value of the Kolmogorov--Smirnov statistic.
#' * `p.value`: approximated p-value of the tests.
#' * `curves`: A list with the \eqn{k} estimators of the distribution functions to
#' be compared.
#' * `groups`: Group indicator
#' * `stat.boot`: the statistics computed from the bootstrap resamples. Included in the output
#' only if `keep.boot = TRUE`.
#'
#' @details
#' In Kiefer (1959), the \eqn{k}-sample version of the Kolmogorov--Smirnov test is defined as
#' \deqn{D_{KS} = \sup_t \sum_{j=1}^k n_j \left( F_{jn_j}(t) - F_n(t)\right)^2,}
#' where \eqn{F_{jn_j}} is the cumulative distribution function estimator of calculated
#' only with the \eqn{j}-th sample, and \eqn{F_n} is given by
#' \deqn{F_n(t) =\sum_{j=1}^k p_j F_{jn_j}(t),}where \eqn{p_1, \ldots, p_k} correspond to the
#' input weights. This function implements the previous test statistic accounting for left
#' truncation and right censoring, as studied in Lago et al. (2026).
#'
#' The \eqn{p}-value is approximated by means of the obvious bootstrap accounting for
#' left truncation and right censoring, studied in Bilker and Wang (1997). The exact
#' algorithm for the \eqn{k}-sample problem can be consulted in Lago et al. (2026).
#'
#' @references
#' Bilker, W.B., and Wang, M.C. (1997). Bootstrapping left truncated
#' and right censored data. Communications in Statistics - Simulation and Computation, 26:141-171.
#'
#' Lago, A., Pardo-Fernández, J.C., and de Uña-Álvarez, J. (2026). Kolmogorov--Smirnov
#' and Cramér--von Mises tests for the k-sample problem for left-truncated and
#' right-censored data. Lifetime Data Analysis: 32(2).
#'
#' Kiefer, J. (1959). K-sample analogues of the Kolmogorov--Smirnov and Cramér--von Mises
#' tests. The Annals of Mathematical Statistics, 30:420-447.
#'
#' @examples
#' if (requireNamespace("KMsurv", quietly = TRUE)) {
#' data("channing", package = "KMsurv")
#' data.channing <- data.frame(channing$ageentry,channing$age,
#' channing$death,channing$gender)
#' aa <- ltrc.pv.ks(data.channing,weights=rep(1/2,2),
#' holes = 'conditional', plot.curves = FALSE, B = 200)
#' plot(aa)
#' aa$pvalue
#' }
#'
#' @export

ltrc.pv.ks <- function(p.sample, weights  = NULL, B = 500, plot.curves = FALSE,
                       holes = c('holes','add','conditional'),
                       keep.boot = FALSE) {

  holes <- match.arg(holes)

  if (!is.matrix(p.sample) && !is.data.frame(p.sample)) {
    stop("'p.sample' must be a matrix or data.frame.")
  }
  p.sample <- as.data.frame(p.sample)

  p.sample <- data.frame(
    U = as.numeric(p.sample[[1]]),
    Y = as.numeric(p.sample[[2]]),
    delta = as.numeric(p.sample[[3]]),
    group = p.sample[[4]]
  )
  if (ncol(p.sample) < 4L) {
    stop("p.sample must have at least 4 columns: (U, Y, delta, group).")
  }
  if (any(p.sample[, 1] > p.sample[, 2])) {
    stop("All observations must satisfy U <= Y.")
  }
  if (any(p.sample[, 3] != 0 & p.sample[, 3] != 1)) {
    stop("delta must be 0/1.")
  }

  ind.group <- sort(unique(p.sample[, 4]))
  K <- length(ind.group)

  if (is.null(weights)) {
    N <- tabulate(factor(p.sample[,4]))
    weights <- N / sum(N)
  } else {
    if (!is.numeric(weights) || length(weights) != K) {
      stop("'weights' must be a numeric vector with one entry per group.",
           call. = FALSE)
    }
    if (any(!is.finite(weights)) || any(weights < 0)) {
      stop("'weights' must contain finite non-negative values.", call. = FALSE)
    }
    if (abs(sum(weights) - 1) > sqrt(.Machine$double.eps)) {
      stop("'weights' must sum to 1.", call. = FALSE)
    }
  }

  detect_hole <- function(x) {
    out <- rep(NA_real_, K)
    for (j in seq_len(K)) {
      sj  <- x[x[, 4] == ind.group[j], 1:3, drop = FALSE]
      est <- ltrc.estimator(o.sample = sj, holes = "holes")
      out[j] <- if (is.null(est$hole)) NA_real_ else est$hole
    }
    out
  }
  hole.ind <- detect_hole(p.sample)
  if (holes == "conditional" && any(!is.na(hole.ind))) {
    largest.hole <- max(hole.ind, na.rm = TRUE)
    p.sample <- p.sample[p.sample[, 2] > largest.hole, , drop = FALSE]
  }
  holes.internal <- if (holes == "add") "add" else "holes"

  stat.0 <- ltrc.ks(p.sample, weights = weights, plot.curves = plot.curves,
                    holes = holes.internal)$statistic

  F0 <- vector("list", K)
  QQ <- vector("list", K)
  GG <- vector("list", K)
  NN <- numeric(K)

  for (j in seq_len(K)) {
    sj <- p.sample[p.sample[, 4] == ind.group[j], , drop = FALSE]
    NN[j] <- nrow(sj)

    F0[[j]] <- ltrc.estimator(sj, conf.int = FALSE, holes = holes.internal)
    GG[[j]] <- ltrc.truncation.estimator(sj, holes = holes)
    QQ[[j]] <- rc.estimator(
      cbind(sj[, 2] - sj[, 1], 1 - sj[, 3])
    )

    if (length(F0[[j]]$fail.time) == 0L) {
      stop(paste("Group", j, "has no observed failure times."))
    }
    if (length(GG[[j]]$trunc.time) == 0L) {
      stop(paste("Group", j, "has no truncation times."))
    }
    if (length(QQ[[j]]$fail.time) == 0L) {
      QQ[[j]]$fail.time <- Inf
      QQ[[j]]$prob <- 1
      warning(paste("Group", j, "has no censoring times."))
    }else if(length(QQ[[j]]$fail.time) == 1L){
      QQ[[j]]$fail.time <- rep(QQ[[j]]$fail.time , 2 )
      QQ[[j]]$prob <- rep(1/2,2)
    }
  }

  stat.boot <- numeric(B)

  for (b in seq_len(B)) {
    resample <- vector("list", K)

    for (j in seq_len(K)) {
      n.j <- NN[j]

      trunc.boot <- sample(GG[[j]]$trunc.time, n.j, TRUE, GG[[j]]$prob)

      ind.evtimes <- sample(seq_len(K), n.j, TRUE, weights)
      evtimes.boot <- numeric(n.j)
      for (k in seq_len(K)) {
        idx <- which(ind.evtimes == k)
        if (length(idx) > 0L) {
          evtimes.boot[idx] <- sample(F0[[k]]$fail.time, length(idx), TRUE, F0[[k]]$prob)
        }
      }

      ind.bad <- trunc.boot > evtimes.boot
      while (any(ind.bad)) {
        m <- sum(ind.bad)
        trunc.boot[ind.bad] <- sample(GG[[j]]$trunc.time, m, TRUE, GG[[j]]$prob)
        ind.evtimes.new <- sample(seq_len(K), m, TRUE, weights)
        evtimes.new <- numeric(m)
        for (k in seq_len(K)) {
          idx <- which(ind.evtimes.new == k)
          if (length(idx) > 0L) {
            evtimes.new[idx] <- sample(F0[[k]]$fail.time, length(idx), TRUE, F0[[k]]$prob)
          }
        }
        evtimes.boot[ind.bad] <- evtimes.new
        ind.bad <- trunc.boot > evtimes.boot
      }

      rc.boot <- sample(QQ[[j]]$fail.time, n.j, TRUE, QQ[[j]]$prob)

      resample[[j]] <- data.frame(
        u = trunc.boot,
        y = pmin(evtimes.boot, trunc.boot + rc.boot),
        delta = as.numeric(evtimes.boot <= trunc.boot + rc.boot),
        group = rep(ind.group[j], n.j)
      )
    }

    resample.df <- do.call(rbind, resample)
    row.names(resample.df) <- NULL

    stat.boot[b] <- ltrc.ks(resample.df, weights = weights, plot.curves = FALSE,
                            holes = holes.internal)$statistic
  }

  pvalue <- mean(c(stat.boot, stat.0) >= stat.0)

  out <- list(statistic = stat.0,
              pvalue = pvalue,
              curves = F0,
              groups = ind.group
              )

  if(keep.boot){
    out$stat.boot <- stat.boot
  }

  class(out) <- "ksamples_test"

  return(out)
}
