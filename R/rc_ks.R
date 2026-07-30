#' Kolmogorov--Smirnov-type statistic for the k-sample problem under right censoring
#'
#' Computes the Kolmogorov--Smirnov statistic adapted to right censoring.
#'
#' @param p.sample Matrix or data frame with columns (Y, Delta, group).
#' @param weights Vector of weights of the estimator under the null hypothesis.
#' If no vector is provided, weights proportional to each sample size are assigned automatically
#' by the function. Not used when `tests="logrank"`.
#' @param plot.curves Logical, FALSE by default. Otherwise, it represents the k distribution
#' functions to be compared.
#'
#' @return A list with the following components:
#' * `statistic`: numerical value of the Kolmogorov--Smirnov statistics.
#' * `time`: time when the maximum is attained.
#'
#' @details
#' In Kiefer (1959), the \eqn{k}-sample version of the Kolmogorov--Smirnov test is defined as
#' \deqn{D_{KS} = \sup_t \sum_{j=1}^k n_j \left( F_{jn_j}(t) - F_n(t)\right)^2,}
#' where \eqn{F_{jn_j}} is the cumulative distribution function estimator of calculated
#' only with the \eqn{j}-th sample, and \eqn{F_n} is given by
#' \deqn{F_n(t) =\sum_{j=1}^k p_j F_{jn_j}(t),}where \eqn{p_1, \ldots, p_k} correspond to the
#' input weights. This function implements the previous test statistic accounting for
#' right censoring.
#'
#' For the particular case of \eqn{k=2}, the test statistic reduces to
#' \deqn{D_{KS} = (n_1p_2^2 + n_2p_1^2) \left(\sup_t \mid F_{1n_1}(t) - F_{2n_2}(t) \mid\right)^2,}
#' thus the test indeed extends the classical two-sample Kolmogorov--Smirnov test. A two-sample
#' Kolmogorov--Smirnov-type test for right-censored data was proposed in Schumacher (1984).
#'
#' @references
#' Kiefer, J. (1959). K-sample analogues of the Kolmogorov--Smirnov and Cramér--von Mises
#' tests. The Annals of Mathematical Statistics, 30:420-447.
#'
#' Schumacher, M. (1984). Two-Sample Tests of Cramér--von Mises-and Kolmogorov--Smirnov-Type
#' for Randomly Censored Data. International Statistical Review/Revue Internationale de
#' Statistique, 263-281.
#'
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#' library(survival)
#'
#' data <- data.frame(y=stanford2$time,
#' delta=stanford2$status,
#' group=as.numeric(stanford2$age<=35)+as.numeric(stanford2$age<=47))
#' rc.ks(p.sample=data,weights=table(data$group)/nrow(data),plot.curves=TRUE)
#' }
#'
#' @export

rc.ks <- function(p.sample, weights = NULL, plot.curves = FALSE) {

  if (!is.matrix(p.sample) && !is.data.frame(p.sample)) {
    stop("p.sample must be a matrix or data.frame.")
  }

  if (ncol(p.sample) < 3L) {
    stop("p.sample must have at least 3 columns: (Y, delta, group).")
  }

  if (any(p.sample[, 2] != 0 & p.sample[, 2] != 1)) {
    stop("delta must be 0/1.")
  }

  p.sample <- data.frame(
    Y = as.numeric(p.sample[[1]]),
    delta = as.numeric(p.sample[[2]]),
    group = p.sample[[3]]
  )

  nodes <- sort(unique(p.sample[p.sample[, 2] == 1, 1]))
  n <- length(nodes)

  if (n == 0L) {
    return(list(statistic = 0, time = NA_real_))
  }

  g <- factor(p.sample[, 3])
  lev <- levels(g)
  k <- length(lev)
  groups <- tabulate(g)

  if (is.null(weights)) {
    weights <- groups / sum(groups)
  } else { if (!is.numeric(weights) || length(weights) != k) {
    stop(" 'weights' must be a numeric vector with one entry per group.")
  }
  if (any(!is.finite(weights)) || any(weights < 0)) {
    stop(" 'weights' must contain finite non-negative values.")
  }
  if (abs(sum(weights) - 1) > sqrt(.Machine$double.eps)) {
    stop(" 'weights' must sum to 1.")
  }
  }

  E <- matrix(1, nrow = k, ncol = n)

  for (i in seq_len(k)) {
    est <- rc.estimator(
      p.sample[g == lev[i], , drop = FALSE]
    )

    if (length(est$fail.time) > 0L) {
      idx <- findInterval(nodes, est$fail.time, rightmost.closed = TRUE)
      E[i, ] <- c(1, est$estimation)[idx + 1L]
    }
  }

  est.pooled <- as.vector(crossprod(weights, E))

  M <- sweep(E, 2L, est.pooled, FUN = "-")^2
  M <- sweep(M, 1L, groups, FUN = "*")

  colM <- colSums(M)
  ind.max <- which.max(colM)

  if (plot.curves) {
    matplot(nodes, t(E), type = "s", lwd = 2, lty = 1,
            ylim = c(0, 1),
            xlab = "t", ylab = "Survival function estimators",
            main = "")
  }

  list(statistic = colM[ind.max], time = nodes[ind.max])
}
