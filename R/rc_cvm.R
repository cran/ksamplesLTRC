#' Cramér--von Mises-type statistic for the k-sample problem under
#' right censoring
#'
#' Computes the Cramér--von Mises statistic adapted to left truncation and right
#' censoring.
#'
#' @param p.sample Matrix or data frame with columns (U, Y, Delta, group).
#' @param weights Vector of weights of the estimator under the null hypothesis.
#' If no vector is provided, weights proportional to each sample size are assigned automatically
#' by the function. Not used when `tests="logrank"`.
#' @param plot.curves Logical, default to FALSE. If TRUE, the \eqn{k} distributions to be compared are depicted
#'
#' @return Numerical value of the test statistic computed from the \eqn{k} samples.
#'
#' @details
#' In Kiefer (1959), the \eqn{k}-sample version of the Cramér--von Mises test is defined as
#' \deqn{D_{CvM} = \sum_{j=1}^k n_j \int \left( F_{jn_j}(t) - F_n(t)\right)^2 dF_n(t),}
#' where \eqn{F_{jn_j}} is the cumulative distribution function estimator of calculated
#' only with the \eqn{j}-th sample, and \eqn{F_n} is given by
#' \deqn{F_n(t) =\sum_{j=1}^k p_j F_{jn_j}(t),}where \eqn{p_1, \ldots, p_k} correspond to the
#' input weights. This function implements the previous test statistic accounting for
#' right censoring.
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
#' rc.cvm(data,weights=table(data$group)/nrow(data))
#' }
#'
#' @export

rc.cvm <- function(p.sample, weights = NULL, plot.curves = FALSE) {

  if (!is.matrix(p.sample) && !is.data.frame(p.sample)) {
    stop("sample must be a matrix or data.frame.")
  }

  if (ncol(p.sample) < 3L) {
    stop("sample must have at least 3 columns: (Y, delta, group).")
  }

  p.sample <- data.frame(
    Y = as.numeric(p.sample[, 1]),
    delta = as.numeric(p.sample[, 2]),
    group = p.sample[, 3]
  )

  if (any(p.sample[, 2] != 0 & p.sample[, 2] != 1)) {
    stop("delta must be 0/1.")
  }

  groups <- factor(p.sample[, 3])
  lev <- levels(groups)
  k <- length(lev)

  if (k < 2L) {
    stop("Need at least two groups.")
  }

  N <- tabulate(groups, nbins = k)
  n_total <- sum(N)

  if (is.null(weights)) {
    weights <- N / n_total
  } else {
    if (!is.numeric(weights) || length(weights) != k) {
      stop("'weights' must be a numeric vector with one entry per group.")
    }

    if (any(!is.finite(weights)) || any(weights < 0)) {
      stop("'weights' must contain finite non-negative values.")
    }

    if (!isTRUE(all.equal(sum(weights), 1))) {
      stop("'weights' must sum to 1.")
    }
  }

  nodes <- sort(unique(p.sample[p.sample[, 2] == 1, 1]))
  L <- length(nodes)

  if (L == 0L) {
    stop("No failures observed.")
  }

  F_mat <- matrix(0, nrow = k, ncol = L)
  dF_mat <- matrix(0, nrow = k, ncol = L)

  for (i in seq_len(k)) {
    est <- rc.estimator(
      p.sample[groups == lev[i], 1:2, drop = FALSE]
    )

    if (length(est$fail.time) > 0L) {
      idx <- findInterval(nodes, est$fail.time)
      F_mat[i, ] <- 1 - c(1, est$estimation)[idx + 1L]

      jump.idx <- match(est$fail.time, nodes)
      dF_mat[i, jump.idx] <- est$prob
    }
  }

  F_pooled <- as.vector(crossprod(weights, F_mat))
  dF_pooled <- as.vector(crossprod(weights, dF_mat))

  diff2 <- sweep(F_mat, 2L, F_pooled, FUN = "-")^2

  M <- rowSums(
    diff2 *
      matrix(
        dF_pooled,
        nrow = k,
        ncol = L,
        byrow = TRUE
      )
  )

  statistic <- sum(N * M)

  if (plot.curves) {
    matplot(
      nodes,
      t(1 - F_mat),
      type = "s",
      lty = 1,
      lwd = 2,
      ylim = c(0, 1),
      xlab = "t",
      ylab = "Survival function estimators"
    )
  }

  list(statistic = statistic)
}

