#' Cramér--von Mises-type statistic for the \eqn{k}-sample problem under left truncation
#'
#' Computes the Cramér--von Mises statistic adapted to left truncation.
#'
#' @param p.sample Matrix or data frame with columns (U, X, group).
#' @param holes Way to handle holes, if present. See [lt.estimator()]
#' for more details.
#' @param plot.curves Logical, default to FALSE. If true, the \eqn{k} distributions to
#' be compared are depicted.
#' @param weights Vector of weights of the estimator under the null hypothesis.
#' If no vector is provided, weights proportional to each sample size are assigned automatically
#' by the function. Not used when `tests="logrank"`.
#'
#' @return Numerical value of the test statistic computed from the \eqn{k} samples.
#'
#' @details
#' In Kiefer (1959), the \eqn{k}-sample version of the Cramér--von Mises test is defined as
#' \deqn{D_{CvM} = \sum_{j=1}^k n_j \int \left( F_{jn_j}(t) - F_n(t)\right)^2 dF_n(t),}
#' where \eqn{F_{jn_j}} is the cumulative distribution function estimator of calculated
#' only with the \eqn{j}-th sample, and \eqn{F_n} is given by
#' \deqn{F_n(t) =\sum_{j=1}^k p_j F_{jn_j}(t),}where \eqn{p_1, \ldots, p_k} correspond to the
#' input `weights`. This function implements the previous test statistic accounting for left
#' truncation.
#'
#' @references
#' Kiefer, J. (1959). K-sample analogues of the Kolmogorov--Smirnov and Cramér--von Mises
#' tests. The Annals of Mathematical Statistics, 30:420-447.
#'
#' @examples
#' if (requireNamespace("mvna", quietly = TRUE)) {
#' data("abortion", package = "mvna")
#' pooled.sample <- abortion[abortion$cause==3, 2:4]
#' head(pooled.sample)
#' aa <- lt.cvm(p.sample=pooled.sample,weights=rep(1/2,2),plot.curves=TRUE)
#' aa
#' }
#'
#' @export

lt.cvm <- function(p.sample, weights = NULL,
                   holes = c("holes", "add", "conditional"),
                   plot.curves = FALSE) {

  holes <- match.arg(holes)

  if (!is.matrix(p.sample) && !is.data.frame(p.sample)) {
    stop("'p.sample' must be a matrix or data.frame.", call. = FALSE)
  }
  p.sample <- as.data.frame(p.sample)
  p.sample <- data.frame(
    U = as.numeric(p.sample[[1]]),
    Y = as.numeric(p.sample[[2]]),
    group = p.sample[[3]]
  )

  if (ncol(p.sample) < 3L) {
    stop("'p.sample' must have at least 3 columns: (U, X, group).",
         call. = FALSE)
  }

  if (any(p.sample[, 1] > p.sample[, 2])) {
    stop("All observations must satisfy U <= X.", call. = FALSE)
  }

  g <- factor(p.sample[, 3])
  groups <- levels(g)
  k <- length(groups)
  N <- tabulate(g)

  if (k < 2L) {
    stop("Need at least two groups.", call. = FALSE)
  }

  if (is.null(weights)) {
    weights <- N / sum(N)
  } else {
    if (!is.numeric(weights) || length(weights) != k) {
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

  est.list <- vector("list", k)
  support.len <- integer(k)

  for (i in seq_len(k)) {
    samp.i <- p.sample[g == groups[i], 1:2, drop = FALSE]
    est.i <- lt.estimator(samp.i, conf.int = FALSE, holes = holes)
    est.list[[i]] <- est.i
    support.len[i] <- length(est.i$fail.time)
  }

  if (sum(support.len) == 0L) {
    stop("No observed target times.", call. = FALSE)
  }

  idx.end <- cumsum(support.len)
  idx.start <- c(1L, head(idx.end, -1L) + 1L)
  n <- sum(support.len)

  nodes <- numeric(n)
  Omega <- numeric(n)
  PP <- numeric(n)

  for (i in seq_len(k)) {
    if (support.len[i] > 0L) {
      rng <- idx.start[i]:idx.end[i]
      nodes[rng] <- est.list[[i]]$fail.time
      Omega[rng] <- est.list[[i]]$prob
      PP[rng] <- weights[i]
    }
  }

  E <- matrix(1, nrow = k, ncol = n)

  for (i in seq_len(k)) {
    est.i <- est.list[[i]]

    if (length(est.i$fail.time) > 0L) {
      idx <- findInterval(nodes, est.i$fail.time, rightmost.closed = TRUE)
      E[i, ] <- c(1, est.i$estimation)[idx + 1L]
    }
  }

  est.pooled <- as.vector(crossprod(weights, E))

  weight.measure <- PP * Omega

  M <- numeric(k)
  for (i in seq_len(k)) {
    M[i] <- sum(weight.measure * (E[i, ] - est.pooled)^2)
  }

  statistic <- sum(N * M)

  if (plot.curves) {
    xlim.range <- range(nodes, finite = TRUE)

    plot(est.list[[1]], lwd = 2, col = 2,
         xlab = "t", ylab = "Survival function estimators",
         main = "", xlim = xlim.range)

    if (k >= 2L) {
      for (i in 2:k) {
        lines(est.list[[i]], lwd = 2, col = i + 1)
      }
    }
  }

  list(statistic = statistic)
}


