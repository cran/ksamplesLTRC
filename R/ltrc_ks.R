#' Kolmogorov--Smirnov-type statistic for the \eqn{k}-sample problem under left truncation
#' and right censoring
#'
#' Computes the Kolmogorov--Smirnov statistic adapted to left truncation and right
#' censoring.
#'
#' @param p.sample Matrix or data frame with columns (U, Y, Delta, group).
#' @param weights Vector of weights of the estimator under the null hypothesis.
#' If no vector is provided, weights proportional to each sample size are assigned automatically
#' by the function. Not used when `tests="logrank"`.
#' @param plot.curves Logical, default to FALSE. If TRUE, the \eqn{k} curves to be
#' compared are depicted.
#' @param holes How to handle holes if present in a sample. One of `holes` (default),
#' `conditional` or `add`. See [ltrc.estimator()] for details.
#'
#' @return A list with the following components:
#' * `statistic`: numerical value of the Kolmogorov--Smirnov statistic.
#' * `time`: time when the maximum is attained.
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
#' @references
#' Lago, A., Pardo-Fernández, J.C., and de Uña-Álvarez, J. (2026) Kolmogorov--Smirnov
#' and Cramér--von Mises tests for the k-sample problem for left-truncated and
#' right-censored data. Lifetime Data Analysis: 32(2).
#'
#' @examples
#' if (requireNamespace("KMsurv", quietly = TRUE)) {
#' data("channing", package = "KMsurv")
#'
#' data.channing <- data.frame(channing$ageentry, channing$age,
#' channing$death, channing$gender)
#' ltrc.ks(data.channing,weights=rep(1/2,2),holes='conditional')
#' }
#'
#' @export

ltrc.ks <- function(p.sample, weights  = NULL, plot.curves = FALSE,
                    holes = c("holes", "add", "conditional")) {
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
    stop("'p.sample' must have at least 4 columns: (U, Y, delta, group).")
  }
  if (any(p.sample[, 1] > p.sample[, 2])) {
    stop("All observations must satisfy U <= Y.")
  }
  if (any(!p.sample[, 3] %in% c(0, 1))) {
    stop("'delta' must be 0/1.")
  }

  if (holes == "conditional") {
    ind.group <- sort(unique(p.sample[, 4]))
    K <- length(ind.group)
    hole.ind <- rep(NA_real_, K)
    for (j in seq_len(K)) {
      sj  <- p.sample[p.sample[, 4] == ind.group[j], 1:3, drop = FALSE]
      est <- ltrc.estimator(o.sample = sj, holes = "holes")
      hole.ind[j] <- if (is.null(est$hole)) NA_real_ else est$hole
    }
    if (any(!is.na(hole.ind))) {
      largest.hole <- max(hole.ind, na.rm = TRUE)
      p.sample <- p.sample[p.sample[, 2] > largest.hole, , drop = FALSE]
    }
    holes <- "holes"
  }

  nodes <- sort(unique(p.sample[p.sample[, 3] == 1, 2]))
  n <- length(nodes)
  if (n == 0L) {
    return(list(statistic = 0, time = NA_real_))
  }

  g <- factor(p.sample[, 4])
  lev <- levels(g)
  k <- length(lev)
  groups <- tabulate(g)

  if (is.null(weights)) {
    weights <- groups / sum(groups)
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

  E <- matrix(1, nrow = k, ncol = n)

  for (i in seq_len(k)) {
    group_sample <- p.sample[g == lev[i], , drop = FALSE]
    est <- ltrc.estimator(group_sample, holes = holes)
    if (length(est$fail.time) > 0L) {
      idx <- findInterval(nodes, est$fail.time, rightmost.closed = TRUE)
      E[i, ] <- c(1, est$estimation)[idx + 1L]
    }
  }

  est.pooled <- as.vector(crossprod(weights, E))
  M <- sweep(E, 2L, est.pooled, FUN = "-")^2
  M <- M * groups
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





