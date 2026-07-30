#' Kolmogorov--Smirnov-type statistic for the \eqn{k}-sample problem under left truncation
#'
#' Computes the Kolmogorov--Smirnov statistic adapted to left truncation.
#'
#' @param p.sample Matrix or data frame with columns (U, X, group).
#' @param holes What to do when a hole is present in a sample. See [lt.estimator()]
#' for more details.
#' @param plot.curves Logical, default to FALSE. If true, the \eqn{k} distributions to
#' be compared are depicted.
#' @param weights Vector of weights of the estimator under the null hypothesis.
#' If no vector is provided, weights proportional to each sample size are assigned automatically
#' by the function. Not used when `tests="logrank"`.
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
#' truncation. The two-sample version of this statistic was studied in Lago et al. (2025).
#' Note that, when \eqn{k=2}, one has
#' \deqn{D_{KS} = (n_1p_2^2 + n_2p_1^2) \left( \sup_t \mid F_{1n_1}(t) - F_{2n_2}(t) \mid \right)^2.}
#'
#' @references
#' Kiefer, J. (1959). K-sample analogues of the Kolmogorov--Smirnov and Cramér--von Mises
#' tests. The Annals of Mathematical Statistics, 30:420-447.
#'
#' Lago, A., de Uña-Álvarez, J., & Pardo-Fernández, J. C. (2025). A Kolmogorov--Smirnov-type
#' test for the two-sample problem with left-truncated data. Test, 34(1), 69-90.
#'
#' @examples
#' if (requireNamespace("mvna", quietly = TRUE)) {
#' data("abortion", package = "mvna")
#' pooled.sample <- abortion[abortion$cause==3, 2:4]
#' head(pooled.sample)
#' aa <- lt.ks(p.sample=pooled.sample,weights=rep(1/2,2),plot.curves=TRUE)
#' aa
#' }
#'
#'
#' @export

lt.ks <- function(p.sample, weights = NULL, plot.curves = FALSE,
                  holes = c("holes", "add", "conditional")) {

  holes <- match.arg(holes)

  if (!is.matrix(p.sample) && !is.data.frame(p.sample)) {
    stop("'p.sample' must be a matrix or data.frame.")
  }
  p.sample <- as.data.frame(p.sample)
  p.sample <- data.frame(
    U = as.numeric(p.sample[[1]]),
    Y = as.numeric(p.sample[[2]]),
    group = p.sample[[3]]
  )

  if (ncol(p.sample) < 3L) {
    stop("'p.sample' must have at least 3 columns: (U, X, group).")
  }

  if (any(p.sample[, 1] > p.sample[, 2])) {
    stop("All observations must satisfy U <= X.")
  }

  nodes <- sort(unique(p.sample[, 2]))
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
    group_sample <- p.sample[g == lev[i], 1:2, drop = FALSE]

    est <- lt.estimator(group_sample, holes = holes)

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

