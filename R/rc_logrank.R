#' Weighted log-rank tests under left truncation
#'
#' Computes weighted log-rank-type tests (Fleming–Harrington family)
#' for k groups with left truncation and right censoring.
#'
#' @param p.sample Matrix or data frame with columns (U, Y, Delta, group).
#' @param p,q Nonnegative weight parameters. If one has length 1 and the other is a vector, the scalar is recycled.
#' @param plot.curves Logical, default to FALSE. If TRUE, the \eqn{k} distributions to be compared are depicted
#'
#' @return A matrix with columns \code{statistic} and \code{p.value}.
#'
#' @details
#' The test statistic is based on the quantity
#' \deqn{L_j = \sum_{i=1}^m \omega_i \left(d_{ji} - d_i \frac{r_{ji}}{r_i}\right),}
#' for \eqn{j \in \{1, \ldots, k\}}. Denote by \eqn{y_1, \ldots, y_m} the distinct
#' observed event times from the \eqn{k} samples. For a time \eqn{y_i} and a group \eqn{j},
#' \eqn{d_{ji}} denotes the total count of individuals that experience the event,
#' and \eqn{r_{ji}} is the number of individuals at risk. In addition,
#' \eqn{d_i = \sum_{j=1}^k d_{ji}} and  \eqn{r_i = \sum_{j=1}^k r_{ji}}, which are
#' the total number of events observed and the total number of individuals at risk
#' at time \eqn{y_i}, respectively.
#'
#' Finally, the weights \eqn{\omega_1, \ldots, \omega_m} are taken as
#' \deqn{w_i = \hat{S}(y_i)^p (1-\hat{S}(y_i))^q,}
#' with \eqn{\hat{S}(x_i)} the pooled survival just before the event time. The choice
#' \eqn{p=0=q} results in the classical version of the log-rank test.
#'
#' The test statistic is given by
#' \deqn{(L_1, \ldots, L_{k-1}) \hat{\Sigma}^{-1}(L_1, \ldots, L_{k-1})^\top,}
#' where \eqn{\hat{\Sigma}} is the matrix with the variance-covariance estimators. Each
#' term in the matrix includes the correction factor \eqn{(n_j-d_j)/(n_j-1)} for tied
#' event times. The test statistic follows an \eqn{\chi^2_{k-1}} distribution under
#' the null hypothesis of equal target distributions in the \eqn{k} populations.
#'
#' @references
#' Klein, J.P., and Moeschberger, M. L. (2003). Survival Analysis: Techniques for Censored
#' and Truncated Data. Springer, New York.
#'
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#' library(survival)
#' data <- data.frame(y=stanford2$time,
#' delta=stanford2$status,
#' group=as.numeric(stanford2$age<=35)+as.numeric(stanford2$age<=47))
#' rc.logrank(data, p = c(0, 1, 2), q = 0)
#' }
#'
#' @export

rc.logrank <- function(p.sample, p = 0, q = 0, plot.curves = FALSE) {

  p <- as.numeric(p)
  q <- as.numeric(q)

  if (length(p) != length(q)) {
    if (length(p) == 1L) {
      p <- rep(p, length(q))
    } else if (length(q) == 1L) {
      q <- rep(q, length(p))
    } else {
      stop("check the weights: `p` and `q` must have the same length, or one of them must have length 1")
    }
  }

  # ---- Basic input checks ----
  if (!is.matrix(p.sample) && !is.data.frame(p.sample)) {
    stop("'p.sample' must be a matrix or data.frame.")
  }
  p.sample <- data.frame(
    Y = as.numeric(p.sample[[1]]),
    delta = as.numeric(p.sample[[2]]),
    group = p.sample[[3]]
  )

  if (ncol(p.sample) < 3L) {
    stop("p.sample must have at least 3 columns: (time, delta, group).")
  }

  if (any(!is.finite(p.sample[, 1]))) {
    stop("Observed times must be finite.")
  }

  if (any(p.sample[, 1] < 0)) {
    stop("Observed times must be non-negative.")
  }

  if (any(p.sample[, 2] != 0 & p.sample[, 2] != 1)) {
    stop("delta must be 0/1.")
  }

  groups <- sort(unique(p.sample[, 3]))
  K <- length(groups)

  if (K < 2L) {
    stop("At least two groups are required.")
  }

  times <- sort(unique(p.sample[p.sample[, 2] == 1, 1]))
  D <- length(times)

  if (D == 0L) {
    stop("There are no observed failures in the sample.")
  }

  d <- n <- matrix(numeric(K * D), nrow = D, ncol = K)

  for (i in 1:K) {
    samp <- p.sample[p.sample[, 3] == groups[i], , drop = FALSE]

    idx <- match(samp[samp[, 2] == 1, 1], times)
    d[, i] <- tabulate(idx, nbins = D)

    # risk set under right censoring:
    # subjects with observed time >= t are at risk just before t
    exit <- findInterval(times - 1e-6, sort(samp[, 1]))
    n[, i] <- nrow(samp) - exit
  }

  d.pooled <- rowSums(d)
  n.pooled <- rowSums(n)

  Z <- numeric(K - 1)
  varcov <- matrix(numeric((K - 1)^2), ncol = K - 1)

  est.pooled <- rc.estimator(p.sample[, 1:2, drop = FALSE], conf.int = FALSE)
  ss <- length(p)

  output <- matrix(numeric(2 * ss), nrow = ss, ncol = 2)
  colnames(output) <- c("statistic", "p.value")
  rownames(output) <- 1:ss

  phi <- c(1, est.pooled$estimation[1:(D - 1)])

  corr <- rep(0, D)
  ok <- n.pooled > 1
  corr[ok] <- (n.pooled[ok] - d.pooled[ok]) / (n.pooled[ok] - 1)

  for (l in 1:ss) {

    Z[] <- 0
    varcov[,] <- 0

    weights <- phi^p[l] * (1 - phi)^q[l]
    rownames(output)[l] <- paste("p=", p[l], "q=", q[l], ";")

    aa <- d.pooled / n.pooled
    bb <- d.pooled / n.pooled^2

    for (i in 1:(K - 1)) {
      Z[i] <- sum(weights * (d[, i] - aa * n[, i]))

      varcov[i, i] <- sum(weights^2 * n[, i] * aa * (1 - n[, i] / n.pooled) * corr)

      if (i < (K - 1)) {
        for (j in (i + 1):(K - 1)) {
          varcov[i, j] <- -sum(weights^2 * bb * n[, i] * n[, j] * corr)
          varcov[j, i] <- varcov[i, j]
        }
      }
    }

    output[l, 1] <- Z %*% solve(varcov) %*% Z
    output[l, 2] <- 1 - pchisq(output[l, 1], df = (K - 1))
  }

  if (plot.curves) {
    xmin <- min(p.sample[, 1])
    xmax <- max(p.sample[, 1])

    plot(
      rc.estimator(p.sample[p.sample[, 3] == groups[1], 1:2, drop = FALSE],
                   conf.int = FALSE),
      col = 2,
      xlim = c(xmin, xmax),
      xlab = "t",
      ylab = "Survival function estimators"
    )

    for (i in 2:K) {
      lines(
        rc.estimator(p.sample[p.sample[, 3] == groups[i], 1:2, drop = FALSE],
                     conf.int = FALSE),
        col = i + 1
      )
    }
  }

  return(output)
}
