#' Weighted log-rank tests under left truncation
#'
#' Computes weighted log-rank-type tests (Fleming–Harrington family)
#' for \eqn{k} groups with left truncation.
#'
#' @param p.sample Matrix or data frame with columns (U, X, group).
#' @param plot.curves Logical, default to FALSE. If true, the the \eqn{k} distribution functions
#' to be compared are depicted.
#' @param holes How the function deals with holes when estimating the survival function
#' with the pooled sample. One of 'add' (default), 'conditional', and 'holes'. If no holes
#' are present in the sample, this argument will not be used.
#' @param p,q Nonnegative tuning parameters. If one has length 1 and the other
#'   is a vector, the scalar is recycled.
#'
#' @return A matrix with columns \code{statistic}, detailing the choices of `p` and
#' `q`, in addition to the corresponding \code{p.value}.
#'
#' @details
#' The test statistic is based on the quantity
#' \deqn{L_j = \sum_{i=1}^m \omega_i \left(d_{ji} - d_i \frac{r_{ji}}{r_i}\right),}
#' for \eqn{j \in \{1, \ldots, k\}}. Denote by \eqn{x_1, \ldots, x_m} the distinct
#' observed event times from the \eqn{k} samples. For a time \eqn{x_i} and a group \eqn{j},
#' \eqn{d_{ji}} denotes the total count of individuals that experience the event,
#' and \eqn{r_{ji}} is the number of individuals at risk. In addition,
#' \eqn{d_i = \sum_{j=1}^k d_{ji}} and  \eqn{r_i = \sum_{j=1}^k r_{ji}}, which are
#' the total number of events observed and the total number of individuals at risk
#' at time \eqn{x_i}, respectively.
#'
#' Finally, the weights \eqn{\omega_1, \ldots, \omega_m} are taken as
#' \deqn{w_i = \hat{S}(x_i)^p (1-\hat{S}(x_i))^q,}
#' with \eqn{\hat{S}(x_i)} the pooled survival just before the event time. The choice
#' \eqn{p=0=q} results in the classical version of the log-rank test.
#'
#' The test statistic is given by
#' \deqn{(L_1, \ldots, L_{k-1}) \hat{\Sigma}^{-1}(L_1, \ldots, L_{k-1})^\top,}
#' where \eqn{\hat{\Sigma}} is the matrix with the variance-covariance estimators. Each
#' term in the matrix includes the correction factor \eqn{(n_j-d_j)/(n_j-1)} for tied
#' event times. The test statistic follows an asymptotic \eqn{\chi^2_{k-1}} distribution under
#' the null hypothesis of equal target distributions in the \eqn{k} populations.
#'
#' @references
#' Klein, J.P., and Moeschberger, M. L. (2003). Survival Analysis: Techniques for Censored
#' and Truncated Data. Springer, New York.
#'
#' @examples
#' if (requireNamespace("mvna", quietly = TRUE)) {
#' data("abortion", package = "mvna")
#' pooled.sample <- abortion[abortion$cause==3, 2:4]
#' head(pooled.sample)
#' lt.logrank(p.sample=pooled.sample,p=c(0,1,0.5,5),q=c(0,0,0.5,1),plot.curves=TRUE)
#' }
#'
#' @export

lt.logrank <- function(p.sample, p = 0, q = 0, plot.curves = FALSE,
                        holes = c("holes", "conditional")) {

  holes <- match.arg(holes)
  p <- as.numeric(p)
  q <- as.numeric(q)

  p.sample <- as.data.frame(p.sample)
  p.sample <- data.frame(
    U = as.numeric(p.sample[[1]]),
    Y = as.numeric(p.sample[[2]]),
    group = p.sample[[3]]
  )

  if (ncol(p.sample) < 3L) {
    stop("`p.sample` must have at least 3 columns: entry, event time, group.")
  }

  if (length(p) != length(q)) {
    if (length(p) == 1L) {
      p <- rep(p, length(q))
    } else if (length(q) == 1L) {
      q <- rep(q, length(p))
    } else {
      stop("`p` and `q` must have the same length, or one of them must have length 1.")
    }
  }

  groups <- unique(as.character(p.sample[, 3]))
  K <- length(groups)

  if (K < 2L) {
    stop("At least two groups are required.")
  }

  largest.holes <- rep(NA_real_, K)

  for (i in seq_len(K)) {
    samp_i <- p.sample[p.sample[, 3] == groups[i], , drop = FALSE]
    est_i <- lt.estimator(samp_i, conf.int = FALSE, holes = holes)

    if (!is.null(est_i$hole)) {
      largest.holes[i] <- est_i$hole
    }
  }

  if (holes == "conditional" && any(!is.na(largest.holes))) {
    t0 <- max(largest.holes, na.rm = TRUE)
    p.sample <- p.sample[p.sample[, 2] > t0, , drop = FALSE]
  }

  times <- sort(unique(p.sample[, 2]))
  D <- length(times)

  if (D == 0L) {
    stop("No observed event times remain after applying the selected hole handling.")
  }

  d <- n <- matrix(0, nrow = D, ncol = K)
  colnames(d) <- colnames(n) <- groups

  eps <- (abs(times) + 1) * .Machine$double.eps

  for (i in seq_len(K)) {
    samp <- p.sample[p.sample[, 3] == groups[i], , drop = FALSE]

    if (nrow(samp) == 0L) next

    idx <- match(samp[, 2], times)
    idx <- idx[!is.na(idx)]
    if (length(idx) > 0L) {
      d[, i] <- tabulate(idx, nbins = D)
    }

    enter <- findInterval(times, sort(samp[, 1]), rightmost.closed = TRUE)
    exit  <- findInterval(times - eps, sort(samp[, 2]), rightmost.closed = TRUE)
    n[, i] <- enter - exit
  }

  d.pooled <- rowSums(d)
  n.pooled <- rowSums(n)

  keep <- n.pooled > 0 & d.pooled > 0
  d <- d[keep, , drop = FALSE]
  n <- n[keep, , drop = FALSE]
  d.pooled <- d.pooled[keep]
  n.pooled <- n.pooled[keep]
  D <- length(d.pooled)

  if (D == 0L) {
    stop("No valid event times remain for the test.")
  }

  est.pooled <- lt.estimator(p.sample, conf.int = FALSE, holes = holes)
  phi <- c(1, est.pooled$estimation[1:(D - 1)])

  corr <- numeric(D)
  ok <- n.pooled > 1
  corr[ok] <- (n.pooled[ok] - d.pooled[ok]) / (n.pooled[ok] - 1)

  ss <- length(p)
  output <- matrix(NA_real_, nrow = ss, ncol = 2)
  colnames(output) <- c("statistic", "p.value")
  rownames(output) <- paste0("p=", p, ", q=", q)

  Z <- numeric(K - 1L)
  varcov <- matrix(0, nrow = K - 1L, ncol = K - 1L)

  aa <- d.pooled / n.pooled
  bb <- d.pooled / (n.pooled^2)

  for (l in seq_len(ss)) {
    Z[] <- 0
    varcov[,] <- 0

    weights <- phi^p[l] * (1 - phi)^q[l]

    for (i in seq_len(K - 1L)) {
      Z[i] <- sum(weights * (d[, i] - aa * n[, i]))

      varcov[i, i] <- sum(weights^2 * n[, i] * aa * (1 - n[, i] / n.pooled) * corr)

      if (i < (K - 1L)) {
        for (j in (i + 1L):(K - 1L)) {
          varcov[i, j] <- -sum(weights^2 * bb * n[, i] * n[, j] * corr)
          varcov[j, i] <- varcov[i, j]
        }
      }
    }

    if (qr(varcov)$rank < (K - 1L)) {
      output[l, ] <- c(NA_real_, NA_real_)
    } else {
      stat <- drop(t(Z) %*% solve(varcov, Z))
      output[l, 1] <- stat
      output[l, 2] <- 1 - pchisq(stat, df = K - 1L)
    }
  }

  if (plot.curves) {
    xmin <- min(p.sample[, 2])
    xmax <- max(p.sample[, 2])

    plot(lt.estimator(p.sample[p.sample[, 3] == groups[1], , drop = FALSE],
                      conf.int = FALSE, holes = holes),
         col = 2, xlim = c(xmin, xmax), xlab = "t",
         ylab = "Survival function estimators")

    if (K >= 2L) {
      for (i in 2:K) {
        lines(lt.estimator(p.sample[p.sample[, 3] == groups[i], , drop = FALSE],
                           conf.int = FALSE, holes = holes),
              col = i + 1)
      }
    }
  }

  return(output)
}

