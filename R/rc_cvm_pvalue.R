#' Cramér--von Mises-type test for the \eqn{k}-sample problem under right censoring
#'
#' Computes the Kolmogorov--Smirnov test adapted to left truncation and right
#' censoring and the corresponding \eqn{p}-value.
#'
#' @param p.sample Matrix or data frame with columns (U, Y, \eqn{\Delta}, group).
#' @param weights Vector of weights of the estimator under the null hypothesis.
#' If no vector is provided, weights proportional to each sample size are assigned automatically
#' by the function. Not used when `tests="logrank"`.
#' @param B Number of bootstrap replications to approximate the \eqn{p}-value.
#' @param keep.boot Logical, default to FALSE. If TRUE, the function returns the evaluation
#' of the Cramér--von Mises statistic on the bootstrap resamples.
#'
#' @return A list with the following:
#' * `statistic`: Cramér--von Mises statistic evaluated on the the samples.
#' * `pvalue`: \eqn{p}-value of the test approximated via bootstrap.
#' * `curves`: A list with the estimators of the \eqn{k} distribution functions to
#' be compared.
#' * `groups`: Group indicator.
#' * `stat.boot`: The Cramér--von Mises statistic evaluated from the bootstrap resamples.
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
#' The \eqn{p}-value is approximated by means of the obvious bootstrap accounting for
#' left truncation and right censoring, studied in Efron (1981).
#'
#' @references
#' Efron, B. (1981). Censored data and the bootstrap. Journal of the American
#' Statistical Association, 76: 312-319.
#'
#' Kiefer, J. (1959). K-sample analogues of the Kolmogorov--Smirnov and Cramér--von Mises
#' tests. The Annals of Mathematical Statistics, 30:420-447.
#'
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#' library(survival)
#' data <- data.frame(y=stanford2$time,
#' delta=stanford2$status,
#' group=as.numeric(stanford2$age<=35)+as.numeric(stanford2$age<=47))
#' aa <- rc.pv.cvm(data,weights=table(data$group)/nrow(data), B = 200)
#' plot(aa)
#' aa$pvalue
#' }
#'
#' @export

rc.pv.cvm <- function(p.sample, weights = NULL, B = 500, keep.boot = FALSE) {

  if (!is.matrix(p.sample) && !is.data.frame(p.sample)) {
    stop("'p.sample' must be a matrix or data.frame.")
  }

  p.sample <- data.frame(
    Y = as.numeric(p.sample[[1]]),
    delta = as.numeric(p.sample[[2]]),
    group = p.sample[[3]]
  )

  if (ncol(p.sample) < 3L) {
    stop("p.sample must have at least 3 columns: (Y, delta, group).")
  }

  if (any(p.sample[, 2] != 0 & p.sample[, 2] != 1)) {
    stop("delta must be 0/1.")
  }

  ind.group <- sort(unique(p.sample[, 3]))
  K <- length(ind.group)

  if (is.null(weights)) {
    N <- tabulate(factor(p.sample[,3]))
    weights <- N / sum(N)
  } else { if (!is.numeric(weights)) {
    stop("'weights' must be numeric.")
  }
  if (length(weights) != K) {
    stop("'weights' must have one entry per group.")
  }
  if (any(!is.finite(weights)) || any(weights <= 0)) {
    stop("'weights' must be finite and strictly positive.")
  }
  if (abs(sum(weights) - 1) > sqrt(.Machine$double.eps)) {
    stop("'weights' must sum to 1.")
  }
  }

  F0 <- vector("list", K)
  QQ <- vector("list", K)
  NN <- numeric(K)

  for (j in seq_len(K)) {
    sj <- p.sample[p.sample[, 3] == ind.group[j], , drop = FALSE]
    NN[j] <- nrow(sj)

    F0[[j]] <- rc.estimator(sj, conf.int = FALSE)
    QQ[[j]] <- rc.estimator(cbind(sj[, 1], 1 - sj[, 2]), conf.int = FALSE)

    if (length(F0[[j]]$fail.time) == 0L) {
      stop(paste("Group", j, "has no observed failure times."))
    }
    if (length(QQ[[j]]$fail.time) == 0L) {
      QQ[[j]]$fail.time <- Inf
      QQ[[j]]$prob <- 1
      warning(paste("Group", j, "has no censoring times: classical k-sample problem."))
    } else if(length(QQ[[j]]$fail.time) == 1L){
      QQ[[j]]$fail.time <- rep(QQ[[j]]$fail.time , 2 )
      QQ[[j]]$prob <- rep(1/2,2)
    }
  }

  stat.boot <- numeric(B)

  for (b in seq_len(B)) {
    resample <- vector("list", K)

    for (j in seq_len(K)) {
      n.j <- NN[j]

      ind.evtimes <- sample(
        seq_len(K),
        size = n.j,
        replace = TRUE,
        prob = weights
      )

      evtimes.boot <- numeric(n.j)
      for (k in seq_len(K)) {
        idx <- which(ind.evtimes == k)
        if (length(idx) > 0L) {
          evtimes.boot[idx] <- sample(
            F0[[k]]$fail.time,
            length(idx),
            replace = TRUE,
            prob = F0[[k]]$prob
          )
        }
      }

      cens.boot <- sample(
        QQ[[j]]$fail.time,
        size = n.j,
        replace = TRUE,
        prob = QQ[[j]]$prob
      )

      resample[[j]] <- data.frame(
        y = pmin(evtimes.boot, cens.boot),
        delta = as.numeric(evtimes.boot <= cens.boot),
        group = rep(ind.group[j], n.j)
      )
    }

    resample.df <- do.call(rbind, resample)
    row.names(resample.df) <- NULL

    stat.boot[b] <- rc.cvm(resample.df, weights = weights)$statistic
  }

  stat.0 <- rc.cvm(p.sample, weights = weights)$statistic
  pvalue <- (1+sum(stat.boot >= stat.0))/(B+1)

  out <- list(
    statistic = stat.0,
    pvalue = pvalue,
    curves = F0,
    groups = ind.group
  )

  if (isTRUE(keep.boot)) {
    out$stat.boot <- stat.boot
  }

  class(out) <- "ksamples_test"

  out
}
