#' Tests for the \eqn{k}-sample problem under right censoring
#'
#' Computes \eqn{k}-sample tests adapted to right
#' censoring and the corresponding \eqn{p}-value.
#'
#' @param p.sample Matrix or data frame with columns (Y, \eqn{\Delta}, group).
#' @param weights Vector of weights of the estimator under the null hypothesis.
#' If no vector is provided, weights proportional to each sample size are assigned automatically
#' by the function. Not used when `tests="logrank"`.
#' @param B Number of bootstrap replications to approximate the \eqn{p}-values of the Kolmogorov--Smirnov
#' and/or Cramér--von Mises tests.
#' @param tests Tests to perform. One or more of Kolmogorov--Smirnov ('ks'), Cramér--von Mises ('cvm'),
#' or weighted log-rank-type tests ('logrank').
#' @param plot.curves Logical, FALSE by default. If set to TRUE, the \eqn{k} survival functions
#' to be compared are depicted.
#' @param p,q Weights of the log-rank-type tests. The choice \eqn{p=0=q} yields the classical
#' log-rank test.
#' @param seed Seed for the bootstrap resamples.
#' @param keep.boot Logical. If `TRUE`, bootstrap replicates are included in the
#'   returned object as `stat.boot`. Default is `FALSE`.
#'
#' @return A list with the following components:
#' * `statistic`: numerical value of the Kolmogorov--Smirnov and/or Cramér--von Mises statistics.
#' * `p.value`: approximated \eqn{p}-value of the tests.
#' * `logrank`: a matrix with weights `p` and `q`, test statistics, and \eqn{p}-values.
#' * `B`: number of bootstrap replications.
#' * `tests`: tests performed, as specified in the input.
#' * `curves`: distribution function estimators of the \eqn{k} samples to be compared.
#' * `groups`: group indicators.
#' * `stat.boot`: evaluation of the tests in the bootstrap resamples. Only when `keep.boot = TRUE`.
#'
#' The `summary` method returns a data frame with one row per test
#' and the following columns:
#'   \describe{
#'     \item{Test}{Name of the test. For the weighted log-rank-type tests, the weights
#'     `p` and `q` are detailed.}
#'     \item{Statistic}{Observed value of the test statistic.}
#'     \item{p.value}{Corresponding \eqn{p}-value.}
#'  }
#'
#'
#' @details
#' This function collects the three tests implemented in this package: Kolmogorov--Smirnov,
#' Cramér--von Mises and weighted log-rank-type tests, (Klein and Moeschberger, 2003)
#' adapted to right censoring. The \eqn{p}-values of the first two are
#' approximated by the obvious bootstrap (Efron, 1981), and the \eqn{p}-value of the
#' log-rank tests is computed by the \eqn{\chi^2_{k-1}} asymptotic distribution
#' of the test statistic, where \eqn{k} is the number of groups.
#'
#' @references
#' Klein, J.P., and Moeschberger, M. L. (2003). Survival Analysis: Techniques for Censored
#' and Truncated Data. Springer, New York.
#'
#' Efron, B. (1981). Censored data and the bootstrap. Journal of the American
#' Statistical Association, 76: 312-319.
#'
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#' library(survival)
#'
#' data <- data.frame(y=stanford2$time,delta=stanford2$status,
#' group=as.numeric(stanford2$age<=35)+as.numeric(stanford2$age<=47))
#' aa <- rc.test(p.sample=data,weights=table(data$group)/nrow(data),
#' plot.curves=TRUE,tests=c('ks','cvm','logrank'),
#' B = 200, p=c(0,1,2), q=0 )
#'
#' plot(aa)
#' summary(aa)
#' }
#'
#' @export

rc.test <- function(p.sample, weights  = NULL,
                    tests = c("ks", "cvm", "logrank"),
                    B = 500,
                    plot.curves = FALSE,
                    p = 0, q = 0,
                    seed = NULL,
                    keep.boot = FALSE) {

  tests <- match.arg(tests, c("ks", "cvm", "logrank"), several.ok = TRUE)

  if (!is.null(seed)) {
    set.seed(seed)
  }


  if (!is.matrix(p.sample) && !is.data.frame(p.sample)) {
    stop("'p.sample' must be a matrix or data.frame.")
  }
  p.sample <- data.frame(
    Y = as.numeric(p.sample[, 1]),
    delta = as.numeric(p.sample[, 2]),
    group = p.sample[, 3]
  )

  if (ncol(p.sample) < 3L) {
    stop("p.sample must have at least 3 columns: (Y, delta, group).")
  }

  if (any(!is.finite(p.sample[, 1]))) {
    stop("Observed times 'Y' must be finite.")
  }

  if (any(p.sample[, 1] < 0)) {
    stop("Observed times 'Y' must be non-negative.")
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

  obs_stats <- list()

  if ("ks" %in% tests) {
    obs_stats$ks <- rc.ks(p.sample, weights = weights,
                          plot.curves = plot.curves)$statistic
    if (plot.curves){
      legend <- numeric(K)
      for(i in seq_len(K)){legend[i] <- paste('Group',ind.group[i])}
      legend('topright',col=1:K,lwd=2,legend=legend)
    }
  }

  if ("cvm" %in% tests) {
    obs_stats$cvm <- rc.cvm(p.sample, weights = weights,
                            plot.curves = FALSE)$statistic
  }

  logrank_results <- NULL

  if ("logrank" %in% tests) {

    lr_out <- rc.logrank(p.sample, p = p, q = q,
                          plot.curves = FALSE)

    logrank_results <- data.frame(
      statistic = lr_out[, 1],
      p.value   = lr_out[, 2]
    )

    rownames(logrank_results) <- paste0("p=", p, ", q=", q)
  }


  F0 <- vector("list", K)
  C0 <- vector("list", K)
  NN <- numeric(K)

  for (j in seq_len(K)) {

    sj <- p.sample[p.sample[, 3] == ind.group[j], 1:3, drop = FALSE]
    NN[j] <- nrow(sj)

    F0[[j]] <- rc.estimator(sj[, 1:2, drop = FALSE], conf.int = FALSE)

    C0[[j]] <- rc.estimator(
      cbind(sj[, 1], 1 - sj[, 2]),
      conf.int = FALSE
    )

    if (length(F0[[j]]$fail.time) == 0L) {
      stop(paste("Group", j, "has no observed failure times."))
    }

    if (length(C0[[j]]$fail.time) == 0L) {
      C0[[j]]$fail.time <- Inf
      C0[[j]]$prob <- 1
      warning(paste("Group", j, "has no censoring times."))
    } else if (length(C0[[j]]$fail.time) == 1L){
      C0[[j]]$fail.time <- rep( C0[[j]]$fail.time , 2)
      C0[[j]]$prob <- rep( 1/2, 2)
    }
  }

  boot_stats <- lapply(obs_stats, function(x) numeric(B))

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
            size = length(idx),
            replace = TRUE,
            prob = F0[[k]]$prob
          )
        }
      }

      # Group-specific censoring bootstrap
      ctimes.boot <- sample(
        C0[[j]]$fail.time,
        size = n.j,
        replace = TRUE,
        prob = C0[[j]]$prob
      )

      y.boot <- pmin(evtimes.boot, ctimes.boot)
      delta.boot <- as.numeric(evtimes.boot <= ctimes.boot)

      resample[[j]] <- data.frame(
        y = y.boot,
        delta = delta.boot,
        group = rep(ind.group[j], n.j)
      )
    }

    resample.df <- do.call(rbind, resample)
    row.names(resample.df) <- NULL

    if ("ks" %in% tests) {
      boot_stats$ks[b] <- rc.ks(resample.df, weights = weights,
                                plot.curves = FALSE)$statistic
    }

    if ("cvm" %in% tests) {
      boot_stats$cvm[b] <- rc.cvm(resample.df, weights = weights,
                                  plot.curves = FALSE)$statistic
    }
  }

  pvalues <- mapply(function(obs, boot) {
    (1 + sum(boot >= obs)) / (B + 1)
  }, obs_stats, boot_stats)

  out <- list(
    statistic = obs_stats,
    p.value = as.list(pvalues),
    logrank = logrank_results,
    B = B,
    tests = tests,
    curves = F0,
    groups = ind.group
  )

  if (isTRUE(keep.boot)) {
    out$stat.boot <- boot_stats
  }

  class(out) <- "ksamples_test"

  out
}


