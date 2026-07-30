#' Tests for the \eqn{k}-sample problem under left truncation
#' and right censoring
#'
#' Computes \eqn{k}-sample tests adapted to left truncation and right
#' censoring and the corresponding \eqn{p}-value.
#'
#' @param p.sample Matrix or data frame with columns (U, Y, \eqn{\Delta}, group).
#' @param weights Vector of weights of the estimator under the null hypothesis.
#' If no vector is provided, weights proportional to each sample size are assigned automatically
#' by the function. Not used when `tests="logrank"`.
#' @param B Number of bootstrap replications to approximate the \eqn{p}-values of the Kolmogorov--Smirnov
#' and/or Cramér--von Mises tests.
#' @param tests Tests to perform. One of Kolmogorov--Smirnov (ks), Cramér--von Mises (cvm),
#' or weighted log-rank-type tests (logrank).
#' @param plot.curves Logical, FALSE by default. If set to TRUE, the \eqn{k} survival functions
#' to be compared are depicted.
#' @param holes One of `holes` (default), `conditional` or `add`. If `holes = "add"` and `logrank`
#' is requested, the log-rank test is computed with `holes = "holes"` and a warning is issued. See
#' [ltrc.estimator()] for details.
#' @param p,q Weights of the log-rank-type tests. The choice \eqn{p=0=q} yields the classical
#' log-rank test.
#' @param seed Seed for the bootstrap resamples.
#' @param keep.boot Logical. If `TRUE`, bootstrap replicates are included in the
#'   returned object as `stat.boot`. Default is `FALSE`.
#'
#'
#' @return A list with the following components:
#' * `statistic`: numerical value of the Kolmogorov--Smirnov and/or Cramér--von Mises statistics.
#' * `p.value`: approximated p-value of the tests.
#' * `stat.boot`: bootstrap statistics for the p-value approximation.
#' * `logrank`: a matrix with weights `p` and `q`, test statistics, and \eqn{p}-values.
#' * `B`: number of bootstrap replications.
#' * `tests`: tests performed, as specified in the input.
#' * `curves`: distribution function estimators of the \eqn{k} samples to be compared.
#' * `groups`: group indicators.
#' * `holes`: method for handling holes.
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
#' @details
#' This function collects the three tests implemented in this package: Kolmogorov--Smirnov,
#' Cramér--von Mises (Lago et al., 2026) and weighted log-rank-type tests, (Klein and Moeschberger, 2003)
#' adapted to left truncation and right censoring. The \eqn{p}-values of the first two are
#' approximated by bootstrap, and the \eqn{p}-value of the log-rank tests is computed by
#' the \eqn{\chi^2_{k-1}} asymptotic distribution of the test statistic, where \eqn{k}
#' is the number of groups.
#'
#' @references
#' Klein, J.P., and Moeschberger, M. L. (2003). Survival Analysis: Techniques for Censored
#' and Truncated Data. Springer, New York.
#'
#' Lago, A., Pardo-Fernández, J.C., and de Uña-Álvarez, J. (2026). Kolmogorov--Smirnov
#' and Cramér--von Mises tests for the k-sample problem for left-truncated and
#' right-censored data. Lifetime Data Analysis: 32(2).
#'
#'
#' @examples
#' if (requireNamespace("KMsurv", quietly = TRUE)) {
#' data("channing", package = "KMsurv")
#'
#' data.channing <- data.frame(channing$ageentry, channing$age,
#' channing$death, channing$gender)
#' aa <- ltrc.test(data.channing,weights=rep(1/2,2),
#' holes='conditional',tests=c('ks','cvm','logrank'), B = 200,
#' p=c(0,1,0.5,0),q=c(0,0,0.5,1))
#' plot(aa)
#' summary(aa)
#' }
#'
#' @export

ltrc.test <- function(p.sample, weights  = NULL,
                      tests = c("ks", "cvm", "logrank"),
                      B = 500,
                      plot.curves = FALSE,
                      holes = c("holes","add", "conditional"),
                      p = 0 , q = 0,
                      seed = NULL,
                      keep.boot = FALSE) {

  holes <- match.arg(holes)
  tests <- match.arg(tests, c("ks", "cvm", "logrank"), several.ok = TRUE)

  if (!is.null(seed)) set.seed(seed)

  hole.flag <- FALSE

  catch_hole_warning <- function(expr) {
    withCallingHandlers(
      expr,
      warning = function(w) {
        if (grepl("hole", conditionMessage(w), ignore.case = TRUE)) {
          hole.flag <<- TRUE
          invokeRestart("muffleWarning")
        }
      }
    )
  }

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
      sj <- x[x[,4]==ind.group[j], 1:3, drop=FALSE]
      est <- catch_hole_warning(ltrc.estimator(o.sample = sj, holes = "holes"))
      out[j] <- if (is.null(est$hole)) NA_real_ else est$hole
    }
    out
  }

  hole.ind <- detect_hole(p.sample)

  if (holes == "conditional" && any(!is.na(hole.ind))) {
    largest.hole <- max(hole.ind, na.rm = TRUE)
    p.sample <- p.sample[p.sample[,2] > largest.hole, , drop=FALSE]
  }

  holes.internal <- if (holes == "add") "add" else "holes"

  obs_stats <- list()

  if ("ks" %in% tests) {
    obs_stats$ks <- catch_hole_warning(
      ltrc.ks(p.sample = p.sample,
              weights = weights,
              plot.curves = plot.curves,
              holes = holes.internal)$statistic
    )

    if (plot.curves) {
      legend("topright", col=seq_len(K), lwd=2,
             legend=paste("Group", ind.group))
    }
  }

  if ("cvm" %in% tests) {
    obs_stats$cvm <- catch_hole_warning(
      ltrc.cvm(p.sample = p.sample,
               weights = weights,
               plot.curves = FALSE,
               holes = holes.internal)$statistic
    )
  }

  logrank_results <- NULL

  if ("logrank" %in% tests) {

    holes.logrank <- if (holes == "add") "holes" else holes.internal

    logrank_results <- catch_hole_warning({
      lr <- ltrc.logrank(p.sample = p.sample,
                          p = p,
                          q = q,
                          plot.curves = FALSE,
                          holes = holes.logrank)
      data.frame(statistic = lr[,1], p.value = lr[,2])
    })

    rownames(logrank_results) <- paste0("p=", p, ", q=", q)
  }

  F0 <- vector("list", K)
  QQ <- vector("list", K)
  GG <- vector("list", K)
  NN <- numeric(K)

  for (j in seq_len(K)) {
    sj <- p.sample[p.sample[,4]==ind.group[j], 1:3, drop=FALSE]
    NN[j] <- nrow(sj)

    F0[[j]] <- catch_hole_warning(
      ltrc.estimator(o.sample = sj, holes = holes.internal)
    )

    GG[[j]] <- catch_hole_warning(
      ltrc.truncation.estimator(o.sample = sj, holes = holes)
    )

    QQ[[j]] <- rc.estimator(
      o.sample = cbind(sj[,2]-sj[,1], 1-sj[,3]),
    )

    if (length(QQ[[j]]$fail.time) == 0L) {
      QQ[[j]]$fail.time <- Inf
      QQ[[j]]$prob <- 1
    } else if (length(QQ[[j]]$fail.time) == 1L) {
      QQ[[j]]$fail.time <- rep( QQ[[j]]$fail.time, 2)
      QQ[[j]]$prob <- rep(1/2,2)
    }

  }

  boot_stats <- lapply(obs_stats, function(x) numeric(B))

  for (b in seq_len(B)) {

    resample <- vector("list", K)

    for (j in seq_len(K)) {

      n.j <- NN[j]

      trunc.boot <- sample(GG[[j]]$trunc.time, n.j, TRUE, GG[[j]]$prob)

      ind.evt <- sample(seq_len(K), n.j, TRUE, weights)
      evt.boot <- numeric(n.j)

      for (k in seq_len(K)) {
        idx <- ind.evt == k
        if (any(idx)) {
          evt.boot[idx] <- sample(F0[[k]]$fail.time, sum(idx), TRUE, F0[[k]]$prob)
        }
      }

      bad <- trunc.boot > evt.boot
      while (any(bad)) {
        m <- sum(bad)
        trunc.boot[bad] <- sample(GG[[j]]$trunc.time, m, TRUE, GG[[j]]$prob)

        ind.new <- sample(seq_len(K), m, TRUE, weights)
        evt.new <- numeric(m)

        for (k in seq_len(K)) {
          idx <- ind.new == k
          if (any(idx)) {
            evt.new[idx] <- sample(F0[[k]]$fail.time, sum(idx), TRUE, F0[[k]]$prob)
          }
        }

        evt.boot[bad] <- evt.new
        bad <- trunc.boot > evt.boot
      }

      rc.boot <- sample(QQ[[j]]$fail.time, n.j, TRUE, QQ[[j]]$prob)

      resample[[j]] <- data.frame(
        u = trunc.boot,
        y = pmin(evt.boot, trunc.boot + rc.boot),
        delta = as.numeric(evt.boot <= trunc.boot + rc.boot),
        group = ind.group[j]
      )
    }

    resample.df <- do.call(rbind, resample)

    if ("ks" %in% tests) {
      boot_stats$ks[b] <- catch_hole_warning(
        ltrc.ks(p.sample = resample.df,
                weights = weights,
                plot.curves = FALSE,
                holes = holes.internal)$statistic
      )
    }

    if ("cvm" %in% tests) {
      boot_stats$cvm[b] <- catch_hole_warning(
        ltrc.cvm(p.sample = resample.df,
                 weights = weights,
                 plot.curves = FALSE,
                 holes = holes.internal)$statistic
      )
    }
  }

  pvalues <- mapply(function(obs, boot) mean(c(boot, obs) >= obs),
                    obs_stats, boot_stats)

  if (hole.flag && holes == "holes") {
    warning("Holes detected in at least one group; inference performed with holes='holes'.")
  }

  out <- list(
    statistic = obs_stats,
    p.value = as.list(pvalues),
    logrank = logrank_results,
    B = B,
    tests = tests,
    curves = F0,
    groups = ind.group,
    holes = holes
  )


  if (isTRUE(keep.boot)) {
    out$stat.boot <- boot_stats
  }

  class(out) <- "ksamples_test"
  out
}


