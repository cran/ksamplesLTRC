#' Cramér--von Mises-type test for the \eqn{k}-sample problem under left truncation
#'
#' Computes the Cramér--von Mises statistic adapted to left truncation,
#' and the corresponding \eqn{p}-value.
#'
#' @param p.sample Matrix or data frame with columns (U, Y, \eqn{\Delta}, group).
#' @param weights Vector of weights of the estimator under the null hypothesis.
#' If no vector is provided, weights proportional to each sample size are assigned automatically
#' by the function. Not used when `tests="logrank"`.
#' @param B Number of bootstrap replications to approximate the \eqn{p}-value.
#' @param plot.curves Logical, default to FALSE. If TRUE, the \eqn{k} curves to be
#' compared are depicted.
#' @param holes How to handle holes if present in a sample. One of `holes` (default),
#' `conditional` or `add`. See details in [ltrc.estimator()].
#' @param keep.boot Logical, default to FALSE. If TRUE, the bootstrap statistics are
#' returned in the output of the function.
#' @param seed Seed for the bootstrap resamples.
#'
#' @return A list with the following elements:
#' * `statistic`: numerical value of the Cramér--von Mises test statistic
#' * `pvalue`: \eqn{p}-value of the test
#' * `curves`: A list with the estimators of the \eqn{k} distribution functions to be
#' compared.
#' * `groups`: Group indicator
#' * `stat.boot`: If `keep.boot = TRUE`, the function returns the statistic evaluated
#' in the bootstrap resamples.
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
#' The \eqn{p}-value is approximated by means of the obvious bootstrap accounting for
#' left truncation, as a particular case of the obvious bootstrap proposed for left-truncated
#' and right-censored data in Bilker and Wang (1997). The exact
#' algorithm for the two-sample problem can be consulted in Lago et al. (2025).
#'
#' @references
#' Bilker, W.B., and Wang, M.C. (1997). Bootstrapping left truncated
#' and right censored data. Communications in Statistics - Simulation and Computation, 26:141-171.
#'
#' Kiefer, J. (1959). K-sample analogues of the Kolmogorov--Smirnov and Cramér--von Mises
#' tests. The Annals of Mathematical Statistics, 30:420-447.
#'
#' Lago, A., de Uña-Álvarez, J., & Pardo-Fernández, J. C. (2025). A Kolmogorov–Smirnov-type test
#' for the two-sample problem with left-truncated data. TEST, 34, 69-90.
#'
#' @examples
#' if (requireNamespace("mvna", quietly = TRUE)) {
#' data("abortion", package = "mvna")
#' pooled.sample <- abortion[abortion$cause==3, 2:4]
#' head(pooled.sample)
#' aa <- lt.pv.cvm(p.sample=pooled.sample,weights=rep(1/2,2),plot.curves=TRUE, B = 200)
#' aa$pvalue
#' plot(aa)
#' }
#'
#'
#' @export

lt.pv.cvm <- function(
    p.sample,
    weights = NULL,
    B = 500,
    plot.curves = FALSE,
    holes = c("holes", "add", "conditional"),
    seed = NULL,
    keep.boot = FALSE
) {

  holes <- match.arg(holes)

  if (!is.null(seed)) {
    set.seed(seed)
  }

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
    group = p.sample[[3]]
  )

  if (ncol(p.sample) < 3L) {
    stop("'p.sample' must have at least 3 columns: (U, X, group).")
  }

  if (any(p.sample[, 1] > p.sample[, 2])) {
    stop("All observations must satisfy U <= X.")
  }

  ind.group <- sort(unique(p.sample[, 3]))
  K <- length(ind.group)

  if (is.null(weights)) {
    N <- tabulate(factor(p.sample[,3]))
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
      sj <- x[
        x[, 3] == ind.group[j],
        1:2,
        drop = FALSE
      ]

      est <- catch_hole_warning(
        lt.estimator(
          sj,
          conf.int = FALSE,
          holes = "holes"
        )
      )

      out[j] <- if (is.null(est$hole)) {
        NA_real_
      } else {
        est$hole
      }
    }

    out
  }

  hole.ind <- detect_hole(p.sample)

  if (holes == "conditional" && any(!is.na(hole.ind))) {
    largest.hole <- max(hole.ind, na.rm = TRUE)

    p.sample <- p.sample[
      p.sample[, 2] > largest.hole,
      ,
      drop = FALSE
    ]
  }

  holes.internal <- if (holes == "add") {
    "add"
  } else {
    "holes"
  }

  stat.0 <- catch_hole_warning(
    lt.cvm(
      p.sample = p.sample,
      weights = weights,
      plot.curves = plot.curves,
      holes = holes.internal
    )$statistic
  )

  if (plot.curves) {
    legend(
      "topright",
      col = seq_len(K),
      lwd = 2,
      legend = paste("Group", ind.group)
    )
  }

  F0 <- vector("list", K)
  GG <- vector("list", K)
  NN <- numeric(K)

  for (j in seq_len(K)) {
    sj <- p.sample[
      p.sample[, 3] == ind.group[j],
      1:2,
      drop = FALSE
    ]

    NN[j] <- nrow(sj)

    F0[[j]] <- catch_hole_warning(
      lt.estimator(
        sj,
        conf.int = FALSE,
        holes = holes.internal
      )
    )

    GG[[j]] <- catch_hole_warning(
      lt.truncation.estimator(
        sj,
        holes = holes
      )
    )

    if (length(F0[[j]]$fail.time) == 0L) {
      stop(paste("Group", j, "has no observed target times."))
    }

    if (length(GG[[j]]$trunc.time) == 0L) {
      stop(paste("Group", j, "has no truncation times."))
    }
  }

  stat.boot <- numeric(B)

  for (b in seq_len(B)) {
    resample <- vector("list", K)

    for (j in seq_len(K)) {
      n.j <- NN[j]

      trunc.boot <- sample(
        GG[[j]]$trunc.time,
        size = n.j,
        replace = TRUE,
        prob = GG[[j]]$prob
      )

      ind.xtimes <- sample(
        seq_len(K),
        size = n.j,
        replace = TRUE,
        prob = weights
      )

      xtimes.boot <- numeric(n.j)

      for (k in seq_len(K)) {
        idx <- ind.xtimes == k

        if (any(idx)) {
          xtimes.boot[idx] <- sample(
            F0[[k]]$fail.time,
            size = sum(idx),
            replace = TRUE,
            prob = F0[[k]]$prob
          )
        }
      }

      bad <- trunc.boot > xtimes.boot

      while (any(bad)) {
        m <- sum(bad)

        trunc.boot[bad] <- sample(
          GG[[j]]$trunc.time,
          size = m,
          replace = TRUE,
          prob = GG[[j]]$prob
        )

        ind.new <- sample(
          seq_len(K),
          size = m,
          replace = TRUE,
          prob = weights
        )

        xtimes.new <- numeric(m)

        for (k in seq_len(K)) {
          idx <- ind.new == k

          if (any(idx)) {
            xtimes.new[idx] <- sample(
              F0[[k]]$fail.time,
              size = sum(idx),
              replace = TRUE,
              prob = F0[[k]]$prob
            )
          }
        }

        xtimes.boot[bad] <- xtimes.new
        bad <- trunc.boot > xtimes.boot
      }

      resample[[j]] <- data.frame(
        u = trunc.boot,
        x = xtimes.boot,
        group = ind.group[j]
      )
    }

    resample.df <- do.call(rbind, resample)
    row.names(resample.df) <- NULL

    tmp <- try(
      catch_hole_warning(
        lt.cvm(
          p.sample = resample.df,
          weights = weights,
          plot.curves = FALSE,
          holes = holes.internal
        )
      ),
      silent = TRUE
    )

    if (inherits(tmp, "try-error")) {
      stop(
        paste0(
          "Bootstrap replicate ",
          b,
          " failed: ",
          as.character(tmp)
        ),
        call. = FALSE
      )
    }

    stat.boot[b] <- tmp$statistic
  }

  pvalue <- (
    1 + sum(stat.boot >= stat.0)
  ) / (B + 1)

  if (hole.flag && holes == "holes") {
    warning(
      "Holes detected in at least one group; inference performed with holes='holes'."
    )
  }

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
