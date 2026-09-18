#' Automatic test selector for survival data
#'
#' Performs tests for data subject to left truncation, right censoring or both.
#'
#' @param time Numeric vector of event or follow-up times.
#' @param status Event indicator (1 = event, 0 = censored).
#' @param entry Optional entry time (for left truncation).
#' @param group Group indicator (factor or vector).
#' @param weights Vector of weights of the estimator under the null hypothesis.
#' If no vector is provided, weights proportional to each sample size are assigned automatically
#' by the function. Not used when `tests="logrank"`.
#' @param tests Character: "ks", "cvm", "logrank".
#' @param B Number of bootstrap replications for approximating the \eqn{p}-value.
#' @param plot.curves Logical, default to FALSE. If TRUE, the function depicts the functions
#' to be compared.
#' @param holes One of 'holes' (default), 'add', or 'conditional'. See [ltrc.estimator()] for details.
#' @param p,q Weights for the log-rank-type tests to be performed. If NULL, the classical log-rank.
#' test is performed.
#' @param seed Seed for the bootstrap statistics.
#' @param keep.boot Logical, default to FALSE. If TRUE, the statistics computed from the bootstrap resamples
#' are stored.
#'
#' @return A list with the following components:
#' * `statistic`: numerical value of the Kolmogorov--Smirnov and/or Cramér--von Mises statistics.
#' * `p.value`: approximated \eqn{p}-value of the tests.
#' * `stat.boot`: bootstrap statistics for the \eqn{p}-value approximation.
#' * `logrank`: a matrix with weights p and q, test statistics, and \eqn{p}-values.
#' * `B`: number of bootstrap replications.
#' * `tests`: tests performed, as specified in the input.
#' * `curves`: distribution function estimators of the \eqn{k} samples to be compared.
#' * `groups`: group indicators.
#' * `holes`: method for handling holes (if appropriate).
#' * `stat.boot`: bootstrap statistics (if `keep.boot = TRUE`).
#'
#' @examples
#' if (requireNamespace("KMsurv", quietly = TRUE)) {
#' data("channing", package = "KMsurv")
#'
#' aa <- ksample.test(time = channing$age,entry = channing$ageentry,
#' status = channing$death, group = channing$gender,
#' weights=rep(1/2,2), B = 200, holes='conditional',
#' tests=c('ks','cvm','logrank'),plot.curves=TRUE,
#' p=c(0,1,0.5,0),q=c(0,0,0.5,1),keep.boot = FALSE)
#'
#' summary(aa)
#' }
#'
#' @export

ksample.test <- function(
    time,
    entry = NULL,
    status = NULL,
    group,
    weights = NULL,
    tests = c("ks", "cvm", "logrank"),
    B = 500,
    plot.curves = FALSE,
    holes = c("holes", "add", "conditional"),
    p = 0,
    q = 0,
    seed = NULL,
    keep.boot = FALSE
) {

  tests <- match.arg(
    tests,
    choices = c("ks", "cvm", "logrank"),
    several.ok = TRUE
  )

  holes <- match.arg(
    holes,
    choices = c("holes", "add", "conditional")
  )

  if (!is.numeric(time))
    stop("'time' must be numeric")

  n <- length(time)

  if (n == 0L)
    stop("'time' cannot be empty")

  if (anyNA(time) || any(!is.finite(time)))
    stop("'time' must contain only finite, nonmissing values")

  if (length(group) != n)
    stop("'time' and 'group' must have the same length")

  if (anyNA(group))
    stop("'group' cannot contain missing values")

  groups <- sort(unique(group))
  K <- length(groups)

  if (K < 2L)
    stop("At least two groups are required")

  if (!is.null(status)) {
    if (length(status) != n)
      stop("'status' and 'time' must have the same length")

    if (anyNA(status))
      stop("'status' cannot contain missing values")

    if (is.factor(status))
      status <- as.character(status)

    suppressWarnings(status.numeric <- as.numeric(status))

    if (anyNA(status.numeric) ||
        any(!status.numeric %in% c(0, 1))) {
      stop("'status' must contain only 0 and 1")
    }

    status <- status.numeric
  }

  if (!is.null(entry)) {
    if (!is.numeric(entry))
      stop("'entry' must be numeric")

    if (length(entry) != n)
      stop("'entry' and 'time' must have the same length")

    if (anyNA(entry) || any(!is.finite(entry)))
      stop("'entry' must contain only finite, nonmissing values")

    if (any(entry > time))
      stop("All observations must satisfy 'entry <= time'")
  }

  B <- as.integer(B)

  if (length(plot.curves) != 1L ||
      is.na(plot.curves) ||
      !is.logical(plot.curves)) {
    stop("'plot.curves' must be TRUE or FALSE")
  }

  if (is.null(weights)) {
    group.index <- match(group, groups)
    g <- tabulate(group.index, nbins = K)
    weights <- g / n
  } else {
    if (!is.numeric(weights) ||
        length(weights) != K ||
        anyNA(weights) ||
        any(!is.finite(weights)) ||
        any(weights <= 0)) {
      stop(
        "'weights' must be finite, positive, and have one value per group"
      )
    }

    if (abs(sum(weights) - 1) >
        sqrt(.Machine$double.eps)) {
      stop("'weights' must sum to 1")
    }
  }

  if (!is.null(entry) && !is.null(status)) {

    p.sample <- data.frame(
      U = entry,
      Y = time,
      delta = status,
      group = group
    )

    return(
      ltrc.test(
        p.sample = p.sample,
        weights = weights,
        tests = tests,
        B = B,
        plot.curves = plot.curves,
        holes = holes,
        p = p,
        q = q,
        seed = seed,
        keep.boot = keep.boot
      )
    )
  }

  if (!is.null(entry) && is.null(status)) {

    p.sample <- data.frame(
      U = entry,
      Y = time,
      group = group
    )

    return(
      lt.test(
        p.sample = p.sample,
        weights = weights,
        tests = tests,
        B = B,
        plot.curves = plot.curves,
        holes = holes,
        p = p,
        q = q,
        seed = seed,
        keep.boot = keep.boot
      )
    )
  }

  if (is.null(entry) && !is.null(status)) {

    p.sample <- data.frame(
      Y = time,
      delta = status,
      group = group
    )

    return(
      rc.test(
        p.sample = p.sample,
        weights = weights,
        tests = tests,
        B = B,
        plot.curves = plot.curves,
        p = p,
        q = q,
        seed = seed,
        keep.boot = keep.boot
      )
    )
  }

  stop(
    paste0(
      "Complete uncensored and nontruncated data are not currently ",
      "supported: provide 'entry', 'status', or both"
    )
  )
}




