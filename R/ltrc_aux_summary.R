#' Summary method for k-sample test objects
#'
#' Produces a table with the test name, test statistic and p-value
#' for the procedures computed by [ltrc.test()].
#'
#' @param object An object of class `"ksamples_test"`.
#' @param ... Additional arguments passed to methods.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{Test}{Name of the test (`"KS"`, `"CvM"`, or the corresponding
#'   `(p, q)` values for weighted rank tests).}
#'   \item{Statistic}{Observed test statistic.}
#'   \item{p-value}{Associated p-value.}
#' }
#'
#' @method summary ksamples_test
#' @export
summary.ksamples_test <- function(object, ...) {

  if (!inherits(object, "ksamples_test")) {
    stop("Object is not of class 'ksamples_test'.", call. = FALSE)
  }

  out_list <- list()

  ## KS
  if (!is.null(object$statistic$ks)) {
    out_list[[length(out_list) + 1L]] <- data.frame(
      Test = "KS",
      Statistic = unname(object$statistic$ks),
      p.value = unname(object$p.value$ks),
      stringsAsFactors = FALSE
    )
  }

  ## CvM
  if (!is.null(object$statistic$cvm)) {
    out_list[[length(out_list) + 1L]] <- data.frame(
      Test = "CvM",
      Statistic = unname(object$statistic$cvm),
      p.value = unname(object$p.value$cvm),
      stringsAsFactors = FALSE
    )
  }

  ## Weighted rank tests
  if (!is.null(object$logrank)) {
    lr <- object$logrank

    test_names <- rownames(lr)
    if (is.null(test_names)) {
      test_names <- rep("logrank", nrow(lr))
    }

    out_list[[length(out_list) + 1L]] <- data.frame(
      Test = test_names,
      Statistic = lr[, "statistic"],
      p.value = lr[, "p.value"],
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }

  if (length(out_list) == 0L) {
    return(data.frame(
      Test = character(0),
      Statistic = numeric(0),
      p.value = numeric(0)
    ))
  }

  out <- do.call(rbind, out_list)
  rownames(out) <- NULL
  out
}
