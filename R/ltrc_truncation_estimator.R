#' Estimator of the truncation distribution function for data subject to
#' left truncation and right censoring
#'
#' It computes the estimator of the truncation distribution function for data subject to
#' left truncation and right censoring
#'
#' @param o.sample Matrix or data frame with columns (U, Y, Delta).
#' @param holes One of 'add', 'conditional' or 'holes' (default). How to deal with holes
#' in the estimation of the target variable survival function. See
#' [ltrc.estimator()] for details.
#'
#' @return An object with the observed truncation times, the estimator at those
#' times, the jumps of the estimator and an estimator of the nontruncation probability
#'
#' @details
#' The estimator of the distribution function of the truncation variable is given by
#' \deqn{ G_n(t) = \sum_i \frac{S_n(u_i)^{-1}}{\sum_j S_n(u_j)^{-1} } \mathbb{I}[u_i \leq t],}
#' where \eqn{u_1, \ldots, u_m} are the distinct truncation times and \eqn{S_n} is the
#' survival function estimator of the target variable.
#'
#' The nontruncation probability is defined as
#' \deqn{ \gamma = \mathbb{P}(U \leq X) = \int G(z)dF(z),}where \eqn{F} is the distribution function of the target
#' variable. Then, it can be estimated as
#' \deqn{ \gamma_n = \sum_{i=1}^m \varphi_i G_n(y_i),}where \eqn{y_1, \ldots, y_m} are
#' the distinct observed event times and \eqn{\varphi_1, \ldots, \varphi_m} are
#' the weights that the estimator \eqn{S_n} assigns to \eqn{y_1, \ldots, y_m}.
#'
#' In practice, it is possible to observe truncation times larger than the largest
#' event time. This would break the formula of the estimator \eqn{G_{n}} defined
#' above. In such a case, following the recommendations in Wang (1991), we reduce
#' the distribution function of the truncation variable to a conditional distribution
#' on the identifiable region, by trimming a small amount of data. This is analogous to
#' conditioning to times largest than the largest hole when estimating the distribution
#' function of the target variable.
#'
#' @references
#' Wang, M.C. (1991). Nonparametric Estimation from Cross-Sectional Survival Data. Journal of
#' the American Statistical Association, 86:130-143.
#'
#' @examples
#' aa <- ltrc.sim(n=5000,rtrunc=rnorm,trunc.args = list(4,1),rtarget=rnorm,
#' target.args=list(4,1),rcens=rexp,cens.args=list(1))
#' head(aa)
#' plot(ltrc.truncation.estimator(aa))
#' curve(pnorm(x,4,1),add=TRUE,col=2)
#'
#' if (requireNamespace("KMsurv", quietly = TRUE)) {
#' data("channing", package = "KMsurv")
#'
#' o.sample <- channing[channing$gender==1,c(3,4,2)] # subgroup of men
#' plot(ltrc.estimator(o.sample)) # hole in the second event time
#' plot(ltrc.truncation.estimator(o.sample,holes='add'))
#' lines(ltrc.truncation.estimator(o.sample,holes='conditional'),col=2)
#' }
#'
#' @export

ltrc.truncation.estimator <- function(o.sample, holes = c("holes", "add", "conditional")) {

  holes <- match.arg(holes)

  if (!is.matrix(o.sample) && !is.data.frame(o.sample)) {
    stop("'o.sample' must be a matrix or data.frame.")
  }
  o.sample <- as.matrix(o.sample)

  if (ncol(o.sample) < 3L) {
    stop("'o.sample' must have at least 3 columns: (U, Y, delta).")
  }

  est <- ltrc.estimator(o.sample, holes = holes)

  if (length(est$fail.time) == 0L) {
    lista <- list(
      trunc.time = numeric(0),
      estimation = numeric(0),
      prob = numeric(0),
      gamma.hat = NA_real_,
      tau = NA_real_
    )
    class(lista) <- "t.estimators"
    return(lista)
  }

  u.ord <- sort(unique(o.sample[, 1]))
  tau <- max(est$fail.time)

  keep_id <- u.ord <= tau
  u.ord.valid <- u.ord[keep_id]

  if (length(u.ord.valid) == 0L) {
    lista <- list(
      trunc.time = numeric(0),
      estimation = numeric(0),
      prob = numeric(0),
      gamma.hat = NA_real_,
      tau = tau
    )
    class(lista) <- "t.estimators"
    return(lista)
  }

  eps <- (abs(u.ord.valid) + 1) * .Machine$double.eps
  ind <- findInterval(u.ord.valid - eps, est$fail.time, rightmost.closed = TRUE)

  s.left <- c(1, est$estimation)[ind + 1L]

  keep_pos <- s.left > 0

  if (!any(keep_pos)) {
    lista <- list(
      trunc.time = numeric(0),
      estimation = numeric(0),
      prob = numeric(0),
      gamma.hat = NA_real_,
      tau = tau
    )
    class(lista) <- "t.estimators"
    return(lista)
  }

  u.ord.valid <- u.ord.valid[keep_pos]
  s.left <- s.left[keep_pos]

  omega.pre <- 1 / s.left
  prob <- omega.pre / sum(omega.pre)
  estimation <- cumsum(prob)
  gamma.hat <- 1 / mean(omega.pre)


  lista <- list(
    trunc.time = u.ord.valid,
    estimation = estimation,
    prob = prob,
    gamma.hat = gamma.hat,
    tau = tau
  )
  class(lista) <- "t.estimators"
  lista
}






