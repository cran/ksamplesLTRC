#' Estimator of the truncation survival function for data subject to
#' left truncation
#'
#' It computes the estimator of the truncation survival function for data subject to
#' left truncation
#'
#' @param o.sample Matrix or data frame with columns (U, X).
#' @param holes Way to handle holes, if present. One of `holes` (default),
#' `add` or `conditional`. See details in [ltrc.estimator()].
#'
#' @return A list with the following:
#' * `trunc.time`: observed truncation times.
#' * `estimation`: estimated survival probability.
#' * `prob`: jumps of the estimator at each observed truncation time.
#' * `gamma`: estimator of the nontruncation probability.
#'
#' @details
#' This estimator is based on the idea that if \eqn{X} is left-truncated by \eqn{U},
#' then \eqn{-U} is left-truncated by \eqn{-X}. In addition, note that
#' \deqn{S_U(u) = 1-S_{_U}(-u),}where \eqn{S_U} is the survival function of \eqn{U}
#' and \eqn{S_{-U}} is the survival function of \eqn{-U}. Thus, since [lt.estimator()]
#' computes the estimator of the survival function,  one can apply the Lynden-Bell
#' estimator to the sample formed by reversing the roles of \eqn{U} and \eqn{X} by
#' merely multiplying by -1.
#'
#' The nontruncation probability is defined as
#' \deqn{\gamma = \mathbb{P} (U \leq X) = \int G(z) dF(z),}
#' where \eqn{U \sim G} and \eqn{X \sim F}. Then, an estimator of \eqn{\gamma}
#' is given by (Woodroofe, 1985)
#' \deqn{ \gamma_n = \sum_{i=1}^m \varphi_i G_n(x_i),}
#' where \eqn{x_1, \ldots, x_m} are the distinct observed event times,
#' \eqn{\varphi_1, \ldots, \varphi_m} are the weights assigned by the Lynden-Bell estimator
#' to the observed event times, and \eqn{G_n(x_1), \ldots, G_n(x_m)} are estimated
#' distribution function of the observed variable evaluated at the distinct observed
#' event times.
#'
#' @references
#' Lynden-Bell, D. (1971). A method of allowing for known observational selection in small
#' samples applied to 3CR quasars. Monthly Notices of the Royal Astronomical Society, 155:95-118.
#'
#' Woodroofe, M. (1985). Estimating a distribution function with truncated data. The
#' Annals of Statistics, 13: 163-177.
#'
#' @examples
#' aa <- lt.sim(500,rnorm,trunc.args=list(4,1),rnorm,target.args = list(4,1))
#' g.est <- lt.truncation.estimator(o.sample=aa)
#' plot(g.est$trunc.time,g.est$estimation,type='s',lwd=2,col=2)
#' curve(1-pnorm(x,4,1),add=TRUE,lwd=2,
#' xlab='t',ylab='Estimated distribution function')
#' # nontruncation probability estimator:
#' g.est$gamma
#'
#' @export

lt.truncation.estimator <- function(o.sample,
                                    holes=c('holes','add','conditional')
                                    ){

  holes <- match.arg(holes)

  o.sample2 <- data.frame(u2=-o.sample[,2],x2=-o.sample[,1])
  est2 <- lt.estimator(o.sample2, conf.int=FALSE, holes=holes)
  trunc.time <- -est2$fail.time
  estimation <- 1-est2$estimation

  est <- lt.estimator(o.sample,conf.int=FALSE,holes=holes)

  ind.x <- findInterval(est$fail.time,rev(trunc.time) )
  g.hat.xi <- rev(estimation)[ind.x]

  gamma.est <- sum(est$prob * g.hat.xi)


  lista <- list(
    trunc.time  = trunc.time,
    estimation  = estimation,
    prob        = est2$prob,
    gamma       = gamma.est
  )
  class(lista) <- "t.estimators"
  return(lista)
}
