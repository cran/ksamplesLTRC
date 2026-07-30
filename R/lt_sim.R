#' Left-truncated data simulation
#'
#' Simulates data under left truncation from prespecified variables for the
#' truncation and target variables.
#'
#' @param n Sample size
#' @param rtrunc A function that simulates data from a distribution
#' @param trunc.args A list with the parameters of the truncation variable
#' @param rtarget A function that simulates data from a distribution
#' @param target.args A list with the parameters of the target distribution
#' @param seed Optional seed.
#'
#'
#' @return
#' A data frame with \eqn{n} rows and two columns: the truncation times are in the first column,
#' and the event times, in the second one.
#'
#' @examples
#' aa <- lt.sim(1000,rnorm,trunc.args=list(4,1),rnorm,target.args = list(4,1))
#' head(aa)
#' plot(lt.estimator(aa),col=2,xlab='t',ylab='Estimated survival probiblity')
#' curve(1-pnorm(x,4,1),add=TRUE)
#'
#' @export


lt.sim <- function(n,
                      rtrunc, trunc.args = list(),
                      rtarget, target.args = list(),
                      seed = NULL) {

  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n <= 0 || n != as.integer(n)) {
    stop("'n' must be one positive integer.")
  }
  n <- as.integer(n)

  if (!is.function(rtrunc)) {
    stop("'rtrunc' must be a function generating truncation times.")
  }

  if (!is.function(rtarget)) {
    stop("'rtarget' must be a function generating target/survival times.")
  }

  if (!is.list(trunc.args)) {
    stop("'trunc.args' must be a list.")
  }

  if (!is.list(target.args)) {
    stop("'target.args' must be a list.")
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  draw_from <- function(rfun, size, args, name) {
    out <- tryCatch(
      do.call(rfun, c(list(n = size), args)),
      error = function(e) {
        stop(sprintf("Error while calling '%s': %s", name, e$message), call. = FALSE)
      }
    )

    if (!is.numeric(out)) {
      stop(sprintf("'%s' must return a numeric vector.", name))
    }

    if (length(out) != size) {
      stop(sprintf("'%s' must return a numeric vector of length 'n'.", name))
    }

    if (any(!is.finite(out))) {
      stop(sprintf("'%s' returned non-finite values.", name))
    }

    out
  }

  batch.size <- max(1000L, n)
  max.draws <- 1e7

  u.keep <- numeric(n)
  x.keep <- numeric(n)

  n.accepted <- 0L
  n.generated <- 0L

  while (n.accepted < n) {

    m <- min(batch.size, max.draws - n.generated)

    if (m <= 0L) {
      stop(
        paste(
          "Maximum number of generated pairs reached before obtaining 'n' observable pairs.",
          "Try using distributions with a larger probability of satisfying U <= X."
        ),
        call. = FALSE
      )
    }

    u.new <- draw_from(rtrunc, size = m, args = trunc.args, name = "rtrunc")
    x.new <- draw_from(rtarget, size = m, args = target.args, name = "rtarget")

    ok <- (u.new <= x.new)
    n.ok <- sum(ok)

    if (n.ok > 0L) {
      take <- min(n - n.accepted, n.ok)
      idx <- which(ok)[seq_len(take)]

      rng <- (n.accepted + 1L):(n.accepted + take)
      u.keep[rng] <- u.new[idx]
      x.keep[rng] <- x.new[idx]

      n.accepted <- n.accepted + take
    }

    n.generated <- n.generated + m
  }

  out <- data.frame(
    u = u.keep,
    x = x.keep
  )

  attr(out, "n.generated") <- n.generated
  attr(out, "acceptance.rate") <- n / n.generated

  return(out)
}
