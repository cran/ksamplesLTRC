library(survival)

test_that("estimator matches Kaplan-Meier under no truncation", {
  set.seed(1)
  x <- rnorm(20,4,1)
  c <- rnorm(20,4,1)
  df <- data.frame(y = pmin(x,c), delta = as.integer(x <= c))

  km <- survival::survfit(survival::Surv(y, delta) ~ 1, data = df)
  km_surv <- summary(km, times = sort(df$y[df$delta == 1]) )$surv

  out <- rc.estimator(df)

  expect_equal(out$estimation, km_surv, tolerance = 1e-8)
})

test_that("estimator matches Lynden-Bell without censoring", {
  set.seed(2)
  u <- rnorm(15,3,1)
  x <- rnorm(15,3.5,1)
  uu <- u[u<= x]
  xx <- x[u <= x]
  delta <- rep(1,length(xx))

  df <- data.frame(u=uu,y=xx,delta=delta)

  lb <- survival::survfit(survival::Surv(u,y,delta)~1,data=df)

  out <- lt.estimator(df)
  lb_surv <- summary(lb, times = out$fail.time )$surv

  expect_equal(out$estimation, lb_surv, tolerance = 1e-8)
})

test_that("estimator matches survfit::survival
          under left truncation and right censoring", {
  set.seed(3)
  u <- rnorm(15,3,1)
  x <- rnorm(15,3.5,1)
  uu <- u[u<= x]
  xx <- x[u <= x]
  rc <- rexp(length(xx),0.5)

  df <- data.frame(u = uu,y = pmin(xx,uu+rc),delta = as.integer(xx <= uu + rc))

  est.survival <- survival::survfit(survival::Surv(u,y,delta) ~ 1,data = df)

  out <- ltrc.estimator(df)
  est_surv <- summary(est.survival, times = out$fail.time )$surv

  expect_equal(out$estimation, est_surv, tolerance = 1e-8)
})

test_that("ignoring a hole collapses the estimator to zero at and after it", {
  df <- data.frame(
    u     = c(0, 0, 0, 5),
    y     = c(3, 4, 5, 8)
    )
  out <- lt.estimator(df, holes = "holes")
  expect_equal(tail(out$estimation, 1), 0)
})

test_that("'add' multiplies survival by exactly 0.5 at a unitary-risk-set hole", {
  df <- data.frame(
    u     = c(0, 3, 3, 3),
    y     = c(2, 8, 8, 6),
    delta = c(1, 1, 1, 1)
  )
  out_add <- lt.estimator(df, holes = "add")
  expect_equal(out_add$estimation[1], 0.5, tolerance = 1e-8)
})

test_that("'conditional' matches re-estimating on data restricted beyond the hole", {
  df <- data.frame(
    u     = c(0, 3, 3, 3),
    y     = c(2, 8, 8, 6),
    delta = c(1, 1, 1, 1)
  )
  out_cond <- lt.estimator(df, holes = "conditional")
  df_restricted <- df[df$y > out_cond$hole, ]
  out_manual <- lt.estimator(df_restricted, holes = "holes")
  expect_equal(out_cond$estimation, out_manual$estimation, tolerance = 1e-8)
})

test_that("holes argument has no effect when sample is hole-free", {
  set.seed(42)
  df <- ltrc.sim(n = 50, rtrunc = rnorm, trunc.args = list(2,1),
                 rtarget = rnorm, target.args = list(5,1),
                 rcens = rexp, cens.args = list(1))
  out_holes <- ltrc.estimator(df, holes = "holes")
  out_add   <- ltrc.estimator(df, holes = "add")
  out_cond  <- ltrc.estimator(df, holes = "conditional")
  expect_equal(out_holes$estimation, out_add$estimation, tolerance = 1e-8)
  expect_equal(out_holes$estimation, out_cond$estimation, tolerance = 1e-8)
})

test_that("truncation estimator for left-truncated and right-censored data
          gives gamma-hat = 1 when there is no truncation", {
  df <- data.frame(u = rep(-10, 30), y = rnorm(30, 5, 1), delta = 1)
  out <- ltrc.truncation.estimator(df)
  expect_equal(out$gamma.hat, 1, tolerance = 1e-8)
})

test_that("truncation estimator for left-truncated data
          gives gamma-hat = 1 when there is no truncation", {
  df <- data.frame(u = rep(-10, 30), x = rnorm(30, 5, 1))
  out <- lt.truncation.estimator(df)
  expect_equal(out$gamma.hat, 1, tolerance = 1e-8)
})

test_that("truncation estimator for left-truncated and right-censored data
          reduces to the estimator of left-truncated data under no censoring",{
      set.seed(30)
      df <- data.frame(lt.sim(50,
                       rnorm,trunc.args=list(4,1),
                       rnorm,target.args = list(4,1)),
                       delta = rep(1,50))
      est.ltrc <- ltrc.truncation.estimator(df)
      est.lt <- lt.truncation.estimator(df)
      expect_equal(est.ltrc$trunc.time, sort(est.lt$trunc.time) )
      expect_equal(est.ltrc$estimation, rev(est.lt$estimation) )
})































