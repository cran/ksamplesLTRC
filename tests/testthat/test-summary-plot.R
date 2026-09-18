test_that("summary.ksamples_test works correctly on a ltrc.test() object", {
  df <- ltrc.sim(n = 60, rtrunc = rnorm, trunc.args = list(3,1),
                 rtarget = rnorm, target.args = list(5,1),
                 rcens = rexp, cens.args = list(1), seed = 1)
  df <- data.frame(df, group = rep(1:2, 30))
  out <- ltrc.test(df, tests = c("ks","cvm","logrank"), B = 0, p = 0, q = 0)

  expect_s3_class(out, "ksamples_test")
  s <- summary(out)
  expect_s3_class(s, "data.frame")
  expect_named(s, c("Test", "Statistic", "p.value"))
  expect_equal(nrow(s), 3)
  expect_no_error(plot(out))
})


test_that("plot.ksamples_test and lines.ksamples_test runs without error
          on an estimator object", {
  df <- ltrc.sim(n = 30, rtrunc = rnorm, trunc.args = list(3,1),
                 rtarget = rnorm, target.args = list(5,1),
                 rcens = rexp, cens.args = list(1), seed = 1)
  out <- ltrc.estimator(df, conf.int = TRUE)
  df2 <- ltrc.sim(n = 30, rtrunc = rnorm, trunc.args = list(3,1),
                 rtarget = rnorm, target.args = list(5,1),
                 rcens = rexp, cens.args = list(1), seed = 1)
  out2 <- ltrc.estimator(df2, conf.int = TRUE)

  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(plot(out))
  expect_no_error(lines(out2))
})

test_that("plot.t.ksamples_test and lines.t.ksamples_test runs without error
          on an estimator object", {
            df <- ltrc.sim(n = 30, rtrunc = rnorm, trunc.args = list(3,1),
                           rtarget = rnorm, target.args = list(5,1),
                           rcens = rexp, cens.args = list(1), seed = 1)
            out <- ltrc.truncation.estimator(df)
            df2 <- ltrc.sim(n = 30, rtrunc = rnorm, trunc.args = list(3,1),
                            rtarget = rnorm, target.args = list(5,1),
                            rcens = rexp, cens.args = list(1), seed = 1)
            out2 <- ltrc.truncation.estimator(df2)

            pdf(NULL)
            on.exit(dev.off())
            expect_no_error(plot(out))
            expect_no_error(lines(out2))
})

test_that("plot.ksamples_test runs without error on a test object with bootstrap curves kept", {
  df <- ltrc.sim(n = 60, rtrunc = rnorm, trunc.args = list(3,1),
                 rtarget = rnorm, target.args = list(5,1),
                 rcens = rexp, cens.args = list(1), seed = 1)
  df <- data.frame(df, group = rep(1:2, 30))
  out <- ltrc.test(df, tests = "ks", B = 20, p = 0, q = 0,
                   plot.curves = TRUE, keep.boot = TRUE, seed = 1)

  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(plot(out))
})
