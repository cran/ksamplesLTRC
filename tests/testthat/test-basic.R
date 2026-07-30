test_that("rc.test runs for log-rank test and returns expected structure", {
  data <- data.frame(
    y = c(1, 2, 3, 4, 5, 6),
    delta = c(1, 1, 0, 1, 0, 1),
    group = c(1, 1, 1, 2, 2, 2)
  )

  out <- rc.test(
    p.sample = data,
    weights = c(0.5, 0.5),
    tests = "logrank",
    B = 0,
    plot.curves = FALSE,
    p = 0,
    q = 0,
    seed = 123
  )

  expect_s3_class(out, "ksamples_test")
  expect_named(
    out,
    c("statistic", "p.value", "logrank", "B", "tests", "curves", "groups"),
    ignore.order = TRUE
  )
  expect_equal(out$B, 0)
  expect_equal(out$tests, "logrank")
  expect_false(is.null(out$logrank))
})
