test_that("ltrc.sim and lt.sim always produces u <= y", {
  out <- ltrc.sim(n = 1000, rtrunc = rnorm, trunc.args = list(3,1),
                  rtarget = rnorm, target.args = list(5,1),
                  rcens = rexp, cens.args = list(1), seed = 1)
  out2 <- lt.sim(n = 1000, rtrunc = rnorm, trunc.args = list(3,1),
                  rtarget = rnorm, target.args = list(5,1),
                  seed = 1)
  expect_true(all(out$u <= out$y))
  expect_true(all(out2$u <= out2$y))
  expect_named(out, c("u","y","delta"))
})

test_that("ltrc.sim, lt.sim and rc.sim return exactly n rows", {
  out <- ltrc.sim(n = 50, rtrunc = rnorm, trunc.args = list(3,1),
                  rtarget = rnorm, target.args = list(5,1),
                  rcens = rexp, cens.args = list(1), seed = 1)
  out2 <- lt.sim(n = 50, rtrunc = rnorm, trunc.args = list(3,1),
                  rtarget = rnorm, target.args = list(5,1),
                  seed = 1)
  out3 <- rc.sim(n = 50,
                 rtarget = rnorm, target.args = list(5,1),
                 rcens = rexp, cens.args = list(1), seed = 1)
  expect_equal(nrow(out), 50)
  expect_equal(nrow(out2), 50)
  expect_equal(nrow(out3), 50)
})

test_that("ltrc.sim's and rc.sim's delta are always exactly 0 or 1", {
  out <- ltrc.sim(n = 100, rtrunc = rnorm, trunc.args = list(3,1),
                  rtarget = rnorm, target.args = list(5,1),
                  rcens = rexp, cens.args = list(1), seed = 1)
  out2 <- rc.sim(n = 100,
                  rtarget = rnorm, target.args = list(5,1),
                  rcens = rexp, cens.args = list(1), seed = 1)
  expect_true(all(out$delta %in% c(0, 1)))
  expect_true(all(out2$delta %in% c(0, 1)))
})


test_that("ltrc.sim, lt.sim and rc.sim are exactly reproducible
          given the same seed", {
  args <- list(n = 30, rtrunc = rnorm, trunc.args = list(3,1),
               rtarget = rnorm, target.args = list(5,1),
               rcens = rexp, cens.args = list(1), seed = 42)
  out1 <- do.call(ltrc.sim, args)
  out2 <- do.call(ltrc.sim, args)
  expect_identical(out1, out2)

  args <- list(n = 30, rtrunc = rnorm, trunc.args = list(3,1),
               rtarget = rnorm, target.args = list(5,1),
               seed = 32)
  out1 <- do.call(lt.sim, args)
  out2 <- do.call(lt.sim, args)
  expect_identical(out1, out2)

  args <- list(n = 30,
               rtarget = rnorm, target.args = list(5,1),
               rcens = rexp, cens.args = list(1), seed = 52)
  out1 <- do.call(rc.sim, args)
  out2 <- do.call(rc.sim, args)
  expect_identical(out1, out2)
})


