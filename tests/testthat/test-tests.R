test_that("log-rank matches survival::survdiff under right censoring only", {
  set.seed(10)
  y1 <- rexp(20, 1);   c1 <- rexp(20, 0.75)
  y2 <- rexp(20, 1.3); c2 <- rexp(20, 0.5)

  df <- data.frame(
    y     = c(pmin(y1, c1), pmin(y2, c2)),
    delta = c(as.integer(y1 <= c1), as.integer(y2 <= c2)),
    group = rep(c(1, 2), each = 20)
  )

  reg  <- survival::survdiff(survival::Surv(y, delta) ~ group, data = df)
  out <- rc.test(df, tests = "logrank", p = 0, q = 0, B = 0)

  expect_equal(out$logrank$statistic, reg$chisq, tolerance = 1e-6)
})

test_that("log-rank matches survival::coxph under left truncation only", {
  set.seed(10)
  u1 <- rexp(20, 1);   x1 <- rexp(20, 0.8)
  u2 <- rexp(20, 1.3); x2 <- rexp(20, 1)

  uu1 <- u1[u1 <= x1 ]; xx1 <- x1[ u1 <= x1 ]
  uu2 <- u2[u2 <= x2 ]; xx2 <- x2[ u2 <= x2 ]

  df <- data.frame(
    u = c(uu1, uu2),
    y = c(xx1, xx2),
    delta = rep(1,length(c(xx1,xx2))),
    group = rep(c(1, 2), c( length(xx1),length(xx2)) )
  )

  reg  <- survival::coxph(survival::Surv(u,y, delta) ~ group, data = df)
  out <- lt.test(df[,-3], tests = "logrank", p = 0, q = 0, B = 0)

  expect_equal(out$logrank$statistic, reg$score, tolerance = 1e-6)
})

test_that("log-rank matches survival::coxph under left truncation
          and right censoring", {
  set.seed(22)
  u1 <- rexp(20, 1);   x1 <- rexp(20, 0.8)
  u2 <- rexp(20, 1.3); x2 <- rexp(20, 1)

  uu1 <- u1[u1 <= x1 ]; xx1 <- x1[ u1 <= x1 ]
  uu2 <- u2[u2 <= x2 ]; xx2 <- x2[ u2 <= x2 ]

  rc1 <- rexp(length(xx1),0.5); rc2 <- rexp(length(xx2),1)

  df <- data.frame(
    u = c(uu1, uu2),
    y = c( pmin( xx1,uu1+rc1 ), pmin( xx2, uu2+rc2 ) ),
    delta = c( as.integer( xx1 <= uu1+rc1 ), as.integer( xx2 <= uu2+rc2 ) ),
    group = rep(c(1, 2), c( length(xx1),length(xx2)) )
  )

  reg  <- survival::coxph(survival::Surv(u,y, delta) ~ group, data = df)
  out <- ltrc.test(df, tests = "logrank", p = 0, q = 0, B = 0)

  expect_equal(out$logrank$statistic, reg$score, tolerance = 1e-6)
})

test_that("KS, CvM and log-rank are exactly zero (and p-values exactly 1) for
          two identical groups", {
  set.seed(30)
  y <- rexp(15, 1); delta <- rbinom(15, 1, 0.7)
  df <- rbind(
    data.frame(y = y, delta = delta, group = 1),
    data.frame(y = y, delta = delta, group = 2)
  )

  out <- rc.test(df, tests = c("ks", "cvm", "logrank"), p = 0, q = 0, B = 10, seed = 1)

  expect_equal(out$statistic$ks,  0, tolerance = 1e-8)
  expect_equal(out$statistic$cvm, 0, tolerance = 1e-8)
  expect_equal(out$logrank$statistic, 0, tolerance = 1e-8)
  expect_equal(out$p.value$ks,  1, tolerance = 1e-8)
  expect_equal(out$p.value$cvm, 1, tolerance = 1e-8)
})

test_that("The statistics under truncation and censoring reduce to the ones under right
          censoring when no truncation is present",{
      set.seed(18)
      x <- rnorm(50,5,1)
      c <- rnorm(50,5.2,1)

      df <- data.frame(
        u = rep(0,50),
        y = pmin(x,c),
        delta = as.integer( x <= c ),
        group = rep(1:2, each=25)
      )

      test.rc <- rc.test(p.sample = df[,-1],
                         tests = c("ks","cvm","logrank"),
                         B = 0, p = 0, q = 0)
      test.ltrc <- ltrc.test(p.sample = df,
                             tests = c("ks","cvm","logrank"),
                             B = 0, p = 0, q = 0)
      expect_equal(test.rc$statistic$ks, test.ltrc$statistic$ks, tolerance = 1e-6)
      expect_equal(test.rc$statistic$cvm, test.ltrc$statistic$cvm, tolerance = 1e-6)
      expect_equal(test.rc$logrank$statistic, test.ltrc$logrank$statistic, tolerance = 1e-6)
})

test_that("The statistics under truncation and censoring reduce to the ones under left
          truncation when no censoring is present",{
      set.seed(50)
      x <- rnorm(60,3,1)
      u <- rnorm(60,2.7,1)

      df <- data.frame(
            u = u[u <= x],
            y = x[u <= x],
            delta = rep(1, sum(u <= x) ),
            group = rep(1:2, c(floor( sum(u <= x)/2 ), ceiling( sum(u <= x)/2 ) ) )
      )
      test.rc <- lt.test(p.sample = df[,-3],
                         tests = c("ks","cvm","logrank"),
                         B = 0, p = 0, q = 0)
      test.ltrc <- ltrc.test(p.sample = df,
                             tests = c("ks","cvm","logrank"),
                             B = 0, p = 0, q = 0)
      expect_equal(test.rc$statistic$ks, test.ltrc$statistic$ks, tolerance = 1e-6)
      expect_equal(test.rc$statistic$cvm, test.ltrc$statistic$cvm, tolerance = 1e-6)
      expect_equal(test.rc$logrank$statistic, test.ltrc$logrank$statistic, tolerance = 1e-6)
})

test_that("The Kolmogorov--Smirnov and Cramér--von Mises statistics
          for left-truncated and right-censored data are proportional to the classical
          implementations when no truncation nor censoring are present",{
        set.seed(73)
        x <- rnorm(50,5,1)
        weights <- c(1/3, 2/3)
        n <- c(30,20)

        df <- data.frame(u = rep(0,50),
                         y = x,
                         delta = rep(1,50),
                         group = rep(1:2, n)
        )
        ks.ord <- ks.test(df$y[df$group == 1], df$y[df$group == 2])
        cvm.ord <- twosamples::cvm_stat(df$y[df$group == 1], df$y[df$group == 2], power = 2)
        ks.ltrc <- ltrc.ks(p.sample = df, weights = weights)
        cvm.ltrc <- ltrc.cvm(p.sample = df, weights = n/sum(n))

        expect_equal(ks.ltrc$statistic, as.numeric( sum(n*rev(weights)^2)*ks.ord$statistic^2 ), tolerance = 1e-6)
        expect_equal(cvm.ltrc$statistic, (n[1] * n[2] / sum(n)^2) * cvm.ord, tolerance = 1e-6 )
})

test_that("weighted log-rank matches survival::survdiff (rho family, q=0) under right censoring", {
  set.seed(1)
  y1 <- rexp(20,1); c1 <- rexp(20,0.7)
  y2 <- rexp(20,1.3); c2 <- rexp(20,0.7)
  df <- data.frame(
    y     = c(pmin(y1,c1), pmin(y2,c2)),
    delta = c(as.integer(y1<=c1), as.integer(y2<=c2)),
    group = rep(1:2, each = 20)
  )
  out.survival  <- survival::survdiff(survival::Surv(y, delta) ~ group, data = df, rho = 1)
  out <- rc.logrank(df, p = 1, q = 0)

  expect_equal(out[1,1], out.survival$chisq, tolerance = 1e-6)
})

test_that("weighted log-rank is exactly zero for identical groups, across the full (p,q) family", {
  set.seed(30)
  y <- rexp(15, 1); delta <- rbinom(15, 1, 0.7)
  df <- rbind(data.frame(y=y, delta=delta, group=1),
              data.frame(y=y, delta=delta, group=2))

  for (pq in list(c(0,1), c(0.5,0.5), c(2,0), c(1,2))) {
    out <- rc.test(df, tests = "logrank", p = pq[1], q = pq[2], B = 0)
    expect_equal(out$logrank$statistic, 0, tolerance = 1e-8)
  }
})

test_that("KS and CvM are exactly zero for three identical groups", {
  set.seed(5)
  y <- rexp(15, 1); delta <- rbinom(15, 1, 0.7)
  df <- rbind(data.frame(y=y, delta=delta, group=1),
              data.frame(y=y, delta=delta, group=2),
              data.frame(y=y, delta=delta, group=3))
  out <- rc.test(df, tests = c("ks","cvm"), B = 0)
  expect_equal(out$statistic$ks,  0, tolerance = 1e-8)
  expect_equal(out$statistic$cvm, 0, tolerance = 1e-8)
})

test_that("'conditional' CvM matches re-computing on the restricted subsample", {
  df <- data.frame(
    u     = c(0,3,3,3, 0,0,0,0),
    y     = c(2,8,8,6, 3,5,7,10),
    delta = c(1,1,1,1, 1,1,1,0),
    group = rep(1:2, each = 4)
  )
  out_cond <- ltrc.cvm(df, holes = "conditional")

  restricted <- df[df$y > 2, ]
  out_manual <- ltrc.cvm(restricted, holes = "holes")
  expect_equal(out_cond$statistic, out_manual$statistic, tolerance = 1e-8)
})

test_that("ltrc.cvm rejects weights that don't sum to 1", {
  df <- data.frame(u=0, y=1:10, delta=1, group=rep(1:2,5))
  expect_error(ltrc.cvm(df, weights = c(0.3, 0.3)))
})

test_that("ltrc.cvm rejects U > Y", {
  df <- data.frame(u=c(5,0), y=c(1,2), delta=1, group=1:2)
  expect_error(ltrc.cvm(df))
})

test_that("ltrc.pv.ks/.cvm give the same statistic as ltrc.ks/.cvm", {
  df <- ltrc.sim(n=100, rtrunc=rnorm, trunc.args=list(3,1),
                 rtarget=rnorm, target.args=list(5,1),
                 rcens=rexp, cens.args=list(1), seed=1)
  df <- data.frame(df, group = rep(1:2, 50))
  expect_equal(ltrc.pv.ks(df, B=1)$statistic,  ltrc.ks(df)$statistic,  tolerance = 1e-10)
  expect_equal(ltrc.pv.cvm(df, B=1)$statistic, ltrc.cvm(df)$statistic, tolerance = 1e-10)
})

test_that("The test statistics in the left truncation and right censoring functions
          are equal",{
          df <- data.frame( ltrc.sim(n=100, rtrunc=rnorm, trunc.args=list(3,1),
                                     rtarget=rnorm, target.args=list(4,1),
                                     rcens=rexp, cens.args=list(1), seed=84),
                             group = rep(1:2, each = 50)
                          )
          bb <- ltrc.test(df, tests = c("ks","cvm"), B = 1)
          cc <- ksample.test(time = df$y, entry = df$u, status = df$delta, group = df$group,
                             tests = c('ks','cvm'), B = 1)
          expect_equal(ltrc.ks(p.sample = df)$statistic, ltrc.pv.ks(df, B = 1)$statistic)
          expect_equal(ltrc.ks(p.sample = df)$statistic, bb$statistic$ks)
          expect_equal(ltrc.ks(p.sample = df)$statistic, cc$statistic$ks)
          expect_equal(ltrc.cvm(p.sample = df)$statistic, ltrc.pv.cvm(df, B = 1)$statistic)
          expect_equal(ltrc.cvm(p.sample = df)$statistic, bb$statistic$cvm)
          expect_equal(ltrc.cvm(p.sample = df)$statistic, cc$statistic$cvm)
})

test_that("The test statistics in the left truncation functions
          are equal",{
            df <- data.frame( lt.sim(n=100, rtrunc=rnorm, trunc.args=list(3,1),
                                       rtarget=rnorm, target.args=list(4,1),
                                       seed=84),
                              group = rep(1:2, each = 50)
            )
            bb <- lt.test(df, tests = c("ks","cvm"), B = 1)
            cc <- ksample.test(time = df$x, entry = df$u, group = df$group,
                               tests = c('ks','cvm'), B = 1)
            expect_equal(lt.ks(p.sample = df)$statistic, lt.pv.ks(df, B = 1)$statistic)
            expect_equal(lt.ks(p.sample = df)$statistic, bb$statistic$ks)
            expect_equal(lt.ks(p.sample = df)$statistic, cc$statistic$ks)
            expect_equal(lt.cvm(p.sample = df)$statistic, lt.pv.cvm(df, B = 1)$statistic)
            expect_equal(lt.cvm(p.sample = df)$statistic, bb$statistic$cvm)
            expect_equal(lt.cvm(p.sample = df)$statistic, cc$statistic$cvm)
})

test_that("The test statistics in the right censoring functions are equal",{
            df <- data.frame( rc.sim(n=100,
                                       rtarget=rweibull, target.args=list(2,3),
                                       rcens=rweibull, cens.args=list(2.5,3), seed=84),
                              group = rep(1:2, each = 50)
            )
            bb <- rc.test(df, tests = c("ks","cvm"), B = 1)
            cc <- ksample.test(time = df$y, status = df$delta, group = df$group,
                               tests = c('ks','cvm'), B = 1)
            expect_equal(rc.ks(p.sample = df)$statistic, rc.pv.ks(df, B = 1)$statistic)
            expect_equal(rc.ks(p.sample = df)$statistic, bb$statistic$ks)
            expect_equal(rc.ks(p.sample = df)$statistic, cc$statistic$ks)
            expect_equal(rc.cvm(p.sample = df)$statistic, rc.pv.cvm(df, B = 1)$statistic)
            expect_equal(rc.cvm(p.sample = df)$statistic, bb$statistic$cvm)
            expect_equal(rc.cvm(p.sample = df)$statistic, cc$statistic$cvm)
})

test_that("ksample.test dispatches to ltrc.test/lt.test/rc.test and matches them exactly", {
  set.seed(1)
  df <- ltrc.sim(n=60, rtrunc=rnorm, trunc.args=list(3,1),
                 rtarget=rnorm, target.args=list(5,1),
                 rcens=rexp, cens.args=list(1), seed=1)
  df <- data.frame(df, group = rep(1:2, 30))

  a1 <- ksample.test(time=df$y, entry=df$u, status=df$delta, group=df$group, tests="logrank", B=0)
  a2 <- ltrc.test(p.sample=data.frame(U=df$u,Y=df$y,delta=df$delta,group=df$group), tests="logrank", B=0)
  expect_equal(a1$logrank$statistic, a2$logrank$statistic)

  b1 <- ksample.test(time=df$y, entry=df$u, group=df$group, tests="logrank", B=0)
  b2 <- lt.test(p.sample=data.frame(U=df$u,Y=df$y,group=df$group), tests="logrank", B=0)
  expect_equal(b1$logrank$statistic, b2$logrank$statistic)

  c1 <- ksample.test(time=df$y, status=df$delta, group=df$group, tests="logrank", B=0)
  c2 <- rc.test(p.sample=data.frame(Y=df$y,delta=df$delta,group=df$group), tests="logrank", B=0)
  expect_equal(c1$logrank$statistic, c2$logrank$statistic)
})

test_that("ksample.test errors when neither entry nor status is given", {
  expect_error(ksample.test(time=1:10, group=rep(1:2,5)), "not currently supported")
})

test_that("ksample.test rejects entry > time and non-0/1 status", {
  expect_error(ksample.test(time=c(1,2), entry=c(5,0), group=1:2, status=c(1,1)), "entry")
  expect_error(ksample.test(time=1:10, status=rep(2,10), group=rep(1:2,5)))
})





