## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(KMsurv)
library(ksamplesLTRC)

## ----quick-example, fig.width=7, fig.height=5---------------------------------
ltrc.data <- data.frame(
  ltrc.sim(n = 500, rtrunc = rnorm, trunc.args = list(4, 1),
           rtarget = rnorm, target.args = list(4, 1),
           rcens = rexp, cens.args = list(1), seed = 1212),
  group = rep(1:2, each = 250)
)

test <- ltrc.test(ltrc.data, tests = "ks", weights = rep(1/2, 2), B = 200, plot.curves = TRUE)
test$statistic
test$p.value

## ----simulation---------------------------------------------------------------
ltrc.data <- ksamplesLTRC::ltrc.sim(n = 500, rtrunc = rnorm, trunc.args = list(4, 1),
                 rtarget = rnorm, target.args = list(4, 1),
                 rcens = rexp, cens.args = list(1), seed = 12)
head(ltrc.data)
mean(ltrc.data[, 1] <= ltrc.data[, 2])
table(ltrc.data$delta)

## ----estimation, fig.width=7, fig.height=5------------------------------------
ltrc.data <- ksamplesLTRC::ltrc.sim(1000, rnorm, trunc.args = list(4, 1),
                                     rnorm, target.args = list(4, 1),
                                     rexp, cens.args = list(1))

est <- ltrc.estimator(o.sample = ltrc.data, conf.int = TRUE, conf.type = 'loglog')
plot(est, col = 'forestgreen', xlab = 't', ylab = 'Estimated survival probability')
curve(1 - pnorm(x, 4, 1), add = TRUE)
str(est)

## ----holes, fig.width=7, fig.height=5-----------------------------------------
data("channing", package = "KMsurv")
o.sample <- data.frame(entry = channing$ageentry, time = channing$age,
                        status = channing$death, group = channing$gender)

o.sample <- subset(o.sample, group == 1)
head(o.sample)
plot(ksamplesLTRC::ltrc.estimator(o.sample, holes = 'holes'), col = 2,
     xlab = 't', ylab = 'Estimated survival probability')
lines(ksamplesLTRC::ltrc.estimator(o.sample, holes = 'add'), col = 3)
lines(ksamplesLTRC::ltrc.estimator(o.sample, holes = 'conditional'), col = 4)

## ----channing-ltrc-test, fig.width=7, fig.height=5----------------------------
data("channing", package = "KMsurv")
data.channing <- data.frame(entry = channing$ageentry, time = channing$age,
                             status = channing$death, group = channing$gender)
aa <- ksamplesLTRC::ltrc.test(data.channing, weights = rep(1/2, 2), holes = 'conditional',
                               tests = c('ks', 'cvm', 'logrank'), plot.curves = TRUE,
                               p = c(0, 1, 0.5, 0), q = c(0, 0, 0.5, 1)) # equal weights for two groups
aa$p.value
aa$logrank

## ----channing-ksample-test, fig.width=7, fig.height=5-------------------------
data("channing", package = "KMsurv")
bb <- ksamplesLTRC::ksample.test(time = channing$age, entry = channing$ageentry,
                                  status = channing$death, group = channing$gender,
                                  weights = rep(1/2, 2), B = 500, holes = 'conditional',
                                  tests = c('ks', 'cvm', 'logrank'), plot.curves = TRUE,
                                  p = c(0, 1, 0.5, 0), q = c(0, 0, 0.5, 1)) # equal weights for two groups
bb$p.value
bb$logrank

## ----single-test-shortcut-----------------------------------------------------
ksamplesLTRC::ltrc.pv.ks(data.channing, weights = rep(1/2, 2), B = 200, holes = 'conditional')$pvalue

