# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

asim <- function(n = 1000, seed = 100,
                 param_list = list(
                   t.max = 12.48,
                   knot.longi = c(0, 1.36, 2.96, 12.48),
                   knot.surv = c(0, 0, 0, 0, 1.05, 1.82, 2.41, 3.13, 3.91, 4.67, 5.56, 6.93, 12.48, 12.48, 12.48, 12.48),
                   betas = c(2.3428173, 0.2790604, 0.6051217, 0.9545555, 0.0162572),
                   sigma.y = 0.1452483,
                   gammas = c(0.4987653, 0.2340489),
                   alpha = matrix(c(0.1275125, 0.4207432, 3.010636, 2.621581), 2, 2),
                   gambh = matrix(c(-6.775537, -4.716391, -2.835004, -1.649737, -1.541211, -1.791381, -1.846360, -1.746593, -1.854139, -2.041937, -2.182281, -2.315315,
                                    -5.755751, -4.994471, -4.431788, -4.261764, -4.356832, -4.474064, -4.600834, -4.686764, -4.778267, -4.917803, -5.076011, -5.209827), 12, 2),
                   mean.Cens = 5.114849,
                   D_c3 = matrix(c(0.48387474, -0.04027017, -0.07481436,  0.01722258,
                                   -0.04027017, 0.77008166,  0.46077603, -0.04049686,
                                   -0.07481436,  0.46077603,  1.37021094,  1.36480027,
                                   0.01722258, -0.04049686,  1.36480027,  2.53766219), 4, 4),
                   age.mean = -0.1248499,
                   age.sd = 6.86713,
                   density.mean = -2.272663,
                   density.sd = 0.6042731,
                   cvisit.sep = 0.25,
                   cvisit.sd = 0.036
                 ), Pcomp = 1, Bcomp = 1, Bsen = 1) {

  # Preparation
  f.org <- JMbayes:::gaussKronrod()$sk
  w <- JMbayes:::gaussKronrod()$wk
  knot.longi <- param_list$knot.longi
  knot.surv <- param_list$knot.surv

  set.seed(seed)

  # Parameters
  t.max <- param_list$t.max
  betas <- param_list$betas
  sigma.y <- param_list$sigma.y
  gammas <- param_list$gammas
  alpha <- param_list$alpha
  gambh <- param_list$gambh
  mean.Cens <- param_list$mean.Cens
  V <- param_list$D_c3
  V <- nearPD(V)$mat
  D <- V

  # Covariates
  cvisit.sep <- param_list$cvisit.sep
  cvisit.sd <- param_list$cvisit.sd
  tm <- t.max %/% cvisit.sep
  times <- c(replicate(n,
                       c(0, seq(cvisit.sep, tm * cvisit.sep, cvisit.sep) + rtruncnorm(tm, -cvisit.sep/2, cvisit.sep/2, 0, cvisit.sd),
                         t.max)))
  density <- rnorm(n, param_list$density.mean, param_list$density.sd)
  age <- rnorm(n, param_list$age.mean, param_list$age.sd)
  K <- length(times)/n
  DF <- data.frame(TimeSince_Dx = times, density = rep(density, each = K), DxAge = rep(age, each = K))
  X <- model.matrix(~ ns(TimeSince_Dx, knots = knot.longi[2:3], B = knot.longi[c(1,4)]) + DxAge  , data = DF)
  Z <- model.matrix(~ ns(TimeSince_Dx, knots = knot.longi[2:3], B = knot.longi[c(1,4)]), data = DF)


  ## Outcome generation
  b <- mvrnorm(n, rep(0, nrow(D)), D)

  id <- rep(1:n, each = K)
  eta.y <- as.vector(X %*% betas + rowSums(Z * b[id, 1:4])) # linear predictor

  y <- rgt(n * K, mu =  eta.y, sigma = sigma.y, df =  3)

  eta.t <- cbind(as.vector(density * gammas[1]),
                 as.vector(density * gammas[2]))

  invS <- function (t, u, i) {
    h <- function (s, k) {
      age.i <- age[i]
      XX <- cbind(1, ns(s, knots = knot.longi[2:3], B = knot.longi[c(1,4)]), age.i)
      XX.back <- cbind(1, ns(s - 1, knots = knot.longi[2:3], B = knot.longi[c(1,4)]), age.i)
      ZZ <- cbind(1, ns(s, knots = knot.longi[2:3], B = knot.longi[c(1,4)]))
      ZZ.back <- cbind(1, ns(s - 1, knots = knot.longi[2:3], B = knot.longi[c(1,4)]))
      f <- as.vector(XX %*% betas + rowSums(ZZ * b[rep(i, nrow(ZZ)), 1:4]))
      f.back <- as.vector(XX.back %*% betas + rowSums(ZZ.back * b[rep(i, nrow(ZZ)), 1:4]))
      bh <- splineDesign(knots = knot.surv, s, ord = 4L, outer.ok = T)

      exp(bh %*% gambh[,k]  + eta.t[i,k] + f * alpha[k,1] + (f - f.back) * alpha[k,2])

    }
    integrate(h, lower = 0, upper = t, k = 1)$value + integrate(h, lower = 0, upper = t, k = 2)$value + log(u)
  }
  u <- runif(n)
  trueTimes <- matrix(NA, n)

  for (i in 1:n) {
    Up <- 50
    tries <- 5
    Root <- try(uniroot(invS, interval = c(1e-05, Up), u = u[i], i = i)$root, TRUE)
    while(inherits(Root, "try-error") && tries > 0) {
      tries <- tries - 1
      Up <- Up + 200
      Root <- try(uniroot(invS, interval = c(1e-05, Up), u = u[i], i = i)$root, TRUE)
    }
    trueTimes[i] <- if (!inherits(Root, "try-error")) Root else NA
  }
  trueTimes[is.na(trueTimes)] <- 1e06

  ## simulate which event is
  # calculate the instantaneous hazard
  h <- function (s, k ,i) {
    age.i <- age[i]
    XX <- cbind(1, ns(s, knots = knot.longi[2:3], B = knot.longi[c(1,4)]), age.i)
    XX.back <- cbind(1, ns(s - 1, knots = knot.longi[2:3], B = knot.longi[c(1,4)]), age.i)
    ZZ <- cbind(1, ns(s, knots = knot.longi[2:3], B = knot.longi[c(1,4)]))
    ZZ.back <- cbind(1, ns(s - 1, knots = knot.longi[2:3], B = knot.longi[c(1,4)]))
    f <- as.vector(XX %*% betas + rowSums(ZZ * b[rep(i, nrow(ZZ)), 1:4]))
    f.back <- as.vector(XX.back %*% betas + rowSums(ZZ.back * b[rep(i, nrow(ZZ)), 1:4]))
    bh <- splineDesign(knots = knot.surv, s, ord = 4L, outer.ok = T)

    return(ifelse(exp(bh %*% gambh[,k]  + eta.t[i,k] + f * alpha[k,1] + (f - f.back) * alpha[k,2]) < 1e-16,
                  1e-16,
                  exp(bh %*% gambh[,k]  + eta.t[i,k] + f * alpha[k,1] + (f - f.back) * alpha[k,2])))

  }

  event.2 <- sapply(1:n, function(i) {
    rbinom(1, 1, h(trueTimes[i], 2, i)/(h(trueTimes[i], 2, i) + h(trueTimes[i], 1, i))) + 1  # in rbinom 1 = event 2
  })

  # simulate censoring times from an exponential distribution,
  # and calculate the observed event times, i.e., min(true event times, censoring times)
  Ctimes <- runif(n, 0, 2 * mean.Cens)
  Time <- pmin(Ctimes, trueTimes)

  event <- sapply(1:n, function(i) {ifelse(trueTimes[i] < Ctimes[i], event.2[i], 0)}) # event indicator

  fixed_visits <- matrix(times, ncol = n)
  fv_idx <- c(1,3,5, seq(9,K, 8))
  if (K %in% fv_idx) {
    fixed_visits <- fixed_visits[fv_idx,]
  } else {
    fixed_visits <- fixed_visits[c(fv_idx, K),]
  }

  Time1 <- sapply(1:n, function(i) {
    max(fixed_visits[fixed_visits[,i] <= Time[i],i])
  })

  Time2 <- sapply(1:n, function(i) {
    ifelse(event[i] == 1, min(fixed_visits[fixed_visits[,i] >= Time[i], i]), Time[i])
  })

  # sum(Time2 < Time1)
  # sum(Time1 <0 & Time2 <0)

  times.mat <- matrix(times, ncol = n)
  ind <- sapply(1:n, function(i) {
    if(event[i] != 1) {
      times.mat[,i] <= rep(Time[i], K)
    } else {
      times.mat[,i] <= rep(Time2[i], K)
    }
  })

  y.obs <- y[ind]

  id.obs <- id[ind]
  id.obs <- match(id.obs, unique(id.obs))

  dat <- DF[ind, ]
  dat$CISNET_ID <- id.obs
  dat$PSAValue <- y.obs

  dat$time.cmp1 <- Time1[id.obs]
  dat$time.cmp2 <- Time2[id.obs]
  dat$time <- Time[id.obs]
  dat$status.cmp <- event[id.obs]
  dat$b1 <- b[id.obs,1]
  dat$b2 <- b[id.obs,2]
  dat$b3 <- b[id.obs,3]
  dat$b4 <- b[id.obs,4]

  dat <- dat[c("CISNET_ID", "TimeSince_Dx", "PSAValue", "time", "time.cmp1",
               "time.cmp2","status.cmp", "density","DxAge",
               "b1", "b2", "b3", "b4")]

  dat.id <- dat[!duplicated(dat$CISNET_ID), ]

  return(list(dat = dat,
              dat.id = dat.id))
}
