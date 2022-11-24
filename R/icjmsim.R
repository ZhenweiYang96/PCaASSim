#' A simulation function based on ICJM fitted on PASS data
#'
#' This function allows you to simulate the outcomes (censoring, progression, early treatment) of patients in AS based on PSA value and baseline age and PSA density.
#' @import MASS
#' @import Matrix
#' @import JMbayes
#' @import splines
#' @import truncnorm
#' @import mathjaxr
#' @param n number of patients on active surveillance.
#' @param seed seed value for reproducibility The default value is 100.
#' @param param_list the estimated parameters from the PASS data and the fitted ICJM: \cr
#' \itemize{
#'   \item \code{t.max} - maximum follow-up time in the PASS data
#'   \item \code{mean.Cens} - mean censoring time in the PASS data
#'   \item \code{knot.longi} - knots used in the natural cubic spline specification of the longitudinal model
#'   \item \code{knot.surv} - knots used in the P-splines of the baseline hazard in the survival model
#'   \item \code{betas} - coefficients for the fixed effects in the longitudinal submodel for PSA
#'   \item \code{sigma.y} - the standard deviation for the residuals in the longitudinal submodel for PSA
#'   \item \code{D_c3} - the covariance matrix for the random effects in the longitudinal submodel for PSA
#'   \item \code{gambh} - coefficients for the P-spline design matrix of the baseline hazard in the survival submodel
#'   \item \code{gammas} - coefficients for the exogenous covariates in the survival submodel;
#'   \item \code{alpha} - coefficients for the impact from PSA on the time-to-event outcomes in the survival submodel
#'   \item \code{age.mean} - mean age observed in the PASS data (centered around 62)
#'   \item \code{age.sd} - standard deviation of age observed in the PASS data
#'   \item \code{density.mean} - mean of the PSA density (in log) observed in the PASS data
#'   \item \code{density.sd} - standard deviation of the PSA density (in log) observed in the PASS data
#'   \item \code{cvisit.sep} - regular clinical visit interval for PSA measurement. The default is 3 months (0.25 years)
#'   \item \code{cvisit.sd} - standard deviation of variation in the clinical visit time
#' }
#' @param Pcomp the compliance rate of PSA measurements. The default value is 1.
#' @param Bcomp the compliance rate of biopsies. The default value is 1.
#' @details The `true event time` is simulated based on the event-specific hazard and bounded with the observed maximum follow-up time (12.48 yrs) in the PASS data. The `censoring time` follows a uniform distribution with the upper bound as the doubled mean of the observed censoring time in the PASS data.
#' @return two datasets. First dataset records all the longitudinal measurements (used in the longitudinal submodel), second records each subject per row (used in the survival submodel). Each datasets contains the following columns: \cr
#' \loadmathjax
#' \itemize{
#'  \item \code{CISNET_ID} - patient ID
#'  \item \code{TimeSince_Dx} - Time since start of AS
#'  \item \code{PSAValue} - transformation of PSA value (ng/ml), \mjeqn{\log_2(\text{PSA} + 1)}{ASCII representation}
#'  \item \code{time} - true event time in practice
#'  \item \code{time.cmp1} - progression-free time
#'  \item \code{time.cmp2} - treatment-free time
#'  \item \code{status.cmp} - observed event indicator (1 = cancer progression; 2 = early treatment; 0 = censoring)
#'  \item \code{density} - baseline PSA density (\mjeqn{\text{ng}/\text{ml}^2}{ASCII representation}) in log transformation, \mjeqn{\log(\text{PSA density})}{ASCII representation}
#'  \item \code{DxAge} - baseline age, centered around 62 years old
#'  \item \code{b1 - b4} - true random effects for PSA
#'  \item \code{time.prg} - true cancer progression time.
#'  \item \code{time.trt} - true early treatment time.
#'  \item \code{time.cen} - true censoring time.
#' }
#' @keywords ICJM Simulation
#' @examples
#' icjmsim(n = 1e5)
#' @export


icjmsim <- function(n = 1000, seed = 100,
                 param_list = list(
                   t.max = 12.48,
                   mean.Cens = 5.114849,
                   knot.longi = c(0, 1.36, 2.96, 12.48),
                   knot.surv = c(0, 0, 0, 0, 1.05, 1.82, 2.41, 3.13, 3.91, 4.67, 5.56, 6.93, 12.48, 12.48, 12.48, 12.48),
                   betas = c(2.3428173, 0.2790604, 0.6051217, 0.9545555, 0.0162572),
                   sigma.y = 0.1452483,
                   D_c3 = matrix(c(0.48387474, -0.04027017, -0.07481436,  0.01722258,
                                   -0.04027017, 0.77008166,  0.46077603, -0.04049686,
                                   -0.07481436,  0.46077603,  1.37021094,  1.36480027,
                                   0.01722258, -0.04049686,  1.36480027,  2.53766219), 4, 4),
                   gambh = matrix(c(-6.775537, -4.716391, -2.835004, -1.649737, -1.541211, -1.791381, -1.846360, -1.746593, -1.854139, -2.041937, -2.182281, -2.315315,
                                    -5.755751, -4.994471, -4.431788, -4.261764, -4.356832, -4.474064, -4.600834, -4.686764, -4.778267, -4.917803, -5.076011, -5.209827), 12, 2),
                   gammas = c(0.4987653, 0.2340489),
                   alpha = matrix(c(0.1275125, 0.4207432, 3.010636, 2.621581), 2, 2),
                   age.mean = -0.1248499,
                   age.sd = 6.86713,
                   density.mean = -2.272663,
                   density.sd = 0.6042731,
                   cvisit.sep = 0.25,
                   cvisit.sd = 0.036
                 ), Pcomp = 1, Bcomp = 1) {

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

  # hazard function for both events
  invS <- function (t, u, i, k) {
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
    integrate(h, lower = 0, upper = t, k = k)$value + log(u)
  }
  u1 <- runif(n) # cancer progression
  u2 <- runif(n) # early treatment
  trueTimes1 <- matrix(NA, n) # true time for progression
  trueTimes2 <- matrix(NA, n) # true time for early treatment

  for (i in 1:n) {
    # event 1 cancer progression
    Root1 <- try(uniroot(invS, interval = c(1e-05, 12.48), u = u1[i], i = i, k = 1)$root, TRUE)
    trueTimes1[i] <- if (!inherits(Root1, "try-error")) Root1 else 12.48 # bounded with the longest follow-up

    # event 2 early treatment
    Root2 <- try(uniroot(invS, interval = c(1e-05, 12.48), u = u2[i], i = i, k = 2)$root, TRUE)
    trueTimes2[i] <- if (!inherits(Root2, "try-error")) Root2 else 12.48 # bounded with the longest follow-up
  }

  # simulate censoring times from an uniform distribution,
  # and calculate the observed event times, i.e., min(observed event times, censoring times)
  Ctimes <- runif(n, 0, 2 * mean.Cens)
  fixed_visits <- matrix(times, ncol = n)
  fv_idx <- c(1,3,5, seq(9,K, 8))
  if (K %in% fv_idx) {
    fixed_visits <- fixed_visits[fv_idx,]
  } else {
    fixed_visits <- fixed_visits[c(fv_idx, K),]
  }
  fixed_visits_cmpl <- matrix(rbinom(length(fixed_visits), 1, Bcomp),
                         nrow(fixed_visits), ncol(fixed_visits)) # the compliance rate of biopsies
  fixed_visits <- fixed_visits * fixed_visits_cmpl # updated biopsies with compliance rate

  Time.prg <- trueTimes1
  Time.prg_obs <- sapply(1:n, function(i) {
    min(fixed_visits[fixed_visits[,i] >= trueTimes1[i],i])
  })
  Time.trt <- trueTimes2
  Time.cen <- Ctimes
  Time <- pmin(Ctimes, Time.prg_obs, Time.trt)

  event <- sapply(1:n, function(i) {which.min(c(Time.cen[i], Time.prg_obs[i], Time.trt[i]))}) - 1 # event indicator

  Time2 <- Time # the observed event time

  Time[event == 1] <- Time.prg[event == 1] # the true event time

  Time1 <- sapply(1:n, function(i) {
    max(fixed_visits[fixed_visits[,i] <= Time[i],i]) # the last biopsy time
  })

  # sum(Time2 < Time1)
  # sum(Time1 <0 & Time2 <0)

  times.mat <- matrix(times, ncol = n)
  ind <- sapply(1:n, function(i) {times.mat[,i] <= rep(Time2[i], K)})

  # add compliance rate for PSA
  ind.cmpl <- rbinom(length(ind), 1, Pcomp)
  ind <- matrix(as.logical(ind * ind.cmpl), nrow(ind), ncol(ind))

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
  dat$time.prg <- Time.prg[id.obs]
  dat$time.trt <- Time.trt[id.obs]
  dat$time.cen <- Time.cen[id.obs]

  dat <- dat[c("CISNET_ID", "TimeSince_Dx", "PSAValue", "time", "time.cmp1",
               "time.cmp2","status.cmp", "density","DxAge",
               "b1", "b2", "b3", "b4", "time.prg", "time.trt", "time.cen")]

  dat.id <- dat[!duplicated(dat$CISNET_ID), ]

  return(list(dat = dat,
              dat.id = dat.id))
}
