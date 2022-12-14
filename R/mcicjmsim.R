#' A simulation function based on MCICJM fitted on PASS data
#'
#' This function allows you to simulate the outcomes (censoring, progression, early treatment) of patients in AS based on PSA value and baseline age and PSA density. \cr
#' In addition to \code{icjmsim()}, this function considers the biopsy sensitivity of 0.5.
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
#'   \item \code{sensitivity} - biopsy sensitivity, set to be 0.5
#' }
#' @param Pcomp the compliance rate of PSA measurements. The default value is 1.
#' @param Bcomp the compliance rate of biopsies. The default value is 1.
#' @details The `true event time` is simulated based on the event-specific hazard and bounded with the observed maximum follow-up time (12.48 yrs) in the PASS data. The `censoring time` follows a uniform distribution with the upper bound as the doubled mean of the observed censoring time in the PASS data.
#' @return three datasets. First dataset records all the longitudinal measurements (used in the longitudinal submodel), second records each subject per row (used in the survival submodel). Each datasets contains the following columns: \cr
#' \loadmathjax
#' \itemize{
#'  \item \code{CISNET_ID} - patient ID
#'  \item \code{TimeSince_Dx} - PSA time since start of AS
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
#' } \cr
#' The third dataset contains all the biopsies each patient participated in and has the following columns: \cr
#' \loadmathjax
#' \itemize{
#'  \item \code{CISNET_ID} - patient ID
#'  \item \code{TimeSince_Dx} - biopsy time since start of AS
#'  \item \code{compliance} - whether the patient participates in this biopsy
#'  \item \code{bioresult} - whether the patient is detected with progression
#'  \item \code{status.cmp} - observed event indicator (1 = cancer progression; 2 = early treatment; 0 = censoring)
#'  \item \code{time} - true event time in practice
#'  \item \code{time.cmp2} - treatment-free time, which is also the time until which the patient is followed
#' }
#' @keywords MCICJM Simulation
#' @examples
#' mcicjmsim(n = 1e5)
#' @export


mcicjmsim <- function(n = 1000, seed = 100,
                    param_list = list(
                      t.max = 12.48,
                      mean.Cens = 5.114849,
                      knot.longi = c(0, 1.36, 2.96, 12.48),
                      knot.surv = c(0, 0, 0, 0, 1.386667, 2.773333, 4.160000, 5.546667, 6.933333, 8.320000, 9.706667, 11.093333, 12.48, 12.48, 12.48, 12.48),
                      betas = c(2.356578697, 0.274391033, 0.612221625, 0.972303239, 0.016394291),
                      sigma.y = 0.1452828,
                      D_c3 = matrix(c(0.484457235 , -0.037087017, -0.086005920,  -0.003621471,
                                      -0.037087017, 0.779744556,  0.424119912, -0.104367978,
                                      -0.086005920,  0.424119912,  1.423327185,  1.474941751,
                                      -0.003621471, -0.104367978,  1.474941751,  2.745178027), 4, 4),
                      gambh = matrix(c(-4.554247885, -2.431267662, -0.726208821, -0.295394857, -0.685145348, -1.199522121, -1.701821794, -2.134330688, -2.382277522, -2.512485057, -2.549171766, -2.552028383,
                                       -5.259334899, -4.639982156, -4.332605036, -4.385041310, -4.562366014, -4.731760452, -4.892164516, -4.995863325, -5.156619515, -5.356645767, -5.542735349, -5.719399707), 12, 2),
                      gammas = c(0.667693823, 0.239728042),
                      alpha = matrix(c(0.078023016, 2.916081305, 0.422728218, 2.131612465), 2, 2),
                      age.mean = -0.1248499,
                      age.sd = 6.86713,
                      density.mean = -2.272663,
                      density.sd = 0.6042731,
                      cvisit.sep = 0.25,
                      cvisit.sd = 0.036,
                      sensitivity = 0.5
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
  Bsens <- param_list$sensitivity
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

      exp(bh %*% gambh[,k]  + eta.t[i,k] + f * alpha[1,k] + (f - f.back) * alpha[2,k]) # alpha, row = assoc; col = nrisk

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
  fv_idx <- c(1, 5, seq(9,K, 8))
  if (K %in% fv_idx) {
    fixed_visits <- fixed_visits[fv_idx,]
  } else {
    fixed_visits <- fixed_visits[c(fv_idx, K),]
  }
  fixed_visits_cmpl <- matrix(rbinom(length(fixed_visits), 1, Bcomp),
                              nrow(fixed_visits), ncol(fixed_visits)) # the compliance rate of biopsies
  fixed_visits_sens <- matrix(rbinom(length(fixed_visits), 1, Bsens),
                              nrow(fixed_visits), ncol(fixed_visits)) # the sensitivity of biopsies
  fixed_visits_obstime <- fixed_visits * fixed_visits_cmpl * fixed_visits_sens # updated biopsies with compliance rate and sensitivity
  fixed_visits_time1 <- fixed_visits * fixed_visits_cmpl

  Time.prg <- trueTimes1
  Time.prg_obs <- sapply(1:n, function(i) {
    ifelse(max(c(trueTimes1[i], fixed_visits_obstime[,i])) == trueTimes1[i],
           t.max + 1,
           min(fixed_visits_obstime[fixed_visits_obstime[,i] >= trueTimes1[i],i]))
  })
  Time.trt <- trueTimes2
  Time.cen <- Ctimes
  Time <- pmin(Ctimes, Time.prg_obs, Time.trt)

  event <- sapply(1:n, function(i) {which.min(c(Time.cen[i], Time.prg_obs[i], Time.trt[i]))}) - 1 # event indicator

  Time2 <- Time # the observed event time

  Time1 <- sapply(1:n, function(i) {
    max(fixed_visits_time1[fixed_visits_time1[,i] < Time[i],i]) # the last biopsy time
  })

  Time[event == 1] <- Time.prg[event == 1] # the true event time

  # sum(Time2 < Time1)
  # sum(Time1 <0 & Time2 <0)

  fixed_visits_senslog <- matrix(0, nrow(fixed_visits), ncol(fixed_visits))
  for (i in 1:n) {
    if (event[i] == 1) {
      fixed_visits_senslog[fixed_visits[,i] == Time.prg_obs[i],i] <- 1
    }
  }


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

  biopsy_full <- data.frame(CISNET_ID = rep(1:n, each = nrow(fixed_visits)),
                            TimeSince_Dx = c(fixed_visits),
                            compliance = as.logical(fixed_visits_time1),
                            bioresult = as.logical(fixed_visits_senslog),
                            status.cmp = rep(event, each = nrow(fixed_visits)),
                            time = rep(Time, each = nrow(fixed_visits)),
                            time.cmp2 = rep(Time2, each = nrow(fixed_visits)))

  biopsy <- biopsy_full[biopsy_full$compliance & biopsy_full$TimeSince_Dx <= biopsy_full$time.cmp2, ]

  return(list(dat = dat,
              dat.id = dat.id,
              biopsy = biopsy))
}
