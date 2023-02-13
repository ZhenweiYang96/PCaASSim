#' A simulation function based on MCICJM fitted on PASS data
#'
#' This function allows you to simulate the outcomes (censoring, progression, early treatment) of patients in AS based on PSA value and baseline age and PSA density. \cr
#' In addition to \code{icjmsim()}, this function considers the biopsy sensitivity (0.6 to 0.9 with a step pf 0.05).
#' @import MASS
#' @import Matrix
#' @import JMbayes
#' @import splines
#' @import truncnorm
#' @import mathjaxr
#' @param n number of patients on active surveillance.
#' @param seed seed value for reproducibility The default value is 100.
#' @param param_list the estimated parameters from the PASS data and the fitted MCICJM: \cr
#' \itemize{
#'   \item \code{t.max} - maximum follow-up time in the PASS data
#'   \item \code{mean.Cens} - mean censoring time in the PASS data
#'   \item \code{knot.longi} - knots used in the natural cubic spline specification of the longitudinal model
#'   \item \code{knot.surv} - knots used in the P-splines of the baseline hazard in the survival model
#'   \item \code{age.mean} - mean age observed in the PASS data (centered around 62)
#'   \item \code{age.sd} - standard deviation of age observed in the PASS data
#'   \item \code{density.mean} - mean of the PSA density (in log) observed in the PASS data
#'   \item \code{density.sd} - standard deviation of the PSA density (in log) observed in the PASS data
#'   \item \code{cvisit.sep} - regular clinical visit interval for PSA measurement. The default is 3 months (0.25 years)
#'   \item \code{cvisit.sd} - standard deviation of variation in the clinical visit time
#' }
#' @param sensitivity biopsy sensitivity, can be chose from (0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9).
#' @param Pcomp the compliance rate of PSA measurements. The default value is 1.
#' @param Bcomp the compliance rate of biopsies. The default value is 1.
#' @details The simulated data are based on seven MCICJMs fitted on the Canary PASS data with fixed biopsy sensitivity ranging from 0.6 to 0.9 (with a step of 0.05). According to the user-specified sensitivity in the above-mentioned range, the following parameters are corresponding chosen \cr
#' \loadmathjax
#' \itemize{
#'   \item \code{betas} - coefficients for the fixed effects in the longitudinal submodel for PSA
#'   \item \code{sigma.y} - the standard deviation for the residuals in the longitudinal submodel for PSA
#'   \item \code{D_c3} - the covariance matrix for the random effects in the longitudinal submodel for PSA
#'   \item \code{gambh} - coefficients for the P-spline design matrix of the baseline hazard in the survival submodel
#'   \item \code{gammas} - coefficients for the exogenous covariates in the survival submodel;
#'   \item \code{alpha} - coefficients for the impact from PSA on the time-to-event outcomes in the survival submodel
#' } \cr
#' The `true event time` is simulated based on the event-specific hazard and bounded with the observed maximum follow-up time (12.48 yrs) in the PASS data. The `censoring time` follows a uniform distribution with the upper bound as the doubled mean of the observed censoring time in the PASS data.
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
#' mcicjmsim(n = 1e5, )
#' @export


mcicjmsim <- function(n = 1000, seed = 100,
                      param_list = list(
                        t.max = 12.48,
                        mean.Cens = 5.114849,
                        knot.longi = c(0, 1.36, 2.96, 12.48),
                        knot.surv = c(0, 0, 0, 0, 1.386667, 2.773333, 4.160000, 5.546667, 6.933333, 8.320000, 9.706667, 11.093333, 12.48, 12.48, 12.48, 12.48),
                        age.mean = -0.1248499,
                        age.sd = 6.86713,
                        density.mean = -2.272663,
                        density.sd = 0.6042731,
                        cvisit.sep = 0.25,
                        cvisit.sd = 0.036
                      ), sensitivity = 0.6, Pcomp = 1, Bcomp = 1) {


  if (!(sensitivity %in% seq(0.6, 0.9, 0.05))) {
    stop("The biopsy sensitivity can only be choosen from {0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9}!")
  }

  if (sensitivity == 0.6) {
    betas = c(2.352129324, 0.272079034, 0.649459222, 1.043220801, 0.015353075)
    sigma.y = 0.1452562
    D_c3 = matrix(c(0.487468820 , -0.037135406, -0.082449747,  0.009012073,
                    -0.037135406, 0.777292878,  0.434979814, -0.089847827,
                    -0.082449747,  0.434979814,  1.432122696,  1.475295280,
                    0.009012073, -0.089847827,  1.475295280,  2.709984093), 4, 4)
    gambh = matrix(c(-2.665400754, -2.271979185, -1.894803896, -1.655702846, -1.641588008, -1.803504514, -2.080335640, -2.423645090, -2.788806058, -3.154227010, -3.537907931, -3.935835573,
                     -5.203236363, -4.593533092, -4.272369528, -4.336942011, -4.515441308, -4.668240880, -4.776437392, -4.894436058, -5.038643556, -5.232683167, -5.448271313, -5.674776358), 12, 2)
    gammas = c(0.431700436, 0.245930039)
    alpha = matrix(c(0.151630901, 1.898640758, 0.398278261, 2.269081065), 2, 2)
  } else if (sensitivity == 0.65) {
    betas = c(2.35916156, 0.28071219, 0.65276735, 1.04312064, 0.01553781)
    sigma.y = 0.1451515
    D_c3 = matrix(c(0.48884643 , -0.03559557, -0.08159353,  0.01088569,
                    -0.03559557, 0.78427933,  0.41951278, -0.11499152,
                    -0.08159353,  0.41951278,  1.47540723,  1.55584195,
                    0.01088569, -0.11499152,  1.55584195,  2.86800354), 4, 4)
    gambh = matrix(c(-2.97732953, -2.55641402, -2.13699510, -1.86784715, -1.84017159, -1.97313631, -2.21714671, -2.50487701, -2.78744611, -3.06344438, -3.32069568, -3.60509819,
                     -5.25539587, -4.62908939, -4.31838007, -4.38518099, -4.55201518, -4.73041562, -4.87841681, -4.99229966, -5.12067110, -5.28471702, -5.46297205, -5.65432410), 12, 2)
    gammas = c(0.38845430, 0.23130498)
    alpha = matrix(c(0.18026827, 1.88920829, 0.40362622, 2.30624195), 2, 2)
  } else if (sensitivity == 0.7) {
    betas = c(2.34896747, 0.26474019, 0.61614647, 0.98789633, 0.01823955)
    sigma.y = 0.1453678
    D_c3 = matrix(c(0.48696044 , -0.04170548, -0.07496125,  0.02288856,
                    -0.04170548, 0.77123411,  0.45075706, -0.05867300,
                    -0.07496125,  0.45075706,  1.35408590,  1.33532107,
                    0.02288856, -0.05867300,  1.33532107,  2.45698367), 4, 4)
    gambh = matrix(c(-2.95256480, -2.53111603, -2.13401507, -1.84315853, -1.76027260, -1.85808493, -2.09596824, -2.35593965, -2.62545995, -2.87886429, -3.13594218, -3.38983808,
                     -5.21628282, -4.58699991, -4.29018145, -4.34477514, -4.51164100, -4.67731307, -4.80494661, -4.90464411, -5.04741293, -5.26605583, -5.47555666, -5.68352060), 12, 2)
    gammas = c(0.41174044, 0.24861549)
    alpha = matrix(c(0.16779912, 1.80523286, 0.40977825, 2.20759527), 2, 2)
  } else if (sensitivity == 0.75) {
    betas = c(2.353551638, 0.269237668, 0.623986138, 1.002590218, 0.015519893)
    sigma.y = 0.1452596
    D_c3 = matrix(c(0.487678959 , -0.036874777, -0.089776421,  -0.003862647,
                    -0.036874777, 0.773825259,  0.432616283, -0.084097818,
                    -0.089776421,  0.432616283,  1.413729969,  1.426883734,
                    -0.003862647, -0.084097818,  1.426883734,  2.595005578), 4, 4)
    gambh = matrix(c(-3.022865564, -2.574942908, -2.165002586, -1.868394333, -1.780851418, -1.866034743, -2.043196689, -2.260330582, -2.506676980, -2.740856927, -2.954514094, -3.146719267,
                     -5.128957296, -4.551501961, -4.263715050, -4.308050583, -4.474338093, -4.645207487, -4.802948164, -4.912553004, -5.113348423, -5.372396769, -5.660759223, -5.967330241), 12, 2)
    gammas = c(0.414315298, 0.251801149)
    alpha = matrix(c(0.162897891, 1.788904469, 0.399214083, 2.222056696), 2, 2)
  } else if (sensitivity == 0.8) {
    betas = c(2.356588481, 0.272948830, 0.644815621, 1.033282724, 0.016234042)
    sigma.y = 0.1452074
    D_c3 = matrix(c(0.486820011 , -0.038621202, -0.082794049,  0.008926739,
                    -0.038621202, 0.770276184,  0.440816481, -0.084097818,
                    -0.082794049,  0.440816481,  1.420830851,  1.446250164,
                    0.008926739, -0.084097818,  1.446250164,  2.644159252), 4, 4)
    gambh = matrix(c(-3.165103149, -2.769622909, -2.370050023, -2.069669115, -1.926431838, -1.955697588, -2.122101152, -2.371421927, -2.609462994, -2.851727184, -3.087294192, -3.325937065,
                     -5.360538803, -4.713752544, -4.378904766, -4.434691999, -4.627262507, -4.800936028, -4.961507735, -5.075348794, -5.217484592, -5.426660617, -5.628647195, -5.812534988), 12, 2)
    gammas = c(0.385386823, 0.220831854)
    alpha = matrix(c(0.188245170, 1.748985371, 0.423016178, 2.261760329), 2, 2)
  } else if (sensitivity == 0.85) {
    betas = c(2.352953430, 0.275257794, 0.597686113, 0.944273172, 0.016290199)
    sigma.y = 0.1452433
    D_c3 = matrix(c(0.488642352 , -0.036540460, -0.092497113,  -0.009357853,
                    -0.036540460, 0.776358730,  0.421805357, -0.110353303,
                    -0.092497113,  0.421805357,  1.447349602,  1.509494709,
                    -0.009357853, -0.110353303,  1.509494709,  2.782744172), 4, 4)
    gambh = matrix(c(-3.226348200, -2.759855160, -2.326239353, -1.990560306, -1.865688229, -1.877771556, -2.032542847, -2.273768270, -2.531882698, -2.829122059, -3.114709292, -3.380424364,
                     -5.408805589, -4.792106879, -4.487409068, -4.549508620, -4.729942716, -4.883242980, -5.024760368, -5.148596877, -5.323965018, -5.529580632, -5.754638580, -5.971977563), 12, 2)
    gammas = c(0.404188401, 0.202586480)
    alpha = matrix(c(0.169192069, 1.761858497, 0.449120817, 2.173414918), 2, 2)
  } else if (sensitivity == 0.9) {
    betas = c(2.353453856, 0.277881395, 0.644051169, 1.021577864, 0.014977841)
    sigma.y = 0.1452473
    D_c3 = matrix(c(0.487862642 , -0.038482497, -0.085932777,  -0.009357853,
                    -0.038482497, 0.767263971,  0.442514087, -0.064342839,
                    -0.085932777,  0.442514087,  1.471572079,  1.504941133,
                    -0.009357853, -0.064342839,  1.504941133,  2.699179164), 4, 4)
    gambh = matrix(c(-3.254276234, -2.803285544, -2.357368545, -2.015589272, -1.854142561, -1.858645160, -1.979729984, -2.211776079, -2.492273565, -2.777049276, -3.067553623, -3.349756364,
                     -5.321748676, -4.713794514, -4.401910430, -4.464901247, -4.634442679, -4.787715027, -4.946660705, -5.070132506, -5.250413918, -5.497500996, -5.745547531, -5.955826158), 12, 2)
    gammas = c(0.406758832, 0.217842503)
    alpha = matrix(c(0.161481731, 1.754116382, 0.422726916, 2.288462247), 2, 2)
  }

  # Preparation
  f.org <- JMbayes:::gaussKronrod()$sk
  w <- JMbayes:::gaussKronrod()$wk
  knot.longi <- param_list$knot.longi
  knot.surv <- param_list$knot.surv

  set.seed(seed)

  # Parameters
  t.max <- param_list$t.max
  # betas <- param_list$betas
  # sigma.y <- param_list$sigma.y
  # gammas <- param_list$gammas
  # alpha <- param_list$alpha
  # gambh <- param_list$gambh
  mean.Cens <- param_list$mean.Cens
  Bsens <- sensitivity
  V <-  D_c3 #param_list$
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
