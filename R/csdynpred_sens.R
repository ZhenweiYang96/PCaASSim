#' A dynamic prediction for cancer progression based on MCICJM also incorporating sensitivity
#'
#' This function allows you to predict patient-specific risk of cancer progression based on the fitted MCICJM given a certain sensitivity value (or a sensitivity prior). \cr
#' The risk refers to \loadmathjax \mjdeqn{\Pr(T^*_i = T^{\small{prg},*}_i < t^{(p)} \mid T^{\small{prg}+}_i > t^{(b)})}{ASCII representation}
#' The prediction can be dynamically updated with new PSA values.
#'
#' @import splines
#' @import mvtnorm
#' @import MASS
#' @import JMbayes
#' @import rjags
#' @import mathjaxr
#' @param sens_fit the sensitivity parameter assumed in the MCICJM fitted on the training (PASS) data. The default value is 0.6.
#' @param sens_pred the fixed sensitivity value used in the prediction. If not specified, will automatically extract the value from posterior distribution of the sens parameter from the fitted model.
#' @param model the self-fitted MCICJM model. If not specified, the default fitted MCICJMs (according to sensitivity) is used.
#' @param t_biopsies the times of prior biopsies (excluding the one at year 0, the start of active surveillance).
#' @param t_start the start time point (in years) to predict, e.g., usually the last biopsy. The default value is 0.
#' @param t_visits the time point(s) of clinical visit, whose risks are of interest. \
#' @param n.step the interval of each point estimates. The default value is 50.
#' @param iter iteration number of the MCMC sampling after adaption. The default value is 250.
#' @param n.adapt iteration number of MCMC sampling in the adaption phase. The default value is 150.
#' @param seed seed value for reproducibility The default value is 100.
#' @param idVar the name of the ID variable in the testset. The default is "CISNET_ID".
#' @param varname the names of related variables, default is the same in the PASS data: \cr
#' \itemize{
#'   \item \code{longitimeVar} - the variable name of PSA measurement times
#'   \item \code{longi.bsVar} - the variable name of baseline covariate (i.e., age centered around 62) in the longitudinal submodel
#'   \item \code{surv.bsVar} - the variable name of baseline covariate (i.e., PSA density in log scale) in the survival submodel
#' }
#' @param newdata the newdata contains the PSA measurements and baseline covariates of only one patient.
#' @param sigma.inits the initial value of \mjeqn{\sigma}{ASCII representation} in Robbins-Monro process. The default value is 0.27.
#' @param sigMat.inits the initial value of \mjeqn{\Sigma}{ASCII representation} in Robbins-Monro process.
#' @param qp the number of Gaussian quadrature point in the numerical approximation for the integral. The default value is 15.
#' @details The prediction of cancer progression are based on seven MCICJMs fitted on the Canary PASS data with fixed biopsy sensitivity ranging from 0.6 to 0.9 (with a step of 0.05). According to the user-specified sensitivity in the above-mentioned range, the corresponding MCICJM is chosen as the training model. \cr
#' This function can only generate predictions based on a dataset containing only ONE patient. This facilitates the users to make prediction in parallel (e.g., using packages \code{doParallel} or \code{future}) if they have larger datasets.
#' @return a list of following components: \cr
#' \loadmathjax
#' \itemize{
#'  \item \code{risk_visits} - the table of predicted risks at the clinical visits of interest, containing four columns: \cr
#'  \loadmathjax
#'  \itemize{
#'    \item \code{t.hor} - time point
#'    \item \code{mean} - mean estimates
#'    \item \code{lower} - lower bound of the 95\% credible interval
#'    \item \code{upper} - lower bound of the 95\% credible interval.
#'  } \cr
#'  \item \code{risk} - the risk estimates (with uncertainty) of all time points. This can be used to draw a risk curve.
#'  \item \code{priorbiopsy_info} - a list includes two tables:
#'  \loadmathjax
#'  \itemize{
#'    \item \code{priorbsp.risk} - risk estimates for each prior biopsy intervals (rows = iterations, cols = biopsy intervals)
#'    \item \code{baseline.bp.risk} - the offsets of the risks attributed from the biopsy sensitivity in all iterations
#'  } \cr
#'  \item \code{longi} - the estimated PSA trajectory (with uncertainty).
#'  \item \code{acc} - the acceptance rate in the Robbins-Monro process. This is used to check convergence.
#'  \item \code{surv} - the vector of survival probability until \code{t_start}.
#'  \item \code{parameters} - the MCMC samples of the parameters from the MCICJM model: \cr
#'  \loadmathjax
#'  \itemize{
#'    \item \code{beta} - fixed effects in the longitudinal submodel
#'    \item \code{b} - patient-specific random effects in the longitudinal submodel
#'    \item \code{gam_bh} - the coefficients in the baseline hazard
#'    \item \code{gamma} - the coefficient of the baseline covariate (PSA density) in the survival submodel
#'    \item \code{alpha} - the association parameters of PSA in different functional forms
#'    \item \code{Sigma} - the variance covariance matrix of the random effects in the longitudinal submodel
#'    \item \code{tau} - the inverse of the residuals variance in the longitudinal submodel
#'    \item \code{sens} - the sensitivity parameter used in each iteration, either fixed or sampled from posterior
#'  } \cr
#'  \item \code{inits} - the initial values for \mjeqn{\sigma}{ASCII representation} and \mjeqn{\Sigma}{ASCII representation} matrix in the Robbins-Monro process.
#'  \item \code{overall.surv} - a table of the estimates overall survival probabilities at all time points.
#'  \item \code{full.results} - a list of calculations from MCMC samples at each time points.
#' }
#' @keywords Dynamic prediction
#' @examples
#' csdypred_sens(sens_pred = 0.8, model = mcicjm, t_biopsy = 1:2, t_start = 2, t_visits = c(1, 2, 4, 6), newdata = test)
#' @export

csdypred_sens <- function(sens_fit = 0.6,
                          sens_pred = NULL,
                          model = NULL,
                          t_biopsies = NULL,
                          t_start = 0,
                          t_visits = NULL,
                          n.step = 50L,
                          iter = 250L, n.adapt = 150L,
                          seed = 100L,
                          idVar = "CISNET_ID",
                          varname = list(longitimeVar = NULL,
                                         #survtimeVar = NULL,
                                         longi.bsVar = NULL,
                                         surv.bsVar = NULL),
                          newdata,
                          sigma.inits = 0.27,
                          sigMat.inits = NULL,
                          qp = 15) {
  if (is.null(newdata)) {
    stop("Test set is missing! \n")
  }

  if (is.null(t_visits)) {
    stop("The clinical visits time points are missing\n")
  }

  if (is.null(t_biopsies) & t_start!=0) {
    stop("The time points of the past biopsies are missing. You have to provide the prior biopsy times!\n")
  }

  if (!is.null(t_biopsies) & max(t_biopsies) != t_start) {
    stop("The last conducted biopsy time is not consistent in 't_start' and 't_biopsies")
  }

  if (length(unique(newdata[[idVar]])) > 1) {
    stop("Dynamic predictions can be generated for only one patient at a time.\n")
  }

  if (is.null(model)) {
    object <- get(paste0("mcicjm.model.", sens_fit*100))
  } else {
    object <- model
  }

  sens_included <- any(names(object$summary$coef) == "sens")
  if (is.null(sens_pred) & !sens_included) {
    stop("The MCICJM does not have a sensitivity prior. Thus, the exact sensitvity value must be specified.")
  }

  # Survtimes <- NULL
  # if (t_start == 0) {
  #   Survtimes[[1]] <- sort(c(t_visits,
  #                            seq(t_start, max(t_visits), length.out = n.step)))
  #   Survtimes[[1]] <- Survtimes[[1]][!duplicated(Survtimes[[1]])] # remove duplicated time points
  # } else {
  #   t_biopsies <- c(0, t_biopsies)
  #   for (i in 1:length(t_biopsies)) {
  #     Survtimes[[i]] <- seq(t_biopsies[i], t_biopsies[i+1],
  #                           length.out = n.step.past)
  #   }
  #   Survtimes[[length(t_biopsies)+1]] <- sort(c(t_visits,
  #                            seq(t_start, max(t_visits), length.out = n.step)))
  #   Survtimes[[length(t_biopsies)+1]] <- Survtimes[[length(t_biopsies)+1]][!duplicated(Survtimes[[length(t_biopsies)+1]])] # remove duplicated time points
  # }

  Survtimes <- sort(c(t_visits,
                      seq(t_start, max(t_visits), length.out = n.step)))
  Survtimes <- Survtimes[!duplicated(Survtimes)]

  u <- max(Survtimes)
  t <- min(Survtimes)

  # time points of the prior biopsies
  if (t_start != 0) {
    t.priorbsp <- c(0, t_biopsies)
  }

  # extract the posterior distribution

  jags_res <- extract_pd(object$mcmc)
  pool.bs_gammas <- jags_res$pool.bs_gammas
  pool.gammas <- jags_res$pool.gammas
  pool.alphas <- jags_res$pool.alphas
  pool.beta <- jags_res$pool.beta
  pool.D <- jags_res$pool.D
  pool.tau <- jags_res$pool.tau
  if (sens_included) {
    pool.sens <- jags_res$pool.sens
  }
  # extract model info from JAGS
  # idVar <- ifelse(is.null(varname$idVar),
  #                 object$model_info$var_names$id,
  #                 varname$idVar)

  test <- newdata
  test.id <- test[!duplicated(test[[idVar]]),]
  df.ns.longi <- sqrt(ncol(pool.D)) - 1
  knot.longi <- object$model_info$knots$knot.longi
  knot.surv <- object$model_info$knots$knot.surv


  # survtimeVar <- ifelse(is.null(varname$survtimeVar),
  #                       object$model_info$var_names$survtime,
  #                       varname$survtimeVar)
  longitimeVar <- ifelse(is.null(varname$longitimeVar),
                         object$model_info$var_names$longitime,
                         varname$longitimeVar)
  longi.bsVar <- ifelse(is.null(varname$longi.bsVar),
                        object$model_info$var_names$longi.bsVar,
                        varname$longi.bsVar)
  surv.bsVar <- ifelse(is.null(varname$surv.bsVar),
                       object$model_info$var_names$surv.bsVar,
                       varname$surv.bsVar)

  # preparation for sample from posterior distribution
  step <- Survtimes
  total <- iter

  beta <- matrix(NA, total, ncol(pool.beta))
  b <- matrix(NA, total, df.ns.longi + 1)
  tau <- matrix(NA, total, 1)
  gam_bh <- array(NA, dim = c(total, ncol(pool.bs_gammas)/2, 2))
  gamma <- array(NA, dim =c(total, ncol(pool.gammas)/2, 2))
  Sigma <- array(NA, dim = c(total, df.ns.longi + 1, df.ns.longi + 1))
  alpha <- array(NA, dim = c(total, ncol(pool.alphas)/2, 2)) # col = risk
  sens <- matrix(NA, total, 1)

  # risk output data frame
  risk <- data.frame(t.hor = step,
                     mean = rep(NA, length(step)),
                     upper = rep(NA, length(step)),
                     lower = rep(NA, length(step)))
  # fitted longitudinal model
  longi <- data.frame(t.hor = seq(0, max(test[[longitimeVar]]), length.out = 50),
                      PSA.pred = rep(NA,50))

  f.org <- JMbayes:::gaussKronrod(qp)$sk
  w <- JMbayes:::gaussKronrod(qp)$wk

  # risk store set
  risk.store <- matrix(NA, iter, length(step))
  numerator.risk.store <- matrix(NA, iter, length(step))
  denomiator.risks.store <- matrix(NA, iter, length(step))

  # sample
  set.seed(seed)
  idx <- sample(1:nrow(pool.beta), total, replace = T)
  beta <-  pool.beta[idx,]

  tau[,1] <- pool.tau[idx]
  if (sens_included) {sens[,1] <- pool.sens[idx]}
  else {sens[,1] <- rep(sens_pred, total)}

  for (i in 1:ncol(pool.bs_gammas)/2) {
    for (j in 1:2) {
      gam_bh[,i,j] <- pool.bs_gammas[idx, i + (j-1) * ncol(pool.bs_gammas)/2]
    }
  }

  for (i in 1:ncol(pool.gammas)/2) {
    for (j in 1:2) {
      gamma[,i,j] <- pool.gammas[idx, i + (j-1) * ncol(pool.gammas)/2]
    }
  }

  for (i in 1:length(idx)) {
    Sigma[i,,] <- pool.D[idx[i],]
  }

  for (i in 1:length(idx)) {
    alpha[i,,] <- pool.alphas[idx[i],]
  }

  #### sample random effects
  n_test <- nrow(test)
  XL <- cbind(1,
              ns(test[[longitimeVar]], knots = tail(head(knot.longi,-1),-1), B = knot.longi[c(1,df.ns.longi+1)]),
              test[[longi.bsVar]])
  ZL <- cbind(1,
              ns(test[[longitimeVar]], knots = tail(head(knot.longi,-1),-1), B = knot.longi[c(1,df.ns.longi+1)]))
  XL.b.now <- matrix(1, qp, ncol(pool.beta))
  ZL.b.now <- matrix(1, qp, df.ns.longi+1)
  XL.b.now[,(df.ns.longi+2):ncol(pool.beta)] <- test.id[[longi.bsVar]]
  XL.b.now[,2:(df.ns.longi+1)] <- ZL.b.now[,2:(df.ns.longi+1)] <- ns(f.org * t/2 + t/2,
                                                                     knots = tail(head(knot.longi,-1),-1), B = knot.longi[c(1,df.ns.longi+1)])

  XL.b.back <- matrix(1, qp, ncol(pool.beta))
  ZL.b.back <- matrix(1, qp, df.ns.longi+1)
  XL.b.back[,(df.ns.longi+2):ncol(pool.beta)] <- test.id[[longi.bsVar]]
  XL.b.back[,2:(df.ns.longi+1)] <- ZL.b.back[,2:(df.ns.longi+1)] <- ns(f.org * t/2 + t/2 - 1,
                                                                       knots = tail(head(knot.longi,-1),-1), B = knot.longi[c(1,df.ns.longi+1)])
  #outer(ns(test.id$time, knots = knot.longi[2:3], B = knot.longi[c(1,4)])/2, f.org + 1)

  X <- test.id[[surv.bsVar]]
  T.b <- splineDesign(knot.surv, outer(t/2, f.org + 1), ord = 4L)


  loglike <- function(b) {
    n_test/2 * log(tau[i,] / (3*pi))  - 2 * sum(sapply(1:n_test, function(n) {
      log(1 + tau[i,]/3 * (test$PSAValue[n] - (XL[n,,drop = F] %*% beta[i,] + ZL[n,, drop = F] %*% b))^2)
    }))  -
      # n_test/2 * log(tau[i,]) - 0.5 * tau[i,] * sum((test$PSAValue - (XL %*% beta[i,] + ZL %*% b))^2) -
      t/2 * sum(sapply(1:2, function(k) {
        exp(gamma[i,,k] * X) * w %*% sapply(1:qp, function(l) {
          exp(T.b[l,] %*% gam_bh[i,,k] + alpha[i,1,k] * (XL.b.now[l,] %*% beta[i,] + ZL.b.now[l,] %*% b) +
                alpha[i,2,k] * (XL.b.now[l,] %*% beta[i,] + ZL.b.now[l,] %*% b -
                                  XL.b.back[l,] %*% beta[i,] - ZL.b.back[l,]%*% b))
        })
      })) -
      0.5 * log(det(Sigma[i,,])) - 0.5 * matrix(b,1,4) %*% solve(Sigma[i,,]) %*% matrix(b,4,1)
  }

  # RWMH method -------------------------------------------------------------
  update.sigma<- function(sigma2, acc, p=p, j) {
    c=((1-1/d)*sqrt(2*pi)*exp(a^2/2)/(2*a) + 1/(d*p*(1-p)))
    Theta=log(sqrt(sigma2))
    Theta=Theta+c*(acc-p)/max(200, i/d)
    return(exp(Theta))
  }


  update.cov<-function(sigMat, j, thetaM, theta){
    #function to recursively update covariance matrix, as part of adaptive MCMC sampling, updating the covariance matrix
    epsilon=1/j
    thetaM2=((thetaM*j)+theta)/(j+1)
    sigMat=(j-1)/j*sigMat + thetaM%*%t(thetaM)-(j+1)/j*thetaM2%*%t(thetaM2)+1/j*theta%*%t(theta) + epsilon*diag(d)
    return(list(sigMat=sigMat, thetaM=thetaM2))
  }


  #############################################
  # Begin
  #initialise the RM section
  #############################################
  #the optimal acceptance probability for the multivariate case
  pstar=0.234
  a=-qnorm(pstar/2)
  n0=round(5/(pstar*(1-pstar)))
  #iMax, is the max number of iterations before the last restart
  iMax=20
  Numbig=0
  Numsmall=0

  #initialise the MCMC program
  #niter=iter + n.adapt      # niter is the number of iterations in the Markov chain.
  # Change it to the value you want.
  d=df.ns.longi+1 #dimension of parameters to be updated
  output<-rep(0, d) #output parameter values
  output.mat<-output #records all MCMC output for the parameters
  # initial values
  sigma <- sigma.inits
  if (is.null(sigMat.inits)) {sigMat <- diag(d)} else {sigMat <- sigMat.inits}
  meanacc.record <- list()
  # end ---------------------------------------------------------------------

  #### End adaption and start MH
  n.adapt <- c(n.adapt, rep(30, iter - 1))
  for (i in 1:iter) {

    #############################################
    # Begin
    #initialise the RM section
    #############################################
    sigma= sigma       #an arbitrary starting value for sigma, equivalent to theta=ln(sigma)=0
    sigma2=sigma^2 #variance
    sigma.start<- sigma
    sigma.vec<- sigma
    sigMat <- sigMat
    acc.vec<-rep(NA, n.adapt[i]) #records the acceptance probabilities from each MCMC iteration
    num.restart=0
    j=1
    #############################################
    # End
    #initialise the RM section
    #############################################

    for (m in 2:n.adapt[i]) {

      #propose a new value of theta
      output.prop<- as.vector(rmvt(1, sigma=sigma2*sigMat, df = 3) + output)
      pi.old.log <- loglike(output)
      pi.new.log <- loglike(output.prop)
      u<-runif(1)
      acc=min(1, exp(pi.new.log - pi.old.log))
      acc.vec=c(acc.vec,acc)
      j=j+1

      if (u < acc) {
        output<-output.prop
      }
      output.mat<-rbind(output.mat, output)


      #################################################
      # Begin
      # update covariance matrix with adaptive MCMC
      #################################################
      if (j > 100) {
        if (j==101) {
          sigMat=cov(output.mat)
          thetaM=apply(output.mat, 2, mean)
        } else
        {
          tmp=update.cov(sigMat, i, thetaM, output)
          sigMat=tmp$sigMat
          thetaM=tmp$thetaM
        }
      }

      #################################################
      # End
      # update covariance matrix with adaptive MCMC
      #################################################



      ###############################################
      # Begin
      # update sigma using RM
      ###############################################
      if (m>n0) {
        sigma<-update.sigma(sigma2, acc, pstar, i)
        sigma2=sigma^2
        j=j+1
        sigma.vec<-c(sigma.vec, sigma)
        if ((j <= (iMax+n0)) && (Numbig<5 || Numsmall<5)) {
          Toobig<- (sigma > (3*sigma.start))
          Toosmall<-(sigma < (sigma.start/3))

          if (Toobig || Toosmall) {
            #restart the algorithm
            cat("restart the program at", i, "th iteration", "\n")
            sigma.restart<-c(sigma.restart, sigma)
            Numbig<- Numbig + Toobig
            Numsmall <- Numsmall + Toosmall
            i<-n0
            sigma.start<-sigma
          }
        } #end iMax
      }
      ###############################################
      # Begin
      # update sigma using RM
      ###############################################


    } #end niter
    #end of MCMC

    acc.vec <- acc.vec[ceiling(length(acc.vec)/2):length(acc.vec)]
    meanacc<-rep(NA, length(acc.vec))
    for (m in c(1:length(acc.vec))) {
      meanacc[m]=     mean(acc.vec[round(m/2) :m], na.rm = T)
    }
    # if(i == 1) {
    #   meanacc.record = meanacc
    # }
    meanacc.record[[i]] <- meanacc
    # plot(meanacc, type = "l", ylim = c(0, 0.7)); abline(h = 0.234)

    # for(m in 1:n.burnin) {
    #   output.prop<-as.vector(rmvt(1, sigma=sigma2*sigMat, df = 3) + output)
    #   pi.old.log <- loglike(output)
    #   pi.new.log <- loglike(output.prop)
    #   u<-runif(1)
    #   acc=min(1, exp(pi.new.log - pi.old.log))
    #   acc.vec=c(acc.vec,acc)
    #   if (u < acc) {
    #     output<-output.prop
    #   }
    # }
    output.prop<-as.vector(rmvt(1, sigma=sigma2*sigMat, df = 3) + output)
    pi.old.log <- loglike(output)
    pi.new.log <- loglike(output.prop)
    u<-runif(1)
    acc=min(1, exp(pi.new.log - pi.old.log))
    acc.vec=c(acc.vec,acc)
    #j=j+1

    if (u < acc) {
      output<-output.prop
    }
    output.mat<-rbind(output.mat, output)
  }

  # end ---------------------------------------------------------------------

  b <- output.mat
  #
  # beta <- beta[(n.adapt+1):total,]
  # b <- b[(n.adapt+1):total,]
  # tau <- tau[(n.adapt+1):total,]
  # alpha <- alpha[(n.adapt+1):total,,,drop=F]
  # gamma <- gamma[(n.adapt+1):total,,,drop=F]
  # Sigma <- Sigma[(n.adapt+1):total,,]
  # gam_bh <- gam_bh[(n.adapt+1):total,,]

  # start calculating the risk ----------------------------------------------


  #### First start with
  surv <- array(NA, iter)
  full.results <- list()
  overall.surv <- matrix(NA, iter, length(step))

  # denominator - survival part
  t.start <- step[1]
  XL.hs.now <- matrix(1,qp, ncol(pool.beta))
  ZL.hs.now <- matrix(1,qp, df.ns.longi + 1)
  XL.hs.now[,2:(df.ns.longi+1)] <- ZL.hs.now[,2:(df.ns.longi+1)] <- t(sapply(1:qp, function(n) {
    ns(f.org[n] * t.start/2 + t.start/2,
       knots = tail(head(knot.longi,-1),-1), B = knot.longi[c(1,df.ns.longi+1)])
  }))
  XL.hs.now[,(df.ns.longi+2):ncol(pool.beta)] <- test.id[[longi.bsVar]]

  XL.hs.back <- matrix(1, qp, ncol(pool.beta))
  ZL.hs.back <- matrix(1, qp, df.ns.longi + 1)
  XL.hs.back[,2:(df.ns.longi+1)] <- ZL.hs.back[,2:(df.ns.longi+1)] <- t(sapply(1:qp, function(n) {
    ns(f.org[n] * t.start/2 + t.start/2 - 1,
       knots = tail(head(knot.longi,-1),-1), B = knot.longi[c(1,df.ns.longi+1)])
  }))
  XL.hs.back[,(df.ns.longi+2):ncol(pool.beta)] <- test.id[[longi.bsVar]]

  T.hs <- splineDesign(knot.surv, (f.org + 1) * t/2, ord = 4L)

  hazard.survival <-array(NA, dim = c(iter, qp, 2))

  for (j in 1:iter) {
    hazard.survival[j,1:qp,1:2] <- exp(T.hs[1:qp,] %*% gam_bh[j,,1:2] + rep(gamma[j,,1:2] * X, each = qp) +
                                         (XL.hs.now[1:qp,] %*% beta[j,] + ZL.hs.now[1:qp,] %*% b[j,]) %*% alpha[j,1,1:2]  +
                                         (XL.hs.now[1:qp,] %*% beta[j,] + ZL.hs.now[1:qp,] %*% b[j,] -
                                            XL.hs.back[1:qp,] %*% beta[j,] - ZL.hs.back[1:qp,] %*% b[j,]) %*% alpha[j,2,1:2])
    surv[j] <-  exp(- t.start/2 * (rep(1,2) %*% (t(hazard.survival[j,1:qp,1:2]) %*% w)))
  }
  surv[surv < 1e-16] <- 1e-16

  #################################################
  # baseline risk in both numerator and denominator
  #################################################
  baseline.bp.risk <- rep(0, iter)
  priorbsp.risk <- NA
  # if there is no prior biopsies, the baseline risk in the denominator and numerator is 0
  if (t.start != 0) {
    # sens matrix (row = iter, col = sens weights for each interval)
    sens_mat <- sapply(1:length(t_biopsies), function(x) {
      (1 - sens[,1])^(length(t_biopsies) - x + 1)
    })

    priorbsp.risk <- matrix(NA, iter, length(t_biopsies))
    for (i in 1:length(t_biopsies)) {
      tb.low <- t.priorbsp[i]
      tb.upp <- t.priorbsp[i+1]
      tb.step <- f.org * (tb.upp-tb.low)/2 + (tb.upp+tb.low)/2
      T.pb.h <- splineDesign(knot.surv, x = f.org * (tb.upp-tb.low)/2 + (tb.upp+tb.low)/2, ord = 4L)

      XL.pb.h.now <- matrix(1, qp, ncol(pool.beta))
      ZL.pb.h.now <- matrix(1, qp, df.ns.longi+1)
      XL.pb.h.now[,2:(df.ns.longi+1)] <- ZL.pb.h.now[,2:(df.ns.longi+1)] <- ns(f.org * (tb.upp-tb.low)/2 + (tb.upp+tb.low)/2,
                                                                         knots = tail(head(knot.longi,-1),-1), B = knot.longi[c(1,df.ns.longi+1)])
      XL.pb.h.now[,(df.ns.longi+2):ncol(pool.beta)] <- test.id[[longi.bsVar]]

      XL.pb.h.back <- matrix(1, qp, ncol(pool.beta))
      ZL.pb.h.back <- matrix(1, qp, df.ns.longi+1)
      XL.pb.h.back[,2:(df.ns.longi+1)] <- ZL.pb.h.back[,2:(df.ns.longi+1)] <- ns(f.org * (tb.upp-tb.low)/2 + (tb.upp+tb.low)/2 - 1,
                                                                           knots = tail(head(knot.longi,-1),-1), B = knot.longi[c(1,df.ns.longi+1)])
      XL.pb.h.back[,(df.ns.longi+2):ncol(pool.beta)] <- test.id[[longi.bsVar]]

      T.pb.s <- array(1, dim = c(qp, ncol(pool.bs_gammas)/2, qp)) # second qp is the second integral
      XL.pb.s.now <- array(1, dim = c(qp, ncol(pool.beta), qp))
      ZL.pb.s.now <- array(1, dim = c(qp, df.ns.longi+1, qp))
      XL.pb.s.now[,(df.ns.longi+2):ncol(pool.beta),] <- test.id[[longi.bsVar]]

      for (j in 1:qp) {
        XL.pb.s.now[j,2:(df.ns.longi+1),] <- ZL.pb.s.now[j,2:(df.ns.longi+1),] <- t(ns(f.org * tb.step[j]/2 + tb.step[j]/2,
                                                                                 knots = tail(head(knot.longi,-1),-1), B = knot.longi[c(1,df.ns.longi+1)]))
      }

      XL.pb.s.back <- array(1, dim = c(qp, ncol(pool.beta), qp))
      ZL.pb.s.back <- array(1, dim = c(qp, df.ns.longi+1, qp))
      XL.pb.s.back[,(df.ns.longi+2):ncol(pool.beta),] <- test.id[[longi.bsVar]]

      for (j in 1:qp) {
        XL.pb.s.back[j,2:(df.ns.longi+1),] <- ZL.pb.s.back[j,2:(df.ns.longi+1),] <- t(ns(f.org * tb.step[j]/2 + tb.step[j]/2 - 1,
                                                                                   knots = tail(head(knot.longi,-1),-1), B = knot.longi[c(1,df.ns.longi+1)]))
      }

      T.pb.s.prep <- splineDesign(knot.surv, x = outer(f.org + 1, tb.step/2), ord = 4L) # 1:qp - step[1]:f.org[1:qp]
      gk <- split(1:225, 1:qp)  # f.org group
      for (j in 1:qp) {
        T.pb.s[,,j] <- T.pb.s.prep[gk[[j]],]
      }

      priorbsp.survival <- matrix(NA, iter, qp)
      priorbsp.hazard <- matrix(NA, iter, qp)

      for (j in 1:iter) {
        for (m in 1:qp) {
          priorbsp.survival[j,m] <- exp(
            - tb.step[m]/2 * sum(w %*% exp(t(T.pb.s[m,,1:qp]) %*% gam_bh[j,,1:2]  + rep(gamma[j,,1:2] * X, each = qp) +
                                            (t(XL.pb.s.now[m,,1:qp]) %*% beta[j,] + t(ZL.pb.s.now[m,,1:qp]) %*% b[j,]) %*% alpha[j,1,1:2] +
                                            (t(XL.pb.s.now[m,,1:qp]) %*% beta[j,] + t(ZL.pb.s.now[m,,1:qp]) %*% b[j,] -
                                               t(XL.pb.s.back[m,,1:qp]) %*% beta[j,] - t(ZL.pb.s.back[m,,1:qp]) %*% b[j,]) %*% alpha[j,2,1:2]))
          )

        }
        priorbsp.hazard[j,1:qp] <- exp(T.pb.h[1:qp,] %*% gam_bh[j,,1] +  gamma[j,,1] * X +
                                alpha[j,1,1] * (XL.pb.h.now[1:qp,] %*% beta[j,] + ZL.pb.h.now[1:qp,] %*% b[j,]) +
                                alpha[j,2,1] * (XL.pb.h.now[1:qp,] %*% beta[j,] + ZL.pb.h.now[1:qp,] %*% b[j,] -
                                                  XL.pb.h.back[1:qp,] %*% beta[j,] - ZL.pb.h.back[1:qp,] %*% b[j,]))
      }
      priorbsp.risk[1:iter,i] <- (tb.upp-tb.low)/2 * ((priorbsp.hazard[1:iter,] * priorbsp.survival[1:iter,]) %*% w)
    }

    baseline.bp.risk <- sapply(1:iter, function(x) {priorbsp.risk[x,] %*% sens_mat[x,]})
  }
  #######################################
  ### END
  #######################################

  # numerator - cumulative hazard part
  for (i in 1:length(step)) {

    u <- step[i]
    t <- ifelse(i == 1, step[i], step[i-1])
    t.step <- f.org * (u-t)/2 + (u+t)/2
    T.h <- splineDesign(knot.surv, x = f.org * (u-t)/2 + (u+t)/2, ord = 4L)

    XL.h.now <- matrix(1, qp, ncol(pool.beta))
    ZL.h.now <- matrix(1, qp, df.ns.longi+1)
    XL.h.now[,2:(df.ns.longi+1)] <- ZL.h.now[,2:(df.ns.longi+1)] <- ns(f.org * (u-t)/2 + (u+t)/2,
                                                                       knots = tail(head(knot.longi,-1),-1), B = knot.longi[c(1,df.ns.longi+1)])
    XL.h.now[,(df.ns.longi+2):ncol(pool.beta)] <- test.id[[longi.bsVar]]

    XL.h.back <- matrix(1, qp, ncol(pool.beta))
    ZL.h.back <- matrix(1, qp, df.ns.longi+1)
    XL.h.back[,2:(df.ns.longi+1)] <- ZL.h.back[,2:(df.ns.longi+1)] <- ns(f.org * (u-t)/2 + (u+t)/2 - 1,
                                                                         knots = tail(head(knot.longi,-1),-1), B = knot.longi[c(1,df.ns.longi+1)])
    XL.h.back[,(df.ns.longi+2):ncol(pool.beta)] <- test.id[[longi.bsVar]]

    T.s <- array(1, dim = c(qp, ncol(pool.bs_gammas)/2, qp)) # second qp is the second integral
    XL.s.now <- array(1, dim = c(qp, ncol(pool.beta), qp))
    ZL.s.now <- array(1, dim = c(qp, df.ns.longi+1, qp))
    XL.s.now[,(df.ns.longi+2):ncol(pool.beta),] <- test.id[[longi.bsVar]]

    for (j in 1:qp) {
      XL.s.now[j,2:(df.ns.longi+1),] <- ZL.s.now[j,2:(df.ns.longi+1),] <- t(ns(f.org * t.step[j]/2 + t.step[j]/2,
                                                                               knots = tail(head(knot.longi,-1),-1), B = knot.longi[c(1,df.ns.longi+1)]))
    }

    XL.s.back <- array(1, dim = c(qp, ncol(pool.beta), qp))
    ZL.s.back <- array(1, dim = c(qp, df.ns.longi+1, qp))
    XL.s.back[,(df.ns.longi+2):ncol(pool.beta),] <- test.id[[longi.bsVar]]

    for (j in 1:qp) {
      XL.s.back[j,2:(df.ns.longi+1),] <- ZL.s.back[j,2:(df.ns.longi+1),] <- t(ns(f.org * t.step[j]/2 + t.step[j]/2 - 1,
                                                                                 knots = tail(head(knot.longi,-1),-1), B = knot.longi[c(1,df.ns.longi+1)]))
    }

    T.s.prep <- splineDesign(knot.surv, x = outer(f.org + 1, t.step/2), ord = 4L) # 1:qp - step[1]:f.org[1:qp]
    gk <- split(1:225, 1:qp)  # f.org group
    for (j in 1:qp) {
      T.s[,,j] <- T.s.prep[gk[[j]],]
    }

    # In addition (overall survival till different time points)
    XL.os.now <- matrix(1,qp, ncol(pool.beta))
    ZL.os.now <- matrix(1,qp, df.ns.longi + 1)
    XL.os.now[,2:(df.ns.longi+1)] <- ZL.os.now[,2:(df.ns.longi+1)] <- t(sapply(1:qp, function(n) {
      ns(f.org[n] * u/2 + u/2,
         knots = tail(head(knot.longi,-1),-1), B = knot.longi[c(1,df.ns.longi+1)])
    }))
    XL.os.now[,(df.ns.longi+2):ncol(pool.beta)] <- test.id[[longi.bsVar]]

    XL.os.back <- matrix(1,qp, ncol(pool.beta))
    ZL.os.back <- matrix(1,qp, df.ns.longi + 1)
    XL.os.back[,2:(df.ns.longi+1)] <- ZL.os.back[,2:(df.ns.longi+1)] <- t(sapply(1:qp, function(n) {
      ns(f.org[n] * u/2 + u/2 - 1,
         knots = tail(head(knot.longi,-1),-1), B = knot.longi[c(1,df.ns.longi+1)])
    }))
    XL.os.back[,(df.ns.longi+2):ncol(pool.beta)] <- test.id[[longi.bsVar]]


    T.os <- splineDesign(knot.surv, (f.org + 1) * u/2, ord = 4L)


    survival <- matrix(NA, iter, qp)
    hazard <- matrix(NA, iter, qp)
    hazard.overall.survival <- array(NA, dim = c(iter, qp, 2))

    for (j in 1:iter) {
      for (m in 1:qp) {
        survival[j,m] <- exp(
          - t.step[m]/2 * sum(w %*% exp(t(T.s[m,,1:qp]) %*% gam_bh[j,,1:2]  + rep(gamma[j,,1:2] * X, each = qp) +
                                          (t(XL.s.now[m,,1:qp]) %*% beta[j,] + t(ZL.s.now[m,,1:qp]) %*% b[j,]) %*% alpha[j,1,1:2] +
                                          (t(XL.s.now[m,,1:qp]) %*% beta[j,] + t(ZL.s.now[m,,1:qp]) %*% b[j,] -
                                             t(XL.s.back[m,,1:qp]) %*% beta[j,] - t(ZL.s.back[m,,1:qp]) %*% b[j,]) %*% alpha[j,2,1:2]))
        )

      }
      hazard.overall.survival[j,1:qp,1:2] <- exp(T.os[1:qp,] %*% gam_bh[j,,1:2] + rep(gamma[j,,1:2] * X, each = qp) +
                                                   (XL.os.now[1:qp,] %*% beta[j,] + ZL.os.now[1:qp,] %*% b[j,]) %*% alpha[j,1,1:2] +
                                                   (XL.os.now[1:qp,] %*% beta[j,] + ZL.os.now[1:qp,] %*% b[j,] -
                                                      XL.os.back[1:qp,] %*% beta[j,] - ZL.os.back[1:qp,] %*% b[j,]) %*% alpha[j,2,1:2])
      hazard[j,1:qp] <- exp(T.h[1:qp,] %*% gam_bh[j,,1] +  gamma[j,,1] * X +
                              alpha[j,1,1] * (XL.h.now[1:qp,] %*% beta[j,] + ZL.h.now[1:qp,] %*% b[j,]) +
                              alpha[j,2,1] * (XL.h.now[1:qp,] %*% beta[j,] + ZL.h.now[1:qp,] %*% b[j,] -
                                                XL.h.back[1:qp,] %*% beta[j,] - ZL.h.back[1:qp,] %*% b[j,]))
    }
    if (i==1) {
      risk.store[1:iter,i] <- (baseline.bp.risk[1:iter] + (u-t)/2 * ((hazard[1:iter,] * survival[1:iter,]) %*% w)) /
        c(surv[1:iter] + baseline.bp.risk[1:iter])
    } else {
      risk.store[1:iter,i] <- (u-t)/2 * ((hazard[1:iter,] * survival[1:iter,]) %*% w) /
        c(surv[1:iter])
    }


    overall.surv[1:iter,i] <- exp(- u/2 * (hazard.overall.survival[1:iter,,1] %*% w +
                                             hazard.overall.survival[1:iter,,2] %*% w))
    overall.surv[overall.surv[,i] < 1e-16,i] <- 1e-16
    risk[i,2] <- mean(rowSums(risk.store[,1:i, drop = F]))
    risk[i,3:4] <- quantile(rowSums(risk.store[,1:i, drop = F]),
                            probs= c(0.975, 0.025))
    full.results[[i]] <- list(rowSums(risk.store[,1:i, drop = F]))
  }

  longi[,2] <- sapply(seq(0, max(test$TimeSince_Dx), length.out = 50), function(n) {
    mean(sapply(1:iter, function(l) {
      c(1, ns(n, knots = knot.longi[2:3], B = knot.longi[c(1,4)]), test.id$DxAge) %*% beta[l,] +
        c(1, ns(n, knots = knot.longi[2:3], B = knot.longi[c(1,4)])) %*% b[l,]
    }))
  })
  return(list(risk_visits = risk[risk$t.hor %in% c(t_start, t_visits),],
              risk = risk,
              priorbiopsy_info = list(priorbsp.risk = priorbsp.risk,
                                      baseline.bp.risk = baseline.bp.risk),
              longi = longi,
              acc = meanacc.record,
              surv = surv,
              parameters = list(beta = beta,
                                b = b,
                                gam_bh = gam_bh,
                                gamma = gamma,
                                alpha = alpha,
                                Sigma = Sigma,
                                tau = tau,
                                sens = sens),
              inits = list(sigma = sigma, sigMat = sigMat),
              overall.surv = overall.surv,
              full.results = full.results))
}
