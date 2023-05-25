# Extract posterior distribution from JAGS ----------------------------

extract_pd <- function(object = NULL, sample_pool_size = 500) {

  if (is.null(object)) {
    stop("JAGS results of posterior distribution is missing!")
  }

  object <- lapply(object, as.matrix)
  col.name <- colnames(object[[1]])

  if (nrow(object[[1]]) < sample_pool_size) {
    stop("No enough samples from Posterior distribution!")
  }

  idx_pd <- (nrow(object[[1]]) - sample_pool_size + 1):nrow(object[[1]])

  # baseline hazard coefficients
  pool.bs_gammas <- do.call("rbind", lapply(object, function(x) {x[idx_pd,
                                                                   grep("gambh", col.name)]}))

  pool.gammas <- do.call("rbind", lapply(object, function(x) {x[idx_pd,
                                                                grep("gamma", col.name)]}))

  pool.alphas <- do.call("rbind", lapply(object, function(x) {x[idx_pd,
                                                                grep("alpha", col.name)]}))

  pool.beta <- do.call("rbind", lapply(object, function(x) {x[idx_pd,
                                                              grep("betaL", col.name)]}))

  pool.D <- do.call("rbind", lapply(object, function(x) {x[idx_pd,
                                                           grep("D", col.name)]}))

  pool.tau <- do.call("c", lapply(object, function(x) {x[idx_pd,
                                                         grep("tau", col.name)]}))

  return(list(pool.bs_gammas = pool.bs_gammas,
              pool.gammas = pool.gammas,
              pool.alphas = pool.alphas,
              pool.beta = pool.beta,
              pool.D = pool.D,
              pool.tau = pool.tau))
}
