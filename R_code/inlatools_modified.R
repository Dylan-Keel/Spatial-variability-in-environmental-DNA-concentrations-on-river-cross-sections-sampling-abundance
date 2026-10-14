## Modified copies of dispersion_check() and distribution_check() from the R
## package inlatools (Thierry Onkelinx, INBO; GPL-3,
## https://github.com/inbo/inlatools).
##
## The package versions read the observed response out of the model object.
## That does not work for models fit through an inla.stack, so these copies
## take the observed counts as an argument instead. Nothing else is changed;
## plotting uses the package's own plot methods.
## The fast_distribution_check() copy that was in the thesis code is not used
## and is not kept here.

library(assertthat)
library(dplyr)
library(tidyr)
library(purrr)

dispersion <- function(observed, fitted, variance) {
  mean((observed - fitted) ^ 2 / variance)
}

# simulate from the fitted model and compare the dispersion of the simulations
# with the dispersion of the data
dispersion_check <- function(object, observed, nsim = 1000) {
  assert_that(is.count(nsim))
  if (length(object$.args$family) > 1) stop("Only single responses are handled")
  mu <- fitted(object)[seq_along(observed)]
  which_na <- is.na(observed)
  mu <- mu[!which_na]
  observed <- observed[!which_na]
  if (object$.args$family == "poisson") {
    dispersion_data <- dispersion(observed = observed, fitted = mu, variance = mu)
    dispersion_model <- apply(
      matrix(rpois(nsim * length(mu), lambda = mu), ncol = nsim), 2,
      function(x) dispersion(observed = x, fitted = mu, variance = mu))
  } else if (object$.args$family == "nbinomial") {
    size <- object$summary.hyperpar[
      paste("size for the nbinomial observations (1", "overdispersion)", sep = "/"), "mean"]
    dispersion_data <- dispersion(observed = observed, fitted = mu, variance = mu + mu ^ 2 / size)
    dispersion_model <- apply(
      matrix(rnbinom(nsim * length(mu), mu = mu, size = size), ncol = nsim), 2,
      function(x) dispersion(observed = x, fitted = mu, variance = mu + mu ^ 2 / size))
  } else if (object$.args$family == "zeroinflatednbinomial1") {
    mu <- exp(object$summary.linear.predictor[!which_na, "0.5quant"])
    size <- object$summary.hyperpar["size for nbinomial zero-inflated observations", "0.5quant"]
    zi <- object$summary.hyperpar["zero-probability parameter for zero-inflated nbinomial_1", "0.5quant"]
    zi_mu <- (1 - zi) * mu
    zi_var <- (1 - zi) * mu * (1 + mu * (zi + 1 / size))
    dispersion_data <- dispersion(observed = observed, fitted = zi_mu, variance = zi_var)
    dispersion_model <- apply(
      matrix(rbinom(nsim * length(mu), size = 1, prob = 1 - zi) *
               rnbinom(nsim * length(mu), mu = mu, size = size), ncol = nsim), 2,
      function(x) dispersion(observed = x, fitted = zi_mu, variance = zi_var))
  } else if (object$.args$family == "zeroinflatedpoisson1") {
    lambda <- exp(object$summary.linear.predictor[!which_na, "0.5quant"])
    zi <- object$summary.hyperpar[1, "0.5quant"]
    zi_mu  <- (1 - zi) * lambda
    zi_var <- (1 - zi) * (lambda ^ 2 + lambda) - zi_mu ^ 2
    dispersion_data <- dispersion(observed = observed, fitted = zi_mu, variance = zi_var)
    dispersion_model <- apply(
      matrix(rbinom(nsim * length(lambda), size = 1, prob = 1 - zi) *
               rpois(nsim * length(lambda), lambda = lambda), ncol = nsim), 2,
      function(x) dispersion(observed = x, fitted = zi_mu, variance = zi_var))
  } else {
    stop(object$.args$family, " is not yet handled")
  }
  output <- list(data = dispersion_data, model = dispersion_model)
  class(output) <- "dispersion_check"
  output
}


# simulate the response from the posterior and compare the empirical
# distribution of the simulations with that of the data
distribution_check <- function(object, observed, nsim = 1000, seed = 0L) {
  assert_that(is.count(nsim))
  if (length(object$.args$family) > 1) stop("Only single responses are handled")
  samples <- INLA::inla.posterior.sample(n = nsim, result = object, seed = seed)
  relevant <- grep("^Predictor:", rownames(samples[[1]]$latent))
  eta <- map_dfc(samples, "latent")[relevant, , drop = FALSE]
  n_sampled <- switch(
    object$.args$family,
    poisson = {
      mutate_all(exp(eta), rpois, n = nrow(eta)) %>%
        gather("run", "x") %>% count(.data$run, .data$x)
    },
    nbinomial = {
      relevant <- grep("size for the nbinomial observations", names(samples[[1]]$hyperpar))
      size <- INLA::inla.hyperpar.sample(n = nsim, result = object)[, relevant]
      mutate_all(exp(eta), function(mu) rnbinom(n = length(mu), mu = mu, size = size)) %>%
        gather("run", "x") %>% count(.data$run, .data$x)
    },
    zeroinflatedpoisson1 = {
      relevant <- grep("zero-inflated poisson_1", names(samples[[1]]$hyperpar))
      mutate_all(exp(eta), function(lambda, prob_zero) {
        n <- length(lambda)
        rbinom(n = n, size = 1, prob = 1 - prob_zero) * rpois(n = n, lambda = lambda)
      }, prob_zero = INLA::inla.hyperpar.sample(n = nsim, result = object)[, relevant]) %>%
        gather("run", "x") %>% count(.data$run, .data$x)
    },
    stop(object$.args$family, " is not yet handled")
  )
  n_observed <- data.frame(x = observed) %>% count(.data$x)
  n_count <- unique(c(n_observed$x, n_sampled$x))
  ecdf <- n_sampled %>%
    complete(run = .data$run, x = n_count, fill = list(n = 0)) %>%
    group_by(.data$run) %>% arrange(.data$x) %>%
    mutate(ecdf = cumsum(.data$n) / sum(.data$n)) %>%
    group_by(.data$x) %>%
    summarise(median = quantile(.data$ecdf, probs = 0.5),
              lcl    = quantile(.data$ecdf, probs = 0.025),
              ucl    = quantile(.data$ecdf, probs = 0.975)) %>%
    inner_join(n_observed %>% complete(x = n_count, fill = list(n = 0)) %>%
                 arrange(.data$x) %>% mutate(ecdf = cumsum(.data$n) / sum(.data$n)),
               by = "x")
  class(ecdf) <- c("distribution_check", class(ecdf))
  ecdf
}

