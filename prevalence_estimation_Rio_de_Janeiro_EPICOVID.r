library(rstan)
library(reshape2)
library(ggplot2)
library(cmdstanr)
### Auxiliary functions
## Some code taken from https://github.com/bob-carpenter/diagnostic-testing
stanfit <- function(fit) rstan::read_stan_csv(fit$output_files())
spin <- function(x, lower=NULL, upper=NULL, conf=0.95){
  x <- sort(as.vector(x))
  if (!is.null(lower)) {
    if (lower > min(x)) stop("lower bound is not lower than all the data")
    else x <- c(lower, x)
  }
  if (!is.null(upper)) {
    if (upper < max(x)) stop("upper bound is not higher than all the data")
    else x <- c(x, upper)
  }
  n <- length(x)
  gap <- round(conf*n)
  width <- x[(gap+1):n] - x[1:(n-gap)]
  index <- min(which(width==min(width)))
  x[c(index, index + gap)]
}
naive_est <- function(x, n, conf = 0.95){
  ## confidence interval is the normal approximation one
  phat <- x/n
  se <- sqrt(phat*(1-phat)/n)
  z <- qnorm(p = conf)
  interval <- phat + c(-1, 1)*se
  out <- c(phat, interval)
  return(out)
}
get_naive <- function(line){
  return(naive_est(x = line$x, n = line$n))
}
RoganGladen_est <- function(x, n, sens, spec, conf = 0.95){
  naive <- naive_est(x = x, n = n, conf = conf)
  denom <- sens + spec - 1
  theta.hat <- (naive[1] - (1-spec))/denom
  theta.lwr <- max(0, (naive[2] - (1-spec))/denom)
  theta.upr <- min(1, (naive[3] - (1-spec))/denom)
  out <- c(theta.hat, theta.lwr, theta.upr)
  return(out)
}
get_RG <- function(line){
  return(RoganGladen_est(x = line$x, n = line$n, sens = fixed.sens, spec = fixed.spec))
}
get_estimates_simple <- function(line){
  EPICOVID.data <- list(
    n_sample = line$n,
    x_positive = line$x,
    alpha_prev = 1,
    beta_prev = 1,
    alpha_sens =  alpha.sens,
    beta_sens =  beta.sens,
    alpha_spec =  alpha.spec,
    beta_spec = beta.spec,
    pop_size =  line$estimated_population + 0.01
  )
  mcmc.raw <- 
    simple_model$sample(data = EPICOVID.data, refresh=0,
                            parallel_chains = 4, iter_warmup = iterations,
                            iter_sampling = iterations, show_messages = FALSE)
  mcmc <- stanfit(mcmc.raw)
  summy <- summary(mcmc)$summary
  out <- summy[1:6, c(1, 4, 6, 8)]
  colnames(out) <- c("mean", "lower", "median", "upper")
  out <- cbind(out, line)
  ## Replace BCI for theta with SPI (Liu et al, 2015)
  sp.theta <- spin(extract(mcmc, 'theta')$theta, lower = 0, upper = 1)
  summy.theta <- out["theta",]
  summy.theta[c(2, 4)] <-   sp.theta
  out["theta_hpd",] <- summy.theta
  ##
  out$parameter <- rownames(out)
  return(out)
}
get_estimates_detection <- function(line){
  EPICOVID.data <- list(
    n_sample = line$n,
    x_positive = line$x,
    alpha_prev = 1,
    beta_prev = 1,
    alpha_detect = alpha.detect,
    beta_detect = beta.detect,
    alpha_sens =  alpha.sens,
    beta_sens =  beta.sens,
    alpha_spec =  alpha.spec,
    beta_spec = beta.spec,
    Y_obs = line$confirmed + 0.01,
    pop_size = line$estimated_population + 0.01
  )
  mcmc.raw <- 
    detection_model$sample(data = EPICOVID.data, refresh=0,
                            parallel_chains = 4, iter_warmup = iterations,
                            iter_sampling = iterations, show_messages = FALSE)
  mcmc <- stanfit(mcmc.raw)
  summy <- summary(mcmc)$summary
  out <- summy[1:7, c(1, 4, 6, 8)]
  colnames(out) <- c("mean", "lower", "median", "upper")
  out <- cbind(out, line)
  ## Replace BCI for theta with SPI (Liu et al, 2015)
  sp.theta <- spin(extract(mcmc, 'theta')$theta, lower = 0, upper = 1)
  summy.theta <- out["theta", ]
  summy.theta[c(2, 4)] <-   sp.theta
  out["theta_hpd",] <- summy.theta
  ##
  out$parameter <- rownames(out)
  return(out)
}
elicit_beta_quantiles <- function(l, u, alpha = .95){
  q0 <- (1-alpha)/2
  q1 <- (1 + alpha)/2
  loss <- function(x){
    a <- x[1]
    b <- x[2]
    l.hat <- qbeta(q0, shape1 = a, shape2 = b)
    u.hat <- qbeta(q1, shape1 = a, shape2 = b)
    # error <- (l.hat - l)^2 + (u.hat-u)^2
    error <- abs(l.hat - l) + abs(u.hat-u) ## L1 norm: better?
    return(error)  
  }
  opt <- suppressWarnings( optim(loss, par = c(1, 1), lower = c(1E-3, 1E-3)) )
  a <- opt$par[1]
  b <- opt$par[2]
  if(a < 0 || b < 0) warning("Warning: at least one of the obtained parameters is not valid")
  return(list(a = a, b = b))
}
############################################
##### Preparation
## Models 
iterations <- 500
simple_model <- cmdstan_model("stan/simple_prevalence_y_correction.stan")
detection_model <- cmdstan_model("stan/simple_prevalence_y_correction_ProbDetection.stan")

## Load data
complete.data <- read.csv("data/EPICOVID_aggregated_complete.csv")
RJ <- subset(complete.data, city == "BREVES")
M <- RJ$estimated_population[1]

alpha.sens <- 312 ## alpha_delta
beta.sens <-  49  ## beta_delta
alpha.spec <- 234 ## alpha_gamma
beta.spec <-  1 ## beta_gamma

fixed.sens <- alpha.sens/(alpha.sens + beta.sens)
fixed.spec <- alpha.spec/(alpha.spec + beta.spec)

### Probability of detection parameters
under.l <- 1/20
under.u <- 1/3
pdetect.pars <- elicit_beta_quantiles(l = under.l, u = under.u)
pdetect.pars$a/(pdetect.pars$a + pdetect.pars$b)

alpha.detect <- pdetect.pars$a
beta.detect <- pdetect.pars$b
qbeta(p = c(.025, .975), shape1 = alpha.detect, shape2 = beta.detect) ## how successful was the elicitation procedure?

##### Running

## Naive
naive.ests <- lapply(1:nrow(RJ), function(i){
  out <- get_naive(RJ[i, ])
  return(out)
} )
lapply(naive.ests, round, 3)

## Rogan-Gladen
RG.ests <- lapply(1:nrow(RJ), function(i){
  out <- get_RG(RJ[i, ])
  return(out)
} )
lapply(RG.ests, round, 3)

## Simple full Bayes model, theta (prevalence)

simple.fits <- lapply(1:nrow(RJ), function(i){
  sink("/dev/null")
  out <- get_estimates_simple(RJ[i, ])
  sink()
  return(out)
} )

simple.results.dt <- do.call(rbind, simple.fits)
rownames(simple.results.dt) <- NULL

theta.simple <- subset(simple.results.dt, parameter == "theta")
theta_hpd.simple <- subset(simple.results.dt, parameter == "theta_hpd")

round(theta.simple[, c(1, 2, 4)], 3)
round(theta_hpd.simple[, c(3, 2, 4)], 3)


## ProbDetection full Bayes model, theta (prevalence)

detection.fits <- lapply(1:nrow(RJ), function(i){
  sink("/dev/null")
  out <- get_estimates_detection(RJ[i, ])
  sink()
  return(out)
} )

detection.results.dt <- do.call(rbind, detection.fits)
rownames(detection.results.dt) <- NULL

theta_hpd.detection <- subset(detection.results.dt, parameter == "theta_hpd")

round(theta_hpd.detection[, c(3, 2, 4)], 3)

p.detection <- subset(detection.results.dt, parameter == "p_detect")

round(p.detection[, c(1, 2, 4)], 2)

underreporting <- subset(detection.results.dt, parameter == "F")

round(underreporting[, c(1, 2, 4)], 3)

### Now let's get the cases

lapply(naive.ests, function(x) round(x*M))

lapply(RG.ests, function(x) round(x*M))


N.simple <- subset(simple.results.dt, parameter == "N")
round(N.simple[, c(1, 2, 4)])

N.detection <- subset(detection.results.dt, parameter == "N")
round(N.detection[, c(1, 2, 4)])

