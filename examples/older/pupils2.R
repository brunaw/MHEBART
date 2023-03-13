# 4 different ways to fit a cross random effects model
# 1. Using lme4
# 2. Using JAGS
# 3. Using raw R code for MCMC
# 4. Using the MHEBART maths and simplifying to stumps and one tree
# All of these should agree (especially 2-4 as they're Bayesian)
# If not we're in trouble

# Packages
library(tidyverse)
library(lme4)
library(R2jags)
library(mvnfast)
library(LaplacesDemon) # For half t prior
devtools::load_all(".")

load("data/pupils.Rdata")
# pupils


# 1. lmer -----------------------------------------------------------------

pupils_lmer <- lmer(
  achievement ~ 1
    + (1 | primary_school_id)
    + (1 | secondary_school_id),
  data = pupils
)
# View(pupils)
pplme <- predict(pupils_lmer, pupils)
rmse_lmer <- sqrt(mean((pplme - pupils$achievement)^2)) # 6.64
rmse_lmer
summary(pupils_lmer)

# pupils |>
#   mutate(pred = pplme) |>
#   group_by(primary_school_id, secondary_school_id) |>
#   summarise(y = mean(achievement),
#             pred = mean(pred),
#             n = n()) |>
#   View()

# Random effects:
#   Groups              Name        Variance Std.Dev.
# primary_school_id   (Intercept) 0.17190  0.4146
# secondary_school_id (Intercept) 0.06666  0.2582
# Residual                        0.51313  0.7163
# Number of obs: 1000, groups:
#   primary_school_id, 50; secondary_school_id, 30
#
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)   6.3486     0.0789   80.46


# 2. Using JAGS -----------------------------------------------------------

# Define the model
model_string <- "
model {
  for (i in 1:N) {
    y[i] ~ dnorm(fit[i], sigma^-2)
    fit[i] = mu + phi1[group1[i]] + phi2[group2[i]]
  }
  mu ~ dnorm(0, 10^-2)
  for (j in 1:len_group1) {
    phi1[j] ~ dnorm(0, sigma1^-2)
  }
  for (j in 1:len_group2) {
    phi2[j] ~ dnorm(0, sigma2^-2)
  }
  sigma ~ dt(0, 10^-2, 1)T(0,)
  sigma1 ~ dt(0, 10^-2, 1)T(0,)
  sigma2 ~ dt(0, 10^-2, 1)T(0,)
}"

# Compile the model
model_data <- list(
  y = pupils$achievement,
  group1 = pupils$primary_school_id,
  group2 = pupils$secondary_school_id,
  N = nrow(pupils),
  len_group1 = length(unique(pupils$primary_school_id)),
  len_group2 = length(unique(pupils$secondary_school_id))
)

# Run the model
model <- jags(
  model.file = textConnection(model_string),
  data = model_data,
  parameters.to.save = c(
    "mu", "phi1", "phi2", "fit",
    "sigma1", "sigma2", "sigma"
  )
)

# plot(model)

# pupils |>
#   mutate(pred_JAGS = as.numeric(model$BUGSoutput$mean$fit),
#          pred_lmer = pplme) |>
#   group_by(primary_school_id, secondary_school_id) |>
#   summarise(y = mean(achievement),
#             pred_JAGS = mean(pred_JAGS),
#             pred_lmer = mean(pred_lmer),
#             n = n()) |>
#   View() # Almost identical fits


# 3. MCMC code ------------------------------------------------------------

# Model is y_i ~ N(fit[i], tau^-1)
# fit[i] = mu + b[group1[i]] + c[group2[i]]
# b[j] ~ N(0, tau_b^-2)
# c[k] ~ N(0, tau_c^-2)
# tau_b, tau_c ~ dt(0, 10, 1)+
# mu ~ N(0, a^-1)
# sigma ~ Ga(d, e)

# Let M1 and M2 be the group allocation matrices
# Can write as y ~ N(mu%*%ones + M1%*%b + M2%*%c, I/tau)

# So need parameter updates for mu, b, c, sigma (all Gibbs)
# And tau_b, tau_c both MH

# Set up everything
set.seed(123)
num_iter <- 1000
y <- pupils$achievement
# Scale the response target variable
y_min <- min(y)
y_max <- max(y)
y_scale <- (y - y_min)/(y_max - y_min) - 0.5
y = y_scale
M1 <- stats::model.matrix(~ as.factor(pupils$primary_school_id) - 1)
M2 <- stats::model.matrix(~ as.factor(pupils$secondary_school_id) - 1)

# Some useful things needed later
tM1xM1 <- crossprod(M1)
tM2xM2 <- crossprod(M2)
n <- length(y)

# Hyper-parameter values
a <- d <- e <- 0.01

# Starting values
tau <- tau_b <- tau_c <- 1
b <- matrix(rep(0, ncol(M1)), ncol = 1)
c <- matrix(rep(0, ncol(M2)), ncol = 1)
mu <- 0

# Storage
storage <- list(
  mu = rep(NA, num_iter),
  b = matrix(NA, ncol = ncol(M1), nrow = num_iter),
  c = matrix(NA, ncol = ncol(M2), nrow = num_iter),
  tau = rep(NA, num_iter),
  tau_b = rep(NA, num_iter),
  tau_c = rep(NA, num_iter),
  fits = matrix(NA, ncol = n, nrow = num_iter)
)

# Stuff for MH
tau_b_sd <- tau_c_sd <- 1

# Progress bar
pb <- utils::txtProgressBar(
  min = 1, max = num_iter,
  style = 3, width = 60,
  title = "Running model..."
)

for (i in 1:num_iter) {
  utils::setTxtProgressBar(pb, i)
  
  # Storage update
  storage$mu[i] <- mu
  storage$b[i,] <- b
  storage$c[i,] <- c
  storage$tau[i] <- tau
  storage$tau_b[i] <- tau_b
  storage$tau_c[i] <- tau_c
  storage$fits[i,] <- mu + M1 %*% b + M2 %*% c
  
  # Update mu
  Rmu <- y - M1 %*% b - M2 %*% c
  prec <- tau * n + a
  mu <- rnorm(1, sum(Rmu) / prec, sqrt(1 / prec))

  # Update b
  Rb <- y - mu - M2 %*% c
  prec <- tau * tM1xM1 + diag(tau_b, ncol(M1))
  b <- t(mvnfast::rmvn(1, solve(prec, tau * crossprod(M1, Rb)), solve(prec)))

  # Update c
  Rc <- y - mu - M1 %*% b
  prec <- tau * tM2xM2 + diag(tau_c, ncol(M2))
  c <- t(mvnfast::rmvn(1, solve(prec, tau * crossprod(M2, Rc)), solve(prec)))

  # Update tau
  S <- sum((y - mu - M1 %*% b - M2 %*% c)^2)
  tau <- rgamma(1,
    shape = d + n / 2,
    rate = e + S / 2
  )

  # Update tau_b
  repeat {
    # Proposal distribution
    new_tau_b <- tau_b + stats::rnorm(1, sd = tau_b_sd)
    if (new_tau_b > 0) {
      break
    }
  }
  log_rat <- stats::pnorm(tau_b, sd = tau_b_sd, log = TRUE) -
    stats::pnorm(new_tau_b, sd = tau_b_sd, log = TRUE)

  post_new <- sum(dnorm(b, 0, 1/sqrt(new_tau_b))) + dhalft(new_tau_b, scale = 10, nu = 1, log = TRUE)
  post_old <- sum(dnorm(b, 0, 1/sqrt(tau_b))) + dhalft(tau_b, scale = 10, nu = 1, log = TRUE)

  log_alpha <- post_new - post_old + log_rat

  accept <- log_alpha >= 0 || log_alpha >= log(stats::runif(1))
  tau_b <- ifelse(accept, new_tau_b, tau_b)

  # Update tau_c
  repeat {
    # Proposal distribution
    new_tau_c <- tau_c + stats::rnorm(1, sd = tau_c_sd)
    if (new_tau_c > 0) {
      break
    }
  }
  log_rat <- stats::pnorm(tau_c, sd = tau_c_sd, log = TRUE) -
    stats::pnorm(new_tau_c, sd = tau_c_sd, log = TRUE)
  
  post_new <- sum(dnorm(b, 0, 1/sqrt(new_tau_c))) + LaplacesDemon::dhalft(new_tau_c, scale = 10, nu = 1, log = TRUE)
  post_old <- sum(dnorm(b, 0, 1/sqrt(tau_c))) + LaplacesDemon::dhalft(tau_c, scale = 10, nu = 1, log = TRUE)
  
  log_alpha <- post_new - post_old + log_rat
  
  accept <- log_alpha >= 0 || log_alpha >= log(stats::runif(1))
  tau_c <- ifelse(accept, new_tau_c, tau_c)
}

# Plot some of the outputs
# plot(1/sqrt(storage$tau))
# plot(storage$mu)
# plot(1/sqrt(storage$tau_b))
# plot(1/sqrt(storage$tau_c))

# Compare models ----------------------------------------------------------

# pupils |>
#   mutate(pred_JAGS = as.numeric(model$BUGSoutput$mean$fit),
#          pred_MCMC = colMeans(storage$fits),
#          pred_lmer = pplme) |>
#   group_by(primary_school_id, secondary_school_id) |>
#   summarise(y = mean(achievement),
#             pred_JAGS = mean(pred_JAGS),
#             pred_MCMC = mean(pred_MCMC),
#             pred_lmer = mean(pred_lmer),
#             n = n()) |>
#   View() # Almost identical fits

sqrt(mean((pplme - pupils$achievement)^2)) # 0.6921818
sqrt(mean((model$BUGSoutput$mean$fit - pupils$achievement)^2)) # 0.6920755
sqrt(mean((colMeans(storage$fits) - pupils$achievement)^2)) # 0.6901241


