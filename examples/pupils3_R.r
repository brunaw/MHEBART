# 3 different ways to fit a cross random effects model
# This model standardises the data so uses no mean effects
# 1. Using lme4
# 2. Using JAGS
# 3. Using raw R code for MCMC
# All of these should agree (especially 2 and 3 as they're Bayesian)
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

# Standardise the pupils data

# 1. lmer -----------------------------------------------------------------

pupils_lmer <- lmer(
  achievement ~ -1
    + (1 | primary_school_id)
    + (1 | secondary_school_id),
  data = pupils
)

pplme <- predict(pupils_lmer, pupils)
summary(pupils_lmer)

# pupils |>
#   mutate(pred = pplme) |>
#   group_by(primary_school_id, secondary_school_id) |>
#   summarise(y = mean(achievement),
#             pred = mean(pred),
#             n = n()) |>
#   View()

# 2. Using JAGS -----------------------------------------------------------

# Define the model
model_string <- "
model {
  for (i in 1:N) {
    y[i] ~ dnorm(fit[i], tau)
    fit[i] = tree1[group1[i]] + tree2[group2[i]]
  }
  mu ~ dnorm(0, 10^-2)
  for (j in 1:len_group1) {
    tree1[j] ~ dnorm(0, sigma1^-2)
  }
  for (j in 1:len_group2) {
    tree2[j] ~ dnorm(0, sigma2^-2)
  }
  tau ~ dgamma(0.01, 0.01)
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
    "mu", "tree1", "tree2", "fit",
    "sigma1", "sigma2", "tau"
  )
)

# plot(model)
# 
# pupils |>
#   mutate(pred_JAGS = as.numeric(model$BUGSoutput$mean$fit),
#          pred_lmer = pplme) |>
#   group_by(primary_school_id, secondary_school_id) |>
#   summarise(y = mean(achievement),
#             pred_JAGS = mean(pred_JAGS),
#             pred_lmer = mean(pred_lmer),
#             n = n()) |>
#   View() # Almost identical fits
# stop()

# 3. MCMC code ------------------------------------------------------------

# Model is y_i ~ N(fit[i], tau^-1)
# fit[i] = tree1[group1[i]] + tree2[group2[i]]
# tree1[j] ~ N(0, sigma1^2)
# tree2[k] ~ N(0, sigma2^2)
# sigma1, sigma2 ~ dt(0, 10, 1)+
# tau ~ Ga(0.01, 0.01)

# Let M1 and M2 be the group allocation matrices
# Can write as y ~ N(M1%*%tree1 + M2%*%tree2, I/tau)

# So need parameter updates for tree1, tree2, tau (all Gibbs)
# And sigma1, sigma2 both MH

y <- pupils$achievement
y_min <- min(y)
y_max <- max(y)
y_scale <- (y - y_min)/(y_max - y_min) - 0.5
n <- length(y_scale)
y = y_scale


# Set up everything
set.seed(123)
num_iter <- 50
#y <- pupils$achievement
M1 <- stats::model.matrix(~ as.factor(pupils$primary_school_id) - 1)
M2 <- stats::model.matrix(~ as.factor(pupils$secondary_school_id) - 1)

# Some useful things needed later
tM1xM1 <- crossprod(M1)
tM2xM2 <- crossprod(M2)
n <- length(y)

# Starting values
tau <- 1
sigma1 <- sigma2 <- 1
tree1 <- matrix(rep(0, ncol(M1)), ncol = 1)
tree2 <- matrix(rep(0, ncol(M2)), ncol = 1)

# Storage
storage <- list(
  tree1 = matrix(NA, ncol = ncol(M1), nrow = num_iter),
  tree2 = matrix(NA, ncol = ncol(M2), nrow = num_iter),
  tau = rep(NA, num_iter),
  sigma1 = rep(NA, num_iter),
  sigma2 = rep(NA, num_iter),
  fits = matrix(NA, ncol = n, nrow = num_iter)
)

# Stuff for MH
sigma1_sd <- sigma2_sd <- 1

# Progress bar
pb <- utils::txtProgressBar(
  min = 1, max = num_iter,
  style = 3, width = 60,
  title = "Running model..."
)

for (i in 1:num_iter) {
  utils::setTxtProgressBar(pb, i)
  
  # Storage update
  storage$tree1[i,] <- tree1
  storage$tree2[i,] <- tree2
  storage$tau[i] <- tau
  storage$sigma1[i] <- sigma1
  storage$sigma2[i] <- sigma2
  storage$fits[i,] <- M1 %*% tree1 + M2 %*% tree2

  # Update tree1
  Rtree1 <- y - M2 %*% tree2
  
  prec <- tau * tM1xM1 + diag(1/(sigma1^2), ncol(M1))
  tree1 <- t(mvnfast::rmvn(1, solve(prec, tau * crossprod(M1, Rtree1)), 
                           solve(prec)))
  
  # # Update tree2
  Rtree2 <- y - M1 %*% tree1
  
  prec <- tau * tM2xM2 + diag(1/(sigma2^2), ncol(M2))
  tree2 <- t(mvnfast::rmvn(1, 
                           solve(prec, tau * crossprod(M2, Rtree2)), 
                           solve(prec)))
  
# # Update tau
  S <- sum((y - M1 %*% tree1 - M2 %*% tree2)^2)
  tau <- rgamma(1,
    shape = 0.01 + n / 2,
    rate = 0.01 + S / 2
  )
  # Update sigma1
  repeat {
    # Proposal distribution
    new_sigma1 <- sigma1 + stats::rnorm(1, sd = sigma1_sd)
    if (new_sigma1 > 0) {
      break
    }
  }
  log_rat <- stats::pnorm(sigma1, sd = sigma1_sd, log = TRUE) -
    stats::pnorm(new_sigma1, sd = sigma1_sd, log = TRUE)

  post_new <- sum(dnorm(tree1, 0, new_sigma1, log = TRUE)) + dhalft(new_sigma1, scale = 10, nu = 1, log = TRUE)
  post_old <- sum(dnorm(tree1, 0, sigma1, log = TRUE)) + dhalft(sigma1, scale = 10, nu = 1, log = TRUE)

  log_alpha <- post_new - post_old + log_rat

  accept <- log_alpha >= 0 || log_alpha >= log(stats::runif(1))
  sigma1 <- ifelse(accept, new_sigma1, sigma1)
  # Update sigma2
  repeat {
    # Proposal distribution
    new_sigma2 <- sigma2 + stats::rnorm(1, sd = sigma2_sd)
    if (new_sigma2 > 0) {
      break
    }
  }
  log_rat <- stats::pnorm(sigma2, sd = sigma2_sd, log = TRUE) -
    stats::pnorm(new_sigma2, sd = sigma2_sd, log = TRUE)

  post_new <- sum(dnorm(tree2, 0, new_sigma2, log = TRUE)) + LaplacesDemon::dhalft(new_sigma2, scale = 10, nu = 1, log = TRUE)
  post_old <- sum(dnorm(tree2, 0, sigma2, log = TRUE)) + LaplacesDemon::dhalft(sigma2, scale = 10, nu = 1, log = TRUE)

  log_alpha <- post_new - post_old + log_rat

  accept <- log_alpha >= 0 || log_alpha >= log(stats::runif(1))
  sigma2 <- ifelse(accept, new_sigma2, sigma2)
  sigma2 <- sigma1 <- tau <- 1
}

# Plot some of the outputs
plot(1/sqrt(storage$tau))
plot(storage$sigma1)
plot(storage$sigma2)

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

mean((pplme - pupils$achievement)^2) # 0.7978818
mean((model$BUGSoutput$mean$fit - pupils$achievement)^2) # 0.797759
mean((inv_scale(colMeans(storage$fits), y_max, y_min) - pupils$achievement)^2) # 0.7974181
# pretty much the same

inv_scale <- function(y, max, min) (y + 0.5) * (max - min) + min


# I've checked that with all things fixed 
# we have the same results