# Create a really simple example that HEBART should do really well on

# Just one tree with one covariate and 3 groups
# Terminal node 1 comes from x <= 0.5, has mu = -5, phi 1 2 3 = -7 -5 -3
# Node 2 comes from x > 0.5 and x <= 0.75, has mu = 10, phi = 8 10 12
# Node 3 comes from x > 0.75 has mu  0 and phi = -2 0 2

# Clear workspace and load in packages
library(tidyverse)
devtools::load_all(".")

# Create this data set
simulate_tree <- function(x){
  #set.seed(123)
  n <- 400
  x <- runif(n)
  groups <- factor(rep(1:3, length = n))
  mu <- rep(0, n)
  mu[x<0.75] <- 10
  mu[x<0.5] <- -5
  phi <- rep(NA, length = n)
  phi_tn1 <- c(-7, -5, -3)
  phi[mu == -5] <- phi_tn1[groups[mu == -5]]
  phi_tn2 <- c(8, 10, 12)
  phi[mu == 10] <- phi_tn2[groups[mu == 10]]
  phi_tn3 <- c(-2, 0, 2)
  phi[mu == 0] <- phi_tn3[groups[mu == 0]]
  res_sd <- 1
  y <- rnorm(n, phi, res_sd)
  
  train <- data.frame(
    y = y[1:200], 
    x = x[1:200],
    groups = groups[1:200]
  )
  test <- data.frame(
    y = y[201:400], 
    x = x[201:400],
    groups = groups[201:400]
  )
  
  
  #set.seed(321)
  x <- runif(n)
  groups <- factor(rep(1:3, length = n))
  mu <- rep(0, n)
  mu[x<0.85] <- -3
  mu[x<0.4] <- 8
  phi <- rep(NA, length = n)
  phi_tn1 <- c(7, -2, -10)
  phi[mu == -3] <- phi_tn1[groups[mu == -3]]
  phi_tn2 <- c(8, 15, 7)
  phi[mu == 8] <- phi_tn2[groups[mu == 8]]
  phi_tn3 <- c(-2, 0, 2)
  phi[mu == 0] <- phi_tn3[groups[mu == 0]]
  # qplot(x, phi, colour = groups)
  res_sd <- 1
  y <- rnorm(n, phi, res_sd)
  
  
  train2 <- data.frame(
    y = y[1:200], 
    x = x[1:200],
    groups = groups[1:200]
  )
  test2 <- data.frame(
    y = y[201:400], 
    x = x[201:400],
    groups = groups[201:400]
  )
  
  train3 <- data.frame(
    y = train$y + train2$y, 
    group_1 = train$groups, 
    group_2 = train2$groups, 
    x1 = train$x,
    x2 = train2$x
  )
  
  test3 <- data.frame(
    y = test$y + test2$y, 
    group_1 = test$groups, 
    group_2 = test2$groups, 
    x1 = test$x,
    x2 = test2$x
  )
  
  return(list(train = train3, test = test3))
}
  
data <- tibble(
  data = map(1:10, simulate_tree), 
  train = map(data, "train"),
  test = map(data, "test"))


# Now run HEBART on it ----------------------------------------------------

run_model <- function(train){
  num_trees <- 20
  
  hb_model <- mhebart(
    formula = y ~ x1 + x2,
    data = train,
    group_variables = c("group_1", "group_2"), 
    num_trees = num_trees,
    priors = list(
      alpha = 0.95, # Prior control list
      beta = 2,
      nu = 2,
      lambda = 0.1,
      tau_mu = 16 * num_trees,
      shape_sigma_phi = 0.5,
      scale_sigma_phi = 1, 
      sample_sigma_phi = TRUE
    ), 
    inits = list(tau = 1,
                 sigma_phi = 0.01),
    MCMC = list(iter = 1500, 
                burn = 500, 
                thin = 1,
                sigma_phi_sd = 1)
  )
  hb_model
}

data_model <- data |> 
  mutate(model = map(train, run_model))

predict_test <- function(model, test){
  pp <- predict_mhebart(newX = test, 
                        group_variables = c("group_1", "group_2"),
                        hebart_posterior  = model, 
                        type = "mean")
  sqrt(mean((pp - test$y)^2))
}

data_model2 <- data_model |> 
  mutate(rmse = map2_dbl(model, test, predict_test))

#cor(pp, test$y)
#qplot(test$y, pp) + geom_abline()

predict_lme <- function(train, test){
  lme_ss <- lme4::lmer(y ~ x1 + x2 + (1|group_1) + (1|group_2), train)
  pplme <- predict(lme_ss, test)
  rmse_lmer <- sqrt(mean((pplme - test$y)^2)) # 6.64
  rmse_lmer
}

data_model3 <- data_model2 |> 
  mutate(rmse_lme = map2_dbl(train, test, predict_lme))

summary(data_model3$rmse)
summary(data_model3$rmse_lme)
#----------------------------------------------------------------------
#----------------------------------------------------------------------

cor(pplme, test$y) # 0.936175
rmse_lmer
