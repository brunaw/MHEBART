# Create a really simple example that HEBART should do really well on

# Just one tree with one covariate and 3 groups
# Terminal node 1 comes from x <= 0.5, has mu = -5, phi 1 2 3 = -7 -5 -3
# Node 2 comes from x > 0.5 and x <= 0.75, has mu = 10, phi = 8 10 12
# Node 3 comes from x > 0.75 has mu  0 and phi = -2 0 2

# Clear workspace and load in packages
library(tidyverse)
library(firatheme)
devtools::load_all(".")

# Create this data set
simulate_tree <- function(x, n = 500){
  #set.seed(123)
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
    y = y[1:(n/2)], 
    x = x[1:(n/2)],
    groups = groups[1:(n/2)]
  )
  test <- data.frame(
    y = y[(n/2+1):n], 
    x = x[(n/2+1):n],
    groups = groups[(n/2+1):n]
  )
  
  
  #set.seed(321)
  x <- runif(n)
  groups <- factor(rep(1:5, length = n))
  mu <- rep(0, n)
  mu[x<0.85] <- -3
  mu[x<0.4] <- 8
  phi <- rep(NA, length = n)
  phi_tn1 <- c(7, -2, -10, 2, 5)
  phi[mu == -3] <- phi_tn1[groups[mu == -3]]
  phi_tn2 <- c(8, 15, 7, 4, -1)
  phi[mu == 8] <- phi_tn2[groups[mu == 8]]
  phi_tn3 <- c(-2, 0, 2, 20, 4)
  phi[mu == 0] <- phi_tn3[groups[mu == 0]]
  # qplot(x, phi, colour = groups)
  res_sd <- 1
  y <- rnorm(n, phi, res_sd)
  
  
  train2 <- data.frame(
    y = y[1:(n/2)], 
    x = x[1:(n/2)],
    groups = groups[1:(n/2)]
  )
  test2 <- data.frame(
    y = y[(n/2+1):n], 
    x = x[(n/2+1):n],
    groups = groups[(n/2+1):n]
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
  rmse <- sqrt(mean((pp - test$y)^2))
  return(list(pred = pp, rmse = rmse))
}

data_model2 <- data_model |> 
  mutate(mh = map2(model, test, predict_test), 
         pred_mh = map(mh, "pred"), 
         rmse_hm = map_dbl(mh, "rmse"))

#cor(pp, test$y)
#qplot(test$y, pp) + geom_abline()

predict_lme <- function(train, test){
  lme_ss <- lme4::lmer(y ~ x1 + x2 + (1|group_1) + (1|group_2), train)
  pplme <- predict(lme_ss, test)
  rmse_lmer <- sqrt(mean((pplme - test$y)^2)) # 6.64
  return(list(pred = pplme, rmse = rmse_lmer))
}

data_model3 <- data_model2 |> 
  mutate(lme = map2(train, test, predict_lme), 
         pred_lme = map(lme, "pred"), 
         rmse_lme = map_dbl(lme, "rmse"))

write_rds(data_model3, file = "results/simple_example.rds")
summary(data_model3$rmse_hm)
summary(data_model3$rmse_lme)
#----------------------------------------------------------------------
#----------------------------------------------------------------------
# Plots 


summary_preds <- preds |> 
  dplyr::select(id, pred_hebart, pred_bart, pred_lme) |> 
  tidyr::unnest(c(pred_hebart, pred_bart, pred_lme), 
                names_repair = "unique")

rmses <- summary_preds |> 
  dplyr::group_by(id) |> 
  dplyr::summarise(
    rmse_hebart = rmse(pred_hebart, real_hebart),
    rmse_bart = rmse(pred_bart, real_bart),
    rmse_lme = rmse(pred_lme, real_lme)
  ) |> 
  tidyr::pivot_longer(cols = c(rmse_hebart, rmse_bart, rmse_lme)) 



data_model3 |>
  select(rmse_hm, rmse_lme) |> 
  tidyr::pivot_longer(cols = c(rmse_hm, rmse_lme)) |> 
  mutate(name = str_remove(name, "rmse\\_"), 
         name = str_replace(name, "hm", "CRF-HEBART"),
         name = str_to_upper(name)) |> 
  ggplot(aes(y = value, x = name)) +
  geom_boxplot(fill = "#F96209", alpha = 0.7) +
  #facet_wrap(~type, scales = 'free') +
  labs(y = "Estimated test RMSE", 
       x = 'Algorithms'
       #title = paste0("RMSE\nHE-BART: ", rss_hbart, ", LME: ", rss_lme)
  ) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  theme_fira()
ggsave(file = "paper/boxplot_simulated.png",
       width = 4, height = 3)
#----------------------------------------------------------------------
#----------------------------------------------------------------------
