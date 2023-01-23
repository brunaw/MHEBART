library(tidyverse)
library(firatheme)
devtools::load_all(".")


url <- "http://www.stat.columbia.edu/~gelman/arm/examples/electric.company/electric.dat"

electric <- read.table (url, header=TRUE) 

electric2 <- data_frame(
  y = c(electric$treated.Posttest, electric$control.Posttest), 
  sup = c(electric$Supplement., electric$Supplement.), 
  treatment = rep(1:2, each = nrow(electric)), 
  group_1 = c(electric$Grade, electric$Grade), 
  group_2 = c(electric$City, electric$City)
)

electric2 |> 
  ggplot(aes(y = y)) +
  geom_boxplot(aes(group = treatment)) +
  facet_wrap(~group_1+group_2)

rows <- sample(1:nrow(electric2), size = 45)
test <- electric2 |> slice(rows)
train <- electric2 |> slice((1:nrow(electric2))[-rows])

lme_ss <- lme4::lmer(y ~ treatment + (1|group_1) + (group_2), train)
pplme <- predict(lme_ss, test)
rmse_lmer <- sqrt(mean((pplme - test$y)^2)) # 6.64
rmse_lmer

num_trees <- 20

hb_model <- mhebart(
  formula = y ~ treatment,
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
  MCMC = list(iter = 2500, 
              burn = 500, 
              thin = 1,
              sigma_phi_sd = 1)
)
hb_model

pp <- predict_mhebart(newX = test, 
                      group_variables = c("group_1", "group_2"),
                      hebart_posterior  = hb_model, 
                      type = "mean")
rmse <- sqrt(mean((pp - test$y)^2))
rmse
