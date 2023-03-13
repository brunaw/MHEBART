library(tidyverse)
library(tidymodels)
library(lme4)
library(firatheme)
library(mhebart)

url <- "https://raw.githubusercontent.com/avehtari/ROS-Examples/master/Earnings/data/earnings.csv"
df <- read_csv(url)
df$earn <- log(df$earn+0.01)
df$age_cat <-  cut(df$age, 
                   breaks=c(0, 25, 35, 45, 60, 100), 
                   labels=c("age_25","age_35","age_45", 
                            "age_60", "age_100"))

df <- na.omit(df)
#df <- df  |> sample_n(500)
data <- vfold_cv(df, 10) |> 
  mutate(
    train = map(splits, training),
    test = map(splits, testing))
train <- data$train[[1]]
test <- data$test[[1]]

run_mhebart <- function(train, test){
  num_trees   <- 12
  
  # Model parameters -----------------------------------
  group_variables <-  c("ethnicity", "age_cat")
  formula        <- earn ~ height + weight + male+ education + age
  
  # Running the model ----------------------------------
  hb_model <- mhebart(formula,
                      data           = train,
                      group_variables = group_variables, 
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
                                   sigma_phi = 1),
                      MCMC = list(iter = 600, 
                                  burn = 100, 
                                  thin = 1,
                                  sigma_phi_sd = 2)
  )
  
  pp <- predict_mhebart(newX = test, group_variables, 
                        hebart_posterior = hb_model, type = "mean")
  
  rmse_mhebart <-  sqrt(mean((pp - test$earn)^2))
  
  test$preds <- pp
  
  
  lme_ss <- lme4::lmer(earn ~ height + weight +
                         (1|ethnicity) + (1|age_cat), train)
  pplme <- predict(lme_ss, test)
  rmse_lmer <- sqrt(mean((pplme - test$earn)^2)) 
  
  return(list(train = train, test = test, 
              hb_model = hb_model, rmse_mhebart = rmse_mhebart, 
              pred_hebart = pp, 
              lme_model = lme_ss, rmse_lme = rmse_lmer, 
              pred_lme = pplme
  )) 
}

runs <- data |> 
  mutate(
  all = map2(train, test, run_mhebart)
)
saveRDS(runs, "model_earnings.rds")
runs2 <-  runs |> 
  mutate(
    rmse_lme = map_dbl(all,"rmse_lme"),
    rmse_mhebart = map_dbl(all,"rmse_mhebart"),
    pred_lme = map(all, "pred_lme"),
    pred_hebart = map(all, "pred_hebart"),
    id = 1:n(), 
    test = map(all, "test")
  )

runs2$rmse_lme
runs2$rmse_mhebart
