library(tidyverse)
library(tidymodels)
library(lme4)
library(firatheme)
library(mhebart)

url <- "https://raw.githubusercontent.com/avehtari/ROS-Examples/master/Earnings/data/earnings.csv"
df <- read_csv(url) |> 
  filter(earn > 0) |> 
  na.omit()
df$earn <- log(df$earn)
df$age_cat <-  cut(df$age, 
                   breaks=c(0, 25, 35, 45, 60, 100), 
                   labels=c("age_25","age_35","age_45", 
                            "age_60", "age_100"))

df |> 
  ggplot(aes(y = earn, x  = education)) +
  geom_point() +
  facet_wrap(~ethnicity+age_cat, scales = 'free') +
  theme_fira()

data <- vfold_cv(df, 10) |> 
  mutate(
    train = map(splits, training),
    test = map(splits, testing))

train <- data$train[[1]]
test <- data$test[[1]]

run_mhebart <- function(train, test){
  num_trees   <- 8
  
  # Model parameters -----------------------------------
  group_variables <-  c("ethnicity", "age_cat")
  formula        <- earn ~ height + weight + male + education
  
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
                      MCMC = list(iter = 100, 
                                  burn = 10, 
                                  thin = 1,
                                  sigma_phi_sd = 2)
  )
  
  pp <- predict_mhebart(newX = test, group_variables, 
                        hebart_posterior = hb_model, type = "mean")
  
  rmse_mhebart <-  sqrt(mean((pp - test$earn)^2))
  
  test$preds <- pp
  
  
  lme_ss <- lme4::lmer(earn ~ height + weight + male + education +
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

saveRDS(runs, "results/model_earnings.rds")
#--------------------------------------------------------
runs <- read_rds("results/model_earnings.rds")
runs2 <-  runs |> 
  mutate(
    rmse_lme = map_dbl(all,"rmse_lme"),
    rmse_mhebart = map_dbl(all,"rmse_mhebart"),
    pred_lme = map(all, "pred_lme"),
    pred_hebart = map(all, "pred_hebart"),
    id = 1:n(), 
    test = map(all, "test")
  )

runs2$rmse_lme |> summary()
runs2$rmse_mhebart |> summary()
#--------------------------------------------------------
predict_other <- function(test, pred, pred_lme){
  others <- which(test$ethnicity == 'White')
  y_other <- test$earn[others]
  pred_other <- pred[others]
  pred_lme_other <- pred_lme[others]
  
  rmse <- sqrt(mean((pred_other - y_other)^2))
  rmse_lme <- sqrt(mean((pred_lme_other - y_other)^2))
  
  return(list(rmse = rmse, rmse_lme = rmse_lme))
}

runs4 <- runs2 |> 
  mutate(rmses_other = 
           pmap(list(test, pred_hebart, pred_lme), predict_other))
runs4$rmses_other |> map_dbl("rmse") |> mean()
runs4$rmses_other |> map_dbl("rmse_lme") |> mean()
#--------------------------------------------------------
predict_max <- function(test, pred, pred_lme){
  y <- test$earn
  max_hb <- which.max((pred - y)^2)
  max_lme <- which.max((pred_lme - y)^2)
  
  return(list(max_hb = test |> slice(max_hb), 
              max_lme = test |> slice(max_lme)))
}

runs4 <- runs2 |> 
  mutate(maxes = 
           pmap(list(test, pred_hebart, pred_lme), predict_max))

runs4$maxes |> map_df("max_hb")  |> 
  View()
runs4$maxes |> map_df("max_lme") |> View()

remove <- function(test, pred, pred_lme){
  y <- test$earn
  max_hb <- which.max((pred - y)^2)
  max_lme <- which.max((pred_lme - y)^2)
  
  y_new <- y[-max_hb]
  rmse_mhebart <-  sqrt(mean((pred[-max_hb] - y_new)^2))
  
  y_new <- y[-max_lme]
  rmse_lme <-  sqrt(mean((pred_lme[-max_lme] - y_new)^2))
  
  
  return(list(rmse_mhebart = rmse_mhebart, 
              rmse_lme = rmse_lme))
}

runs4 <- runs2 |> 
  mutate(new_rmse = 
           pmap(list(test, pred_hebart, pred_lme), remove))

runs4$new_rmse |> map_dbl("rmse_mhebart") |> mean()

runs4$new_rmse |> map_dbl("rmse_lme") |> mean()

BottleRocket2 = c("#FAD510", "#CB2314", "#0E86D4",
                           "#1E1E1E", "#18A558")

#-------------------------------------------------------------
# RMSE per groups



#-------------------------------------------------------------
# runs <- models_gapminder
runs3 <-  runs |> 
  mutate(
    rmse_lme = map_dbl(all,"rmse_lme"),
    rmse_mhebart = map_dbl(all,"rmse_mhebart"),
    pred_lme = map(all, "pred_lme"),
    pred_hebart = map(all, "pred_hebart"),
    id = 1:n(), 
    test = map(all, "test")
  )  |> 
  unnest(pred_lme, pred_hebart, test) |> 
  select(earn, pred_hebart, pred_lme, rmse_lme, rmse_mhebart,
         height, weight, male, education, age_cat, ethnicity) |> 
  mutate(
    error = (pred_hebart - earn)^2, 
    error_lme = (pred_lme - earn)^2, 
  ) |> 
  ungroup() |> 
  group_by(age_cat, ethnicity) |> 
  summarise(
    rmse = sqrt(mean(error)),
    rmse_lme = sqrt(mean(error_lme)))

runs3 |> 
  mutate(id = rmse < rmse_lme) |> 
  arrange(id)

#-------------------------------------------------------------
# summarise(
#     y = mean(earn), 
#     pred_lme_mean = mean(pred_lme), 
#     pred_mhebart_mean = mean(pred_hebart), 
#     low_lme = pred_lme_mean - 1.95 * sd_lme,
#     upp_lme = pred_lme_mean + 1.95 * sd_lme,
#     low_mhebart = pred_mhebart_mean - 1.95 * sd_mhebart,
#     upp_mhebart = pred_mhebart_mean + 1.95 * sd_mhebart
#   )
# 
# 
# 
# runs3 |>  
#   filter(ethnicity %in% c("White", "Hispanic")) |> 
#   ggplot(aes(x = education, y = pred_mhebart_mean)) +
#   facet_wrap(~ethnicity+age_cat, ncol = 2, scales = 'free_y') +
#   geom_ribbon(aes(ymin=low_mhebart, ymax=upp_mhebart),
#               fill = BottleRocket2[3], alpha = 0.5) + 
#   geom_ribbon(aes(ymin=low_lme, ymax=upp_lme),
#               fill = BottleRocket2[2], alpha = 0.3) + 
#   geom_line(aes(colour = BottleRocket2[2]), size = 0.7) +
#   geom_line(aes(y = pred_lme_mean), 
#             colour = BottleRocket2[2], size = 0.7) +
#   geom_line(aes(colour = BottleRocket2[3]), size = 0.7) +
#   geom_line(colour = BottleRocket2[3], size = 0.7) +
#   geom_point(aes(x = education, y = y, colour =  'black'), size = 0.25) + 
#   geom_point(aes(x = education, y = y), 
#              colour =  'black', size = 0.25) + 
#   scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) +
#   scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
#   labs(y = "Life expectancy (in years)", 
#        x = 'Passage of time', 
#        title = "Average predictions per group for MHEBART and LME") + 
#   theme_linedraw(12) +
#   scale_colour_manual(
#     name="Source:",
#     values = c("black", BottleRocket2[3], BottleRocket2[2]), 
#     labels = c("Data", "MHEBART", "LME"), 
#     
#     guide = guide_legend(override.aes = list(
#       size = c(2, 2, 2)))) +
#   theme(panel.spacing.x = unit(0.5, "lines"), 
#         legend.position = "bottom")
# 
# ggsave(file = "images/predictions_plot_earnings_mhebart.png",
#        width = 8, height = 8)


