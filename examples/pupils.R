library(tidyverse)
library(tidymodels)
library(lme4)
library(brms)
library(firatheme)
library(ggokabeito)
library(ggrepel)
devtools::load_all(".")

load("data/pupils.Rdata")
# pupils

# y <- pupils$achievement
# y_min <- min(y)
# y_max <- max(y)
# y_scale <- (y - y_min)/(y_max - y_min) - 0.5
# n <- length(y_scale)
# pupils$achievement = y_scale
# mean(y_scale)


pupils_crossed = lmer(
  achievement ~ 1
  + (1|primary_school_id) 
  + (1|secondary_school_id),
  data = pupils
)
#View(pupils)
pplme <- predict(pupils_crossed, pupils)
rmse_lmer <- mean((pplme - pupils$achievement)^2) # 6.64
rmse_lmer

# pupils |> 
#   mutate(pred = pp) |> 
#   group_by(primary_school_id, secondary_school_id) |> 
#   summarise(y = mean(achievement), 
#             pred = mean(pred),
#             n = n()) |> 
#   View()

num_trees <- 4

hb_model <- mhebart(
  formula = achievement ~ sex + ses,
  data = pupils,
  group_variables = c("primary_school_id", "secondary_school_id"), 
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
  MCMC = list(iter = 100, 
              burn = 0, 
              thin = 1,
              sigma_phi_sd = 1), 
  stumps = TRUE
)
#hb_model$trees[[2]]
#hb_model$mse

pp <- predict_mhebart(newX = pupils, 
                      group_variables = c("primary_school_id", "secondary_school_id"),
                      hebart_posterior  = hb_model, 
                      type = "mean")

rmse <- mean((pp - pupils$achievement)^2)

rmse <- mean((pp - y_scale)^2)
rmse
y <- pupils$achievement
y_min <- min(y)
y_max <- max(y)
y_scale <- (y - y_min)/(y_max - y_min) - 0.5

# mse -> 0.014
# curr_trees = hb_model$trees[[20]]
# 
# predictions_all <- rep(0, length = length(y))
# 
# for(n_g in 1:2){
#   
#   groups <- data[, grouping_variables_names[n_g]]
#   if(!is.vector(groups)){
#     groups <- dplyr::pull(groups, !!grouping_variables_names[n_g])
#   }
#   
#   
#   preds <- get_group_predictions(
#     trees = curr_trees[[n_g]], 
#     X, 
#     groups, 
#     single_tree = num_trees == 1,
#     old_groups = groups
#   )
#   
#   predictions_all <- predictions_all + preds 
#   
# }
# mean((y_scale - predictions_all)^2)
# mean((pp - predictions_all)^2)



# pupils |>
#   mutate(pred = pp) |>
#   filter(secondary_school_id %in% (1:5), primary_school_id %in% (1:5)) |>
#   ggplot(aes(y = achievement, x = ses)) +
#   geom_point() +
#   geom_point(aes(y = pred), colour = 'red') +
#   facet_wrap(~secondary_school_id+primary_school_id)

#-----------------------------------------------------------------------
