library(tidyverse)
library(tidymodels)
library(lme4)
library(brms)
library(firatheme)
library(ggokabeito)
library(ggrepel)
devtools::load_all(".")

load("data/pupils.Rdata")
pupils

pupils_crossed = lmer(
  achievement ~ 
    sex 
  + ses 
  + (1|primary_school_id) 
  + (1|secondary_school_id),
  data = pupils
)

pplme <- predict(pupils_crossed, pupils)
rmse_lmer <- sqrt(mean((pplme - pupils$achievement)^2)) # 6.64
rmse_lmer

num_trees <- 6

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
  MCMC = list(iter = 50, 
              burn = 10, 
              thin = 1,
              sigma_phi_sd = 1)
)

pp <- predict_mhebart(newX = pupils, 
                      group_variables = c("primary_school_id", "secondary_school_id"),
                      hebart_posterior  = hb_model, 
                      type = "mean")
rmse <- sqrt(mean((pp - pupils$achievement)^2))


#-----------------------------------------------------------------------
# bayes_seed <- 1000
# model_country_continent_year <- brm(
#   bf(lifeExp ~ year + (1 + year | continent / country)),
#   data = gapminder,
#   chains = 2, seed = bayes_seed,
#   iter = 1000  # Double the number of iterations to help with convergence
# )

model_country_continent_year %>%
  emmeans(~ year + continent:country,
          at = list(year = c(0), country = countries$country),
          nesting = "country %in% continent",
          epred = TRUE, re_formula = NULL, allow_new_levels = TRUE)


gapminder <- gapminder |> 
  dplyr::mutate_if(is.factor, ~gsub('[[:punct:] ]+',' ',.x)) |> 
  dplyr::mutate_if(is.character, stringr::str_squish) |> 
  dplyr::mutate_if(is.character, stringr::str_squish) |> 
  dplyr::mutate_if(is.character, stringr::str_to_lower) |>
  dplyr::mutate_if(is.character, abjutils::rm_accent) |> 
  dplyr::mutate_if(is.character, ~stringr::str_replace_all(.x, " ", "_")) |> 
  dplyr::mutate_if(is.factor, ~stringr::str_remove_all(.x, "\\,|\\."))
table(gapminder$country)

data <- vfold_cv(gapminder, 10) |> 
  mutate(
    train = map(splits, training),
    test = map(splits, testing))
write_rds(data, paste0("results/data_gapminder.rds")) 

fit_lme <- function(train, test){
  # gapminder <- gapminder |> sample_n(1000)
  lme_ss <- lme4::lmer(
    lifeExp ~ year + (1 + year | continent / country), 
    train)
  pplme <- predict(lme_ss, test)
  rmse_lmer <- sqrt(mean((pplme - test$lifeExp)^2)) # 6.64
  
  return(list(pred = pplme, rmse = rmse_lmer))
}


data_model2 <- data |> 
  mutate(lme = map2(train, test, fit_lme), 
         pred_lme = map(lme, "pred"), 
         rmse_lme = map_dbl(lme, "rmse"))
data_model2$rmse_lme

train = data$train[[1]]
test = data$test[[1]]
fit_mhebart <- function(train, test, id){
  num_trees <- 20
  
  hb_model <- mhebart(
    formula = lifeExp ~ year,
    data = train,
    group_variables = c("continent", "country"), 
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
    MCMC = list(iter = 500, 
                burn = 100, 
                thin = 1,
                sigma_phi_sd = 1)
  )
  #write_rds(hb_model, paste0("results/hb_model_", i, ".rds")) 
  
  pp <- predict_mhebart(newX = test, 
                        group_variables = c("continent", "country"),
                        hebart_posterior  = hb_model, 
                        type = "mean")
  rmse <- sqrt(mean((pp - test$lifeExp)^2))
  
  return(list(pred = pp, rmse = rmse))
}