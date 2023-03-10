#-------------------------------------------------------
# A model fit for the gapminder data using the year as a 
# categorical feature ----------------------------------
#-------------------------------------------------------
# Exemplifying:
# Package loading  ----------------------------------
library(magrittr)
library(ggplot2)
library(tidymodels)
library(firatheme)
library(mhebart)
load("data/gapminder_recent_g20.RData")

# Dataset split  ------------------------------------
set.seed(2022)

# Account for random things that happened in each year
df_real     <- gapminder_recent_g20 %>% 
  select(year, country, lifeExp, year0, decade0, continent, gdpPercap) |> 
  set_names(c('X1', 'country', 'y', "X2", "X3", "continent", "X4"))
View(df_real)
df_real$X4 <- log(df_real$X4)
conts <- c("Americas", "Europe", "Africa")
df_real <- df_real |> filter(continent %in% conts)
years       <- unique(df_real$X1)

df_real$year_cat <-  cut(df_real$X1, 
                         breaks=seq(1949, 2020, by = 10),
                         labels=paste("year", seq(1950, 2018, by = 10)))
df_real$X1 |> table()
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------

run_gapminder <- function(x){
  to_remove   <- sample(years, 15)
  train       <- df_real |> filter(!(X1 %in% to_remove))
  test        <- df_real |> filter(X1 %in% to_remove)
  num_trees   <- 10
  
  # Model parameters -----------------------------------
  group_variables <-  c("country", "year_cat")
  formula        <- y ~ X2 + X4 + X3
  
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
                      MCMC = list(iter = 500, 
                                  burn = 100, 
                                  thin = 1,
                                  sigma_phi_sd = 2)
  )

  pp <- predict_mhebart(newX = test, group_variables, 
                        hebart_posterior = hb_model, type = "mean")
  
  rmse_mhebart <-  sqrt(mean((pp - test$y)^2))
  test$preds <- pp
  
  lme_ss <- lme4::lmer(y ~ X2 + X4 + X3 + (1|country) + (1|year_cat), train)
  pplme <- predict(lme_ss, test)
  rmse_lmer <- sqrt(mean((pplme - test$y)^2)) # 3.991818
  
  return(list(train = train, test = test, 
              hb_model = hb_model, rmse_mhebart = rmse_mhebart, 
              pred_hebart = pp, 
              lme_model = lme_ss, rmse_lme = rmse_lmer, 
              pred_lme = pplme
  )) 
}

runs <- tibble(
  all = map(1:10, run_gapminder)
)

saveRDS(runs, "models_gapminder.rds")
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# runs <- models_gapminder

new_lme <- function(train, test){
  
  lme_ss <- me_ss <- lme4::lmer(y ~ X2 + X4 + X3 + (X2 + X4 + X3|country) + (X2 + X4 + X3|year_cat), train)
  pplme <- predict(lme_ss, test)
  rmse_lmer <- sqrt(mean((pplme - test$y)^2)) # 3.991818
  
  return(list(lme_model_new = lme_ss, 
              rmse_lme_new = rmse_lmer, 
              pred_lme_new = pplme
  )) 
}


runs2 <- runs |> 
  mutate(train = map(all, "train"), 
         test = map(all, "test"), 
         new_lme = map2(train, test, new_lme)
  ) 

runs2 <-  runs2 |> 
  mutate(
    rmse_lme = map_dbl(new_lme,"rmse_lme_new"),
    rmse_mhebart = map_dbl(all,"rmse_mhebart"),
    pred_lme = map(new_lme, "pred_lme_new"),
    pred_hebart = map(all, "pred_hebart"),
    id = 1:n(), 
    test = map(all, "test")
  )
runs3$rmse_lme
runs3$rmse_mhebart

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
runs3 <-  runs2  |> 
  unnest(pred_lme, pred_hebart, test) |> 
  dplyr::select(y, pred_hebart, pred_lme, rmse_lme, rmse_mhebart, country, X1, year_cat) |> 
  mutate(sd_lme = sd(pred_lme)/sqrt(n()), 
         sd_mhebart = sd(pred_hebart)/sqrt(n())) |> 
  group_by(X1, country, year_cat) |> 
  summarise(
    y = mean(y), 
    pred_lme_mean = mean(pred_lme), 
    pred_mhebart_mean = mean(pred_hebart), 
    low_lme = pred_lme_mean - 1.95 * sd_lme,
    upp_lme = pred_lme_mean + 1.95 * sd_lme,
    low_mhebart = pred_mhebart_mean - 1.95 * sd_mhebart,
    upp_mhebart = pred_mhebart_mean + 1.95 * sd_mhebart
  )

# -------------------------------------------------------------
# Plots -------------------------------------------------------
avg_rmse_lme <- mean(runs2$rmse_lme)
upp_lme <- avg_rmse_lme + 1.96 * sd(runs2$rmse_lme)/sqrt(10)
low_lme <- avg_rmse_lme - 1.96 * sd(runs2$rmse_lme)/sqrt(10)


avg_rmse <- mean(runs2$rmse_mhebart)
upp <- avg_rmse + 1.96 * sd(runs2$rmse_mhebart)/sqrt(10)
low <- avg_rmse - 1.96 * sd(runs2$rmse_mhebart)/sqrt(10)


BottleRocket2 = c("#FAD510", "#CB2314", "#0E86D4",
                           "#1E1E1E", "#18A558")
                           
selected_countries <- c(
  "South Africa", "Russia", "China", 
  "Turkey", "Indonesia", "Brazil",
  "Argentina", "France", "United States"
)


runs3 |>  
  filter(country %in% selected_countries) |> 
  ggplot(aes(x = X1, y = pred_mhebart_mean)) +
  facet_wrap(~country, ncol = 2, scales = 'free_y') +
  geom_ribbon(aes(ymin=low_mhebart, ymax=upp_mhebart),
              fill = BottleRocket2[3], alpha = 0.5) + 
  geom_ribbon(aes(ymin=low_lme, ymax=upp_lme),
              fill = BottleRocket2[2], alpha = 0.3) + 
  geom_line(aes(colour = BottleRocket2[2]), size = 0.7) +
  geom_line(aes(y = pred_lme_mean), 
            colour = BottleRocket2[2], size = 0.7) +
  geom_line(aes(colour = BottleRocket2[3]), size = 0.7) +
  geom_line(colour = BottleRocket2[3], size = 0.7) +
  geom_point(aes(x = X1, y = y, colour =  'black'), size = 0.25) + 
  geom_point(aes(x = X1, y = y), 
             colour =  'black', size = 0.25) + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  labs(y = "Life expectancy (in years)", 
       x = 'Passage of time', 
       title = 
         
         paste0("Average predictions per group for Crossed-RE HEBART and LME \n", 
                "Crossed-RE HEBART Test RMSE: ", 
               paste0(round(avg_rmse, 2), " [", round(low, 2), ",", round(upp, 2), "]"), "\n", 
               "LME Test RMSE: ", 
               paste0(round(avg_rmse_lme, 2), " [", round(low_lme, 2), ",", round(upp_lme, 2), "]"), "\n"
                )) + 
  theme_linedraw(12) +
  scale_colour_manual(
    name="Source:",
    values = c("black", BottleRocket2[3], BottleRocket2[2]), 
    labels = c("Data", "Crossed-RE HEBART", "LME"), 
    
    guide = guide_legend(override.aes = list(
      size = c(2, 2, 2)))) +
  theme(panel.spacing.x = unit(0.5, "lines"), 
        legend.position = "bottom")

ggsave(file = "images/predictions_plot_gapminder_mhebart.png",
       width = 8, height = 8)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------




