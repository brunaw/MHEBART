# In case you're installing, building, or removing the package:
# remove.packages("hebartBase")
# devtools::document()
# devtools::check()
# devtools::install()

# Exemplifying:
# Package loading  ----------------------------------
library(magrittr)
library(ggplot2)
library(tidymodels)
library(firatheme)
#library(mhebart)
devtools::load_all(".")
load("data/gapminder_recent_g20.RData")

# Dataset split  ------------------------------------
set.seed(2022)

# Account for random things that happened in each year
df_real     <- gapminder_recent_g20 %>% 
  select(year, country, lifeExp, year0, decade0, continent, gdpPercap) |> 
  set_names(c('X1', 'country', 'y', "X2", "X3", "continent", "X4"))

df_real$X4 <- log(df_real$X4)
conts <- c("Americas", "Europe", "Africa")
df_real <- df_real |> filter(continent %in% conts)
years       <- unique(df_real$X1)

df_real$year_cat <-  cut(df_real$X1, 
                   breaks=seq(1949, 2020, by = 10),
                   labels=paste("year", seq(1950, 2018, by = 10)))
df_real$X1 |> table()
#df_real$year_cat <- as.factor(df_real$X1)


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
  
  saveRDS(hb_model, paste0("models_gapminder", x, ".rds")) 
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

runs2 <-  runs |> 
  mutate(
    rmse_lme = map_dbl(all,"rmse_lme"),
    rmse_mhebart = map_dbl(all,"rmse_mhebart"),
    pred_lme = map(all, "pred_lme"),
    pred_hebart = map(all, "pred_hebart"),
    id = 1:n(), 
    test = map(all, "test")
  )  |> 
  unnest(pred_lme, pred_hebart, test) |> 
  select(y, pred_hebart, pred_lme, rmse_lme, rmse_mhebart, country, X1, year_cat) |> 
  mutate(sd_lme = sd(pred_lme)/sqrt(n()), 
         sd_mhebart = sd(pred_hebart)/sqrt(n()))|> 
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

BottleRocket2 = c("#FAD510", "#CB2314", "#0E86D4",
                           "#1E1E1E", "#18A558")
                           
selected_countries <- c(
  "South Africa", "Russia", "China", 
  "Turkey", "Indonesia", "Brazil",
  "Argentina", "France", "United States"
)

runs2 |> View()
runs2 |>  
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
       title = "Average predictions per group for MHEBART and LME") + 
  theme_linedraw(12) +
  scale_colour_manual(
    name="Source:",
    values = c("black", BottleRocket2[3], BottleRocket2[2]), 
    labels = c("Data", "MHEBART", "LME"), 
    
    guide = guide_legend(override.aes = list(
      size = c(2, 2, 2)))) +
  theme(panel.spacing.x = unit(0.5, "lines"), 
        legend.position = "bottom")
    
ggsave(file = "images/predictions_plot_gapminder_mhebart.png",
       width = 8, height = 8)


  test |> 
    ggplot(aes(X1, y)) +
    geom_point() +
    geom_line(aes(y = preds), colour = 'red') + 
    geom_point(aes(y = preds), colour = 'red') +
  geom_line(aes(y = pplme), colour = 'blue') + 
  geom_point(aes(y = pplme), colour = 'blue') + 
  facet_wrap(~country+continent) +
  ggtitle(paste0("MHEBART RMSE:", round(rmse_mhebart, 2), 
                 ",\nLMER RMSE:", round(rmse_lmer, 2)), 
          "\n red dots: MHEBART prediction, blue dots: LMER prediction")

# Comparison to BART --------------------------
# bart_0 = dbarts::bart2(y ~ X1 + X4 + X3, 
#                        data = train,
#                        test = test,
#                        keepTrees = TRUE)
# ppbart <- bart_0$yhat.test.mean
# sqrt(mean((ppbart - test$y)^2)) # 7.944524
# cor(ppbart, test$y) #    0.698455
# 
# ppbart <- bart_0$yhat.train.mean
# sqrt(mean((ppbart - train$y)^2)) # 0.8950348- 100 trees

# BART+Group
# bart_0 = dbarts::bart2(y ~ X1 + X2 + X3 + group, 
#                        data = train,
#                        test = test,
#                        keepTrees = TRUE)
# ppbart <- bart_0$yhat.test.mean
# sqrt(mean((ppbart - test$y)^2)) # 0.9425852
# cor(ppbart, test$y) #    0.99683
# 
# ppbart <- bart_0$yhat.train.mean
# sqrt(mean((ppbart - train$y)^2)) # 0.3275252

# Comparison to LME --------------------------
# Average predictions 
preds_y <- data.frame(test, pred = pp, pred_lme = pplme)

preds_y |> 
  filter(group %in% c("China", "South Africa", "Russia", 
                      "United States", "Mexico", 
                      "Canada")) |> 
  ggplot(aes(x = X1, y = y)) +
  geom_point(colour = "gray") +
  geom_line(aes(y = pred, colour= "#75E6DA"), size = 1.3) +
  geom_line(aes(y = pred), colour= "#75E6DA", size = 1,
            alpha = 0.7) +
  geom_line(aes(y = pred_lme, colour= "#F96209"),  size = 1.3) +
  geom_line(aes(y = pred_lme), colour= "#F96209",  size = 1, 
            alpha = 0.7) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) +
  labs(x = "Life expectancy", 
       y = 'Predictions', 
       #title = paste0("RMSE\nHE-BART: ", rss_hbart, ", LME: ", rss_lme)
  ) + 
  facet_wrap(~group, scales = "free", ncol = 2) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  theme_linedraw(14) +
  scale_colour_manual(
    name="Source:",
    values=c(Data="gray", 
             `HEBART Prediction`="#75E6DA", 
             `LME Prediction`= '#F96209'), 
    guide = guide_legend(override.aes = list(
      size = c(3, 3, 3), shape = c(16, 16, 16)))) + 
  theme_fira()
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------




