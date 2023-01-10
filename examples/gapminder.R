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
library(mhebart)
load("data/gapminder_recent_g20.RData")

# Dataset split  ------------------------------------
set.seed(2022)

df_real     <- gapminder_recent_g20 %>% 
  select(year, country, lifeExp, year0, decade0, continent, gdpPercap) |> 
  set_names(c('X1', 'country', 'y', "X2", "X3", "continent", "X4"))

df_real$X4 <- log(df_real$X4)
conts <- c("Americas", "Europe", "Africa")
df_real <- df_real |> filter(continent %in% conts)
years       <- unique(df_real$X1)
to_remove   <- sample(years, 10)
train       <- df_real |> filter(!(X1 %in% to_remove))
test        <- df_real |> filter(X1 %in% to_remove)
num_trees   <- 20

# Model parameters -----------------------------------
group_variables <-  c("country", "continent")
formula        <- y ~ X1 + X4 + X3

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
                   MCMC = list(iter = 1250, 
                               burn = 50, 
                               thin = 1,
                               sigma_phi_sd = 2)
                   )
hb_model
pp <- predict_mhebart(newX = test, group_variables, 
                      hebart_posterior = hb_model, type = "mean")
rmse_mhebart <-  sqrt(mean((pp - test$y)^2))

test$preds <- pp

rmse_mhebart <-  sqrt(mean((pp - test$y)^2))

sqrt(mean((pp - test$y)^2)) # 1.907164
cor(pp, scale(test$y))  # 0.9881025

lme_ss <- lme4::lmer(y ~ X1 + X4 + X3 + (1|country) + (1|continent), train)
pplme <- predict(lme_ss, test)
rmse_lmer <- sqrt(mean((pplme - test$y)^2)) # 3.991818
cor(pplme, test$y) # 0.936175

test$pplme <- pplme
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




