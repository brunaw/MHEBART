library(foreign)
library(tidyverse)
library(mhebart)
library(lme4)
heights <- read.dta ("examples/heights.dta")
attach(heights)

# define variables
age <- 90 - yearbn                     # survey was conducted in 1990
age[age<18] <- NA
age.category <- ifelse (age<35, 1, ifelse (age<50, 2, 3))
eth <- ifelse (race==2, 1, ifelse (hisp==1, 2, ifelse (race==1, 3, 4)))
male <- 2 - sex

# (for simplicity) remove cases with missing data
ok <- !is.na (earn+height+sex+age.category+eth) & earn>0 & yearbn>25
heights.clean <- as.data.frame (cbind (earn, height, sex, race, hisp, ed, age,
                                       age.category, eth, male)[ok,])
n <- nrow (heights.clean)
height.jitter.add <- runif (n, -.2, .2)

# rename variables
y <- log(earn[ok])
x <- height[ok]
n <- length(y[ok])
n.age <- 3
n.eth <- 4
age <- age.category
age.ok <- age[ok]
eth.ok <- eth[ok]

## Regression centering the predictors
##M1 <- lmer (y ~ x.centered + (1 + x.centered | eth) + (1 + x.centered | age) + (1 + x.centered | eth:age))
x.centered <- x - mean(x)
df <- data.frame(
  x = x.centered, y, eth = eth.ok, age = age.ok
)
group_variables <-  c("eth", "age")
formula        <- y ~ x


data = df

#dataList.2 <- list(N=length(y),y=y,x=x.centered,n_eth=n.eth,n_age=n.age,eth=eth.ok,age=age[ok])

M1 <- lmer (y ~ x.centered + (1 + x.centered | eth.ok) + (1 + x.centered | age.ok))
pp <- predict(M1)
mse <- mean((pp - y)^2)
rmse_lmer <- sqrt(mse)


hb_model <- mhebart(formula = y ~ x,
                   data = df,
                   group_variables = c("eth", "age"), 
                   num_trees = 10,
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
                   MCMC = list(iter = 50, 
                               burn = 25, 
                               thin = 1,
                               sigma_phi_sd = 0.5)
)
pp <- predict_mhebart(newX = df, c("eth", "age"), 
                     hebart_posterior = hb_model, type = "mean")
rmse_mhebart <-  sqrt(mean((pp - df$y)^2))
df$preds <- pp
df |> 
  ggplot(aes(x.centered, y)) +
  geom_point() +
  geom_line(aes(y = preds), colour = 'red') + 
  geom_point(aes(y = preds), colour = 'red') + 
  facet_wrap(~age+eth) +
  ggtitle(paste0("MHEBART RMSE:", round(rmse_mhebart, 2), 
                 ",\nLMER RMSE:", round(rmse_lmer, 2)))
  
