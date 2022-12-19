library(foreign)
library(hebartBase)
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