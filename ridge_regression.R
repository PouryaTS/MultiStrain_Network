library(tidyverse)

df <- read_csv("chiara_collab/USCeoffDataforRegression_alpha_period.csv") |> 
  filter(! is.na(IM)) |> 
  mutate(EstimatedSelCoeff_znorm = scale(EstimatedSelCoeff),
         MD_znorm = scale(MD),
         SD_znorm = scale(SD),
         CV_znorm = scale(CV),
         IM_znorm = scale(IM))

# MD_log_znorm is achieved by doing scale(log(MD))
# I am creating variables without log in case we want them

# first let's look at collinearity, which we know should exist
model1  <- lm(EstimatedSelCoeff ~ MD + SD + IM, data = df)
summary(model1)
library(car)
vif(model1)
# VIFs for MD and SD really high, as we expected

# is this any different if we use scaled variables?
model2  <- lm(EstimatedSelCoeff_znorm ~ MD_znorm + SD_znorm + IM_znorm, data = df)
summary(model2)
vif(model2)
# same results
model3  <- lm(EstimatedSelCoeff_log_znorm ~ MD_log_znorm + SD_log_znorm + IM_log_znorm, data = df)
summary(model3)
vif(model3)
# VIF even worse with log znorm

# look at variables to check linearity
plot(df %>% select(EstimatedSelCoeff, MD, SD, IM))

# time to move to penalized regression
# first we have to create a matrix of covariates
# I'm thinking I will use the basic znormalized ones, not the log
library(glmnet)
# X: each row is an observation vector, a matrix but essentially a df
# -1 to make sure we don't double count the intercept
# not sure we actually need to znormalize the response
X         <- model.matrix(EstimatedSelCoeff ~ -1 + MD_znorm + SD_znorm + IM_znorm, 
                          data = df)
# repeat response here
Y         <- df$EstimatedSelCoeff_znorm
model_ridge  <- glmnet(x = X, y = Y, family = "gaussian", alpha = 0, 
                       intercept = TRUE)
# I think we will just use their default lambda, the lambda.1se

model_lasso  <- glmnet(x = X, y = Y, family = "gaussian", alpha = 1, 
                       intercept = TRUE)

model_elastic  <- glmnet(x = X, y = Y, family = "gaussian", alpha = 0.5, 
                         intercept = TRUE)
# It's useful to do these where we don't set lambda so we see how coefficient values change
# as lambda varies for these diagnostic plots
par(mfrow = c(1, 3))
plot(model_ridge, xvar = "lambda", label = TRUE)
plot(model_lasso, xvar = "lambda", label = TRUE)
plot(model_elastic, xvar = "lambda", label = TRUE)
# confirm that ridge keeps things in until it can't

model_ridge_cv   <- cv.glmnet(x = X, y = Y, family = "gaussian", alpha = 0,
                              type.measure = "mse", nfolds = nrow(df),
                              grouped = FALSE)

model_lasso_cv   <- cv.glmnet(x = X, y = Y, family = "gaussian", alpha = 1,
                              type.measure = "mse", nfolds = nrow(df),
                              grouped = FALSE)

model_elastic_cv_5   <- cv.glmnet(x = X, y = Y, family = "gaussian", alpha = 0.5,
                                  type.measure = "mse", nfolds = nrow(df),
                                  grouped = FALSE)

par(mfrow = c(1, 1))
plot(model_ridge_cv) # best model has all predictors, -log(lambda) around 5, meaning lambda around 0.01
plot(model_lasso_cv) # best model has 1 predictor, -log(lambda) around 6
plot(model_elastic_cv_5) # best model has 2 predictors, -log(lambda) around 6, I think MD and SD

print(model_ridge_cv) # 3 predictors
print(model_lasso_cv) # 1 predictor
print(model_elastic_cv_5) # 2 predictors

# these are the coefficients for the crossvalidated one
matrix(c(as.numeric(coef(model_ridge_cv)),
         as.numeric(coef(model_lasso_cv)), 
         as.numeric(coef(model_elastic_cv_5))),
       nrow = 3, byrow = TRUE, dimnames = list(Models = c('Ridge', 'Lasso', 'Elastic 0.5'),
                                               Variables = dimnames(coef(model_ridge_cv))[[1]])) %>% 
  t()

# next steps could be a cross validation on alpha (the degree to which it is ridge and lasso)

coef(model_ridge_cv, s = "lambda.min")
coef(model_ridge_cv, s = "lambda.1se")
