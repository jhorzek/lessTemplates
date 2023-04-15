## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(lessTemplates)
library(lavaan)
dataset <- lessTemplates::simulateRICLPM(seed = 123)

## ----echo = FALSE-------------------------------------------------------------
library(lessTemplates)
library(lavaan)
dataset <- simulateRICLPM(seed = 123)
head(dataset)

## -----------------------------------------------------------------------------
model <- "
# autoregressive and cross-lagged effects
eta1_(u) ~ eta1_(u-1) + eta2_(u-1)
eta2_(u) ~ eta1_(u-1) + eta2_(u-1)

# covariances between latent variables
eta1_(u) ~~ eta1_(u) + 0*eta2_(u)
eta2_(u) ~~ eta2_(u)

# measurements
eta1_(u) =~ 1*y1_(u)
eta2_(u) =~ 1*y2_(u)

# measurement error variances
y1_(u) ~~ 0*y1_(u)
y2_(u) ~~ 0*y2_(u)
"

## -----------------------------------------------------------------------------
clpm <- lessTemplates::CLPM(model = model, 
                            data = dataset)

## -----------------------------------------------------------------------------
fit <- lavaan::sem(model = clpm$model,
                   data = clpm$data)
coef(fit)

## -----------------------------------------------------------------------------
summary(fit)

## -----------------------------------------------------------------------------
model <- "
# autoregressive and cross-lagged effects
eta1_(u) ~ a11*eta1_(u-1) + a12*eta2_(u-1)
eta2_(u) ~ a21_(u)*eta1_(u-1) + a22*eta2_(u-1)

# covariances between latent variables
eta1_(u) ~~ v11*eta1_(u) + 0*eta2_(u)
eta2_(u) ~~ v22*eta2_(u)

# measurements
eta1_(u) =~ 1*y1_(u)
eta2_(u) =~ 1*y2_(u)

# measurement error variances
y1_(u) ~~ 0*y1_(u)
y2_(u) ~~ 0*y2_(u)
"

## -----------------------------------------------------------------------------
clpm <- lessTemplates::CLPM(model = model, 
                            data = dataset)
fit <- lavaan::sem(model = clpm$model,
                   data = clpm$data)
coef(fit)

## -----------------------------------------------------------------------------
summary(fit)

## -----------------------------------------------------------------------------
model <- "
# autoregressive and cross-lagged effects
eta1_(u) ~ a11*eta1_(u-1) + a12*eta2_(u-1)
eta2_(u) ~ a21*eta1_(u-1) + a22*eta2_(u-1)

# measurement occasion 4 differs from the other
# measurement occasions:
eta1_(4) ~ a11_4*eta1_(3) + a12_4*eta2_(3)

# covariances between latent variables
eta1_(u) ~~ v11*eta1_(u) + 0*eta2_(u)
eta2_(u) ~~ v22*eta2_(u)

# measurements
eta1_(u) =~ 1*y1_(u)
eta2_(u) =~ 1*y2_(u)

# measurement error variances
y1_(u) ~~ 0*y1_(u)
y2_(u) ~~ 0*y2_(u)
"

## -----------------------------------------------------------------------------
clpm <- lessTemplates::CLPM(model = model, 
                            data = dataset)
fit <- lavaan::sem(model = clpm$model,
                   data = clpm$data)
coef(fit)

## -----------------------------------------------------------------------------
summary(fit)

## -----------------------------------------------------------------------------
model <- "
# autoregressive and cross-lagged effects
eta1_(u) ~ a11*eta1_(u-1) + a12*eta2_(u-1)
eta2_(u) ~ a21*eta1_(u-1) + a22*eta2_(u-1)

# covariances between latent variables
eta1_(u) ~~ v11*eta1_(u) + 0*eta2_(u)
eta2_(u) ~~ v22*eta2_(u)

# measurements
eta1_(u) =~ 1*y1_(u)
eta2_(u) =~ 1*y2_(u)

# measurement error variances
y1_(u) ~~ 0*y1_(u)
y2_(u) ~~ 0*y2_(u)

## Random Intercepts

RI1 =~ 1*y1_(u)
RI2 =~ 1*y2_(u)

RI1 ~~ RI1 + RI2
RI2 ~~ RI2
"

## -----------------------------------------------------------------------------
riclpm <- lessTemplates::CLPM(model = model, 
                              data = dataset)
fit <- lavaan::sem(model = clpm$model,
                   data = clpm$data)
coef(fit)

## -----------------------------------------------------------------------------
summary(fit)

## -----------------------------------------------------------------------------
dataset <- simulateCLPM(seed = 123)
head(dataset)

## -----------------------------------------------------------------------------
model <- "
# autoregressive and cross-lagged effects
eta1_(u) ~ a11*eta1_(u-1) + a12*eta2_(u-1)
eta2_(u) ~ a21*eta1_(u-1) + a22*eta2_(u-1)

# covariances between latent variables
eta1_(u) ~~ v11*eta1_(u) + 0*eta2_(u)
eta2_(u) ~~ v22*eta2_(u)

# measurements
eta1_(u) =~ 1*y1_(u) + l21_(u)*y2_(u) + l31_(u)*y3_(u)
eta2_(u) =~ 1*y4_(u) + l52_(u)*y5_(u) + l62_(u)*y6_(u)

# measurement error variances
y1_(u) ~~ v1*y1_(u)
y2_(u) ~~ v2*y2_(u)
y3_(u) ~~ v3*y3_(u)
y4_(u) ~~ v4*y4_(u)
y5_(u) ~~ v5*y5_(u)
y6_(u) ~~ v6*y6_(u)
"

## -----------------------------------------------------------------------------
clpm <- lessTemplates::CLPM(model = model, 
                            data = dataset)
fit <- lavaan::sem(model = clpm$model,
                   data = clpm$data)
coef(fit)

## -----------------------------------------------------------------------------
summary(fit)

## -----------------------------------------------------------------------------
clpmTransform <- transformCLPM(CLPM = clpm, 
                               parameters = c("l21_(u)", "l31_(u)",
                                              "l52_(u)", "l62_(u)"),
                               transformation = "measurementInvariance")

## -----------------------------------------------------------------------------
cat(clpmTransform$transformation)

## ----include=FALSE------------------------------------------------------------
library(lessSEM)
lassoFit <- lasso(lavaanModel = fit, 
                  regularized = clpmTransform$regularized, 
                  nLambdas = 50,
                  modifyModel = modifyModel(transformations = clpmTransform$transformation))

## ----eval = FALSE-------------------------------------------------------------
#  library(lessSEM)
#  lassoFit <- lasso(lavaanModel = fit,
#                    regularized = clpmTransform$regularized,
#                    nLambdas = 50,
#                    modifyModel = modifyModel(transformations = clpmTransform$transformation))

## -----------------------------------------------------------------------------
coef(lassoFit, criterion = "BIC")

## -----------------------------------------------------------------------------
model <- "
# autoregressive and cross-lagged effects
eta1_(u) ~ a11_(u)*eta1_(u-1) + a12_(u)*eta2_(u-1)
eta2_(u) ~ a21_(u)*eta1_(u-1) + a22_(u)*eta2_(u-1)

# covariances between latent variables
eta1_(u) ~~ v11*eta1_(u) + 0*eta2_(u)
eta2_(u) ~~ v22*eta2_(u)

# measurements
eta1_(u) =~ 1*y1_(u) + l21*y2_(u) + l31*y3_(u)
eta2_(u) =~ 1*y4_(u) + l52*y5_(u) + l62_(u)*y6_(u)

# measurement error variances
y1_(u) ~~ v1*y1_(u)
y2_(u) ~~ v2*y2_(u)
y3_(u) ~~ v3*y3_(u)
y4_(u) ~~ v4*y4_(u)
y5_(u) ~~ v5*y5_(u)
y6_(u) ~~ v6*y6_(u)
"

## -----------------------------------------------------------------------------
clpm <- lessTemplates::CLPM(model = model, 
                            data = dataset)
fit <- lavaan::sem(model = clpm$model,
                   data = clpm$data)
coef(fit)

## -----------------------------------------------------------------------------
summary(fit)

## -----------------------------------------------------------------------------
clpmTransform <- transformCLPM(CLPM = clpm, 
                               parameters = c("a11_(u)", "a12_(u)",
                                              "a21_(u)", "a22_(u)"),
                               transformation = "changepoint")

## -----------------------------------------------------------------------------
cat(clpmTransform$transformation)

## ----include=FALSE------------------------------------------------------------
library(lessSEM)
lassoFit <- lasso(lavaanModel = fit, 
                  regularized = clpmTransform$regularized, 
                  nLambdas = 50,
                  modifyModel = modifyModel(transformations = clpmTransform$transformation))

## ----eval = FALSE-------------------------------------------------------------
#  library(lessSEM)
#  lassoFit <- lasso(lavaanModel = fit,
#                    regularized = clpmTransform$regularized,
#                    nLambdas = 50,
#                    modifyModel = modifyModel(transformations = clpmTransform$transformation))

## -----------------------------------------------------------------------------
coef(lassoFit, criterion = "BIC")

