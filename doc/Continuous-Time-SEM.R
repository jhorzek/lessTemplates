## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, include=FALSE-----------------------------------------------------
library(lessTemplates)
library(lavaan)
library(lessSEM)
library(ctsemOMX)
set.seed(123)

## -----------------------------------------------------------------------------
library(ctsemOMX)
data("AnomAuth")
head(AnomAuth)

## -----------------------------------------------------------------------------
# from long to wide
data <- ctWideToLong(datawide = AnomAuth,
                     Tpoints= 5,
                     n.manifest=2,
                     manifestNames = c("Y1", "Y2"))
data <- ctDeintervalise(datalong = data, 
                        id='id', 
                        dT='dT')
colnames(data) <- c("person", "time", "Y1", "Y2")
data <- as.data.frame(data)
data <- data[!(is.na(data$Y1) & is.na(data$Y2)),]
head(data)

## -----------------------------------------------------------------------------
model <- "
# Specify the latent dynamics in differential equation model notation:
d_eta1(t) ~ eta1(t) + eta2(t)
d_eta2(t) ~ eta1(t) + eta2(t)

# Covariances (Wiener process)
d_eta1(t) ~~ d_eta1(t) + d_eta2(t)
d_eta2(t) ~~ d_eta2(t)

# Latent intercepts
eta1(0) ~ 1
eta2(0) ~ 1

# Measurement model
eta1(t) =~ 1*Y1(t)
eta2(t) =~ 1*Y2(t)

# Manifest variances and covariances
Y1(t) ~~ 0*Y1(t) + 0*Y2(t)
Y2(t) ~~ 0*Y2(t)

# Manifest intercepts
Y1(t) ~ m1*1
Y2(t) ~ m2*1
"

## -----------------------------------------------------------------------------
ctsem <- CTSEM(model = model, 
               data = data)
names(ctsem)

## -----------------------------------------------------------------------------
summary(ctsem$lavaanModel)

## -----------------------------------------------------------------------------
cat(ctsem$transformation)
# the following list contains elements referenced in the transformations:
print(ctsem$transformationList)

## ----include=FALSE------------------------------------------------------------
library(lessSEM)

fit <- bfgs(lavaanModel = ctsem$lavaanModel, 
            modifyModel = modifyModel(transformations = ctsem$transformation,
                                      transformationList = ctsem$transformationList))

## ----eval=FALSE---------------------------------------------------------------
#  library(lessSEM)
#  
#  fit <- bfgs(lavaanModel = ctsem$lavaanModel,
#              modifyModel = modifyModel(transformations = ctsem$transformation,
#                                        transformationList = ctsem$transformationList))

## -----------------------------------------------------------------------------
coef(fit)

## -----------------------------------------------------------------------------
# DRIFT:
print(ctsem$transformationList$DRIFT)

# DIFFUSION:
print(ctsem$transformationList$DIFFUSION)

## -----------------------------------------------------------------------------
AnomAuthmodel <- ctModel(LAMBDA = matrix(c(1, 0, 0, 1), 
                                         nrow = 2, 
                                         ncol = 2), 
                         Tpoints = 5, 
                         n.latent = 2, 
                         n.manifest = 2, 
                         MANIFESTVAR=diag(0, 2), 
                         TRAITVAR = NULL) 
AnomAuthfit <- ctFit(AnomAuth, 
                     AnomAuthmodel)
summary(AnomAuthfit)

## -----------------------------------------------------------------------------
modelRI <- "
# Specify the latent dynamics in differential equation model notation:
d_eta1(t) ~ eta1(t) + eta2(t)
d_eta2(t) ~ eta1(t) + eta2(t)

# Covariances (Wiener process)
d_eta1(t) ~~ d_eta1(t) + d_eta2(t)
d_eta2(t) ~~ d_eta2(t)

# Latent intercepts
eta1(0) ~ 1
eta2(0) ~ 1

# Measurement model
eta1(t) =~ 1*Y1(t)
eta2(t) =~ 1*Y2(t)

# Manifest variances and covariances
Y1(t) ~~ 0*Y1(t) + 0*Y2(t)
Y2(t) ~~ 0*Y2(t)

# Manifest intercepts
Y1(t) ~ m1*1
Y2(t) ~ m2*1

# Add random intercepts
RI1 =~ 1*Y1(t)
RI2 =~ 1*Y2(t)

RI1 ~~ RI1 + RI2
RI2 ~~ RI2
"

## ----include=FALSE------------------------------------------------------------
ctsemRI <- CTSEM(model = modelRI, 
                 data = data)
fitRI <- bfgs(lavaanModel = ctsemRI$lavaanModel, 
              modifyModel = modifyModel(transformations = ctsemRI$transformation,
                                        transformationList = ctsemRI$transformationList),
              control = controlBFGS(breakOuter = 1e-10))

## ----eval=FALSE---------------------------------------------------------------
#  ctsemRI <- CTSEM(model = modelRI,
#                   data = data)
#  fitRI <- bfgs(lavaanModel = ctsemRI$lavaanModel,
#                modifyModel = modifyModel(transformations = ctsemRI$transformation,
#                                          transformationList = ctsemRI$transformationList),
#                control = controlBFGS(breakOuter = 1e-10))

## -----------------------------------------------------------------------------
coef(fitRI)

# DRIFT:
print(ctsemRI$transformationList$DRIFT)

# DIFFUSION:
print(ctsemRI$transformationList$DIFFUSION)

## -----------------------------------------------------------------------------
AIC(fit)
AIC(fitRI)

## -----------------------------------------------------------------------------
AnomAuthmodel <- ctModel(LAMBDA = matrix(c(1, 0, 0, 1), 
                                         nrow = 2, 
                                         ncol = 2), 
                         Tpoints = 5, 
                         n.latent = 2, 
                         n.manifest = 2, 
                         MANIFESTVAR=diag(0, 2), 
                         TRAITVAR = NULL,
                         MANIFESTTRAITVAR = "auto" # add random intercepts
) 
AnomAuthfit <- ctFit(AnomAuth, 
                     AnomAuthmodel)
print(AnomAuthfit$mxobj$DRIFT$values)

