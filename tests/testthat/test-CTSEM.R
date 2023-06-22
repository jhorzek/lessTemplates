test_that("CTSEM works", {

  library(ctsemOMX)
  library(lessTemplates)
  library(lessSEM)
  library(lavaan)

  set.seed(123)

  data(AnomAuth)
  data <- ctWideToLong(datawide = AnomAuth,
                       Tpoints= 5,
                       n.manifest=2,
                       manifestNames = c("Y1", "Y2"))

  data <- ctDeintervalise(datalong = data, id='id', dT='dT')
  colnames(data) <- c("person", "time", "Y1", "Y2")
  data <- as.data.frame(data)
  data <- data[!(is.na(data$Y1) & is.na(data$Y2)),]

  model <- "
d_eta1(t) ~ eta1(t) + eta2(t)
d_eta2(t) ~ eta1(t) + eta2(t)

d_eta1(t) ~~ d_eta1(t) + d_eta2(t)
d_eta2(t) ~~ d_eta2(t)

eta1(t) =~ 1*Y1(t)
eta2(t) =~ 1*Y2(t)

Y1(t) ~~ 0*Y1(t)
Y2(t) ~~ 0*Y2(t)

Y1(t) ~ m1*1
Y2(t) ~ m2*1

eta1(0) ~ 1
eta2(0) ~ 1
"

  ctsem <- lessTemplates::CTSEM(model = model,
                                data = data)

  fit <- bfgs(lavaanModel = ctsem$lavaanModel,
              modifyModel = modifyModel(transformations = ctsem$transformation,
                                        transformationList = ctsem$transformationList))

  # comparison
  AnomAuthmodel <- ctModel(LAMBDA = matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2),
                           Tpoints = 5, n.latent = 2, n.manifest = 2, MANIFESTVAR=diag(0, 2), TRAITVAR = NULL)
  AnomAuthfit <- ctFit(AnomAuth, AnomAuthmodel, useOptimizer = TRUE)

  testthat::expect_equal(abs(AnomAuthfit$mxobj$fitfunction$result[[1]] - fit@fits$m2LL[1])/
                           AnomAuthfit$mxobj$fitfunction$result[[1]] < 1e-5, TRUE)
  testthat::expect_equal(all(abs(AnomAuthfit$mxobj$DRIFT$values - ctsem$transformationList$DRIFT) < 1e-1), TRUE)
  testthat::expect_equal(all(abs(AnomAuthfit$mxobj$DIFFUSION$result - ctsem$transformationList$DIFFUSION) < 1e-1), TRUE)

  ## Testing fixed values:
  model <- "
d_eta1(t) ~ .2*eta1(t) + -0.4*eta2(t)
d_eta2(t) ~ eta1(t) + eta2(t)

d_eta1(t) ~~ 1*d_eta1(t) + .1*d_eta2(t)
d_eta2(t) ~~ d_eta2(t)

eta1(t) =~ 1*Y1(t)
eta2(t) =~ 1*Y2(t)

Y1(t) ~~ 0*Y1(t)
Y2(t) ~~ 0*Y2(t)

Y1(t) ~ m1*1
Y2(t) ~ m2*1

eta1(0) ~ 1
eta2(0) ~ 1
"
  ctsem <- lessTemplates::CTSEM(model = model,
                                data = data)

  fit <- bfgs(lavaanModel = ctsem$lavaanModel,
              modifyModel = modifyModel(transformations = ctsem$transformation,
                                        transformationList = ctsem$transformationList))

  testthat::expect_equal(ctsem$transformationList$DRIFT[1,1], .2)
  testthat::expect_equal(ctsem$transformationList$DRIFT[1,2], -.4)
  testthat::expect_equal(ctsem$transformationList$DIFFUSION[1,1], 1)
  testthat::expect_equal(ctsem$transformationList$DIFFUSION[1,2], .1)
})
