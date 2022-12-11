test_that("definition variables work", {
  #'
  # We will demonstrate the use of definition variables with a latent growth curve
  # model. We assume that the six measurment occasions took place at subject-specific
  # time points.


  #### Data simulation #####
  # You can ignore the following; we just simulate some data to be used
  # for our model.
  ## Population parameters ##
  set.seed(123)
  intercept_mu <- 0
  intercept_sigma <- 1
  slope_mu <- .3
  slope_sigma <- 1

  ## data set ##
  N <- 50
  intercepts <- rnorm(n = N,
                      mean = intercept_mu,
                      sd = intercept_sigma)
  slopes <- rnorm(n = N,
                  mean = slope_mu,
                  sd = slope_sigma)
  times <- as.data.frame(
    matrix(seq(0,5,1),
           nrow = N,
           ncol = 6,
           byrow = TRUE,
           dimnames = list(NULL, paste0("t", 0:5))) +
      cbind(0,matrix(round(runif(n = N*5, min = -.2,max = .2),2),
                     nrow = N,
                     ncol = 5,
                     byrow = TRUE)) # we add some jitter to make the times person-specific
  )

  lgcData <- matrix(NA,
                    nrow = N,
                    ncol = ncol(times),
                    dimnames = list(NULL, paste0("x", 0:5)))

  for(i in 1:N){
    lgcData[i,] <- intercepts[i] + unlist(times[i,])* slopes[i] + rnorm(ncol(lgcData),0,.3)
  }
  lgcData <- as.data.frame(lgcData)

  #### Let's have a look at the data ####
  # lgcData is a data.frame with the observed values at the six measurement
  # occasions:
  head(lgcData)

  # times is a data.frame with the time points at which the six measurements
  # took place for each person. Note that these differ between individuals.
  # These are going to be our definition variables.
  head(times)

  #### Set up the base model ####

  model <- "
int =~ 1*x0 + 1*x1 + 1*x2 + 1*x3 + 1*x4 + 1*x5
slope =~ t0*x0 + t1*x1 + t2*x2 + t3*x3 + t4*x4 + t5*x5

int ~ intMean*1
slope ~ slopeMean*1

int ~~ intVar*int + intSlopeCov*slope
slope ~~ slopeVar*slope

x0 ~~ v*x0; x1 ~~ v*x1; x2 ~~ v*x2; x3 ~~ v*x3; x4 ~~ v*x4; x5 ~~ v*x5

x0 ~ 0*1; x1 ~ 0*1; x2 ~ 0*1; x3 ~ 0*1; x4 ~ 0*1; x5 ~ 0*1
"

  # The most important part of our model definition is that we used the
  # modifiers t0, t1, t2, t3, t4, and t5. These have the same names as our
  # defintion variables in the times variable. lessTemplates with
  # replace the modifiers with the values found in the times variable.

  #### Set up the model with definition variables ####

  dvarModel <- SEMWithDefinitionVariables(lavaanSyntax = model,
                                          definitionVariables = times,
                                          data = lgcData)

  # fit with lessSEM
  library(lessSEM)
  dvarFit <- lessSEM::bfgs(lavaanModel = dvarModel)
  dvarFit@parameters

  library(OpenMx)
  latents <- c("int", "slope")
  manifests <- colnames(lgcData)

  mxMod <- mxModel(type = "RAM",
                   manifestVars = manifests,
                   latentVars = latents,
                   mxPath(from = "int", to = manifests, values = 1, free = FALSE),
                   mxPath(from = "slope", to = manifests, labels = paste0("data.t", 0:5), free = FALSE),
                   mxPath(from = latents, to = latents, arrows = 2, free = TRUE, connect = "unique.pairs", labels = c("intVar", "intSlopeCov", "slopeVar")),
                   mxPath(from = manifests, to = manifests, labels = "v", arrows = 2),
                   mxPath(from = "one", to = latents, values = 0, free = TRUE, labels = c("intMean", "slopeMean")),
                   mxPath(from = "one", to = manifests, values = 0, free = FALSE),
                   mxData(cbind(lgcData, times), type = "raw"))
  fitMod <- mxRun(mxMod)
  mxPar <- omxGetParameters(fitMod)

  testthat::expect_equal(all(abs(unlist(dvarFit@parameters[,names(mxPar)]) - mxPar) < 1e-1), TRUE)
  testthat::expect_equal(abs(fitMod@fitfunction$result[[1]] - dvarFit@fits$m2LL[1]) < 1e-2, TRUE)

})
