test_that("make ARCL", {
  library(mvtnorm)
  set.seed(1234)

  ### simulate data ###
  tps <- 10
  N <- 100
  A <- matrix(c(.6,.25,
                -.3,.65),2,2,TRUE)

  L <- matrix(0, nrow = 6, ncol = 2)
  L[1:3,1] <- c(1,.9,.9)
  L[4:6,2] <- c(1,.9,.9)

  mVar <- matrix(0, nrow(L), nrow(L))
  diag(mVar) <- 1.2 - apply(L^2,1,sum)

  # create data set
  data <- expand.grid(occasion = 1:tps,
                      id = 1:N
  )
  data <- cbind(data,
                matrix(NA,
                       nrow = nrow(data),
                       ncol = 6,
                       dimnames = list(NULL, paste0("y", 1:6)))
  )

  for(i in 1:N){

    for(tp in 1:tps){

      if(tp == 1){
        lvCov <- diag(1, nrow(A))

        par_lCov <- c(lvCov[1,1], lvCov[2,2])
        names(par_lCov) <- paste0(c("lvar1_", "lvar2_"), tp)

        eta <- t(mvtnorm::rmvnorm(n = 1,
                                  mean = c(0,0),
                                  sigma = lvCov))

      }else{

        lvCov <- matrix(0, nrow(A), nrow(A))
        diag(lvCov) <- 1-apply(A^2,1,sum)

        eta <- A%*%eta + t(mvtnorm::rmvnorm(n = 1,
                                            mean = c(0,0),
                                            sigma = lvCov))
      }
      data[data$occasion==tp & data$id == i,paste0("y",1:6)] <- L%*%eta + t(mvtnorm::rmvnorm(n = 1,
                                                                                             mean = rep(0,6),
                                                                                             sigma = mVar))
    }
  }

  library(lavaan)
  nmanifest <- 6
  nlatent <- 2

  manifestIntercept <- matrix(paste0("mInt", 1:nmanifest), nmanifest, 1)
  loadings <- matrix("0", nmanifest, nlatent)
  loadings[1:3,1] <- c("1", "l21","l31")
  loadings[4:6,2] <- c("1", "l52","l62")

  measurementError <- matrix("0", nmanifest, nmanifest)
  diag(measurementError) <- paste0("mvar", 1:nmanifest)

  latentIntercept <- matrix("0", nlatent, 1)
  ARCL <- matrix(c("a11", "a12",
                   "a21", "a22"),2,2,T)

  latentError <- matrix(c("lvar1", "0",
                          "0", "lvar2"),2,2,T)

  # estimate model with all parameters being equal over time

  arclModel <- makeARCL(data = data,
                        latentIntercepts = latentIntercept,
                        loadings = loadings,
                        measurementError = measurementError,
                        manifestIntercepts = manifestIntercept,
                        ARCL = ARCL,
                        latentResiduals = latentError)

  lavaanFit <- sem(model = arclModel$model,
                   data = arclModel$data)

  coef(lavaanFit)[unique(names(coef(lavaanFit)))]

  # make ARCL occasion specific

  data$lag <- data$occasion
  arclModel2 <- makeARCL(data = data,
                         latentIntercepts = latentIntercept,
                         loadings = loadings,
                         measurementError = measurementError,
                         manifestIntercepts = manifestIntercept,
                         ARCL = ARCL,
                         latentResiduals = latentError,
                         lagSpecific = "ARCL")

  lavaanFit2 <- sem(model = arclModel2$model,
                    data = arclModel2$data)
  coef(lavaanFit2)[unique(names(coef(lavaanFit2)))]



})
