#' CTSEM
#'
#' Set up a continuous time structural equation model with lessSEM
#'
#' @param model syntax to specify the model
#' @param data data set in long format with person and time for each observation
#' @examples
#' library(ctsemOMX)
#' library(lessTemplates)
#' library(lessSEM)
#' library(lavaan)
#'
#' data(AnomAuth)
#' data <- ctWideToLong(datawide = AnomAuth,
#'                      Tpoints= 5,
#'                      n.manifest=2,
#'                      manifestNames = c("Y1", "Y2"))
#'
#' data <- ctDeintervalise(datalong = data, id='id', dT='dT')
#' colnames(data) <- c("person", "time", "Y1", "Y2")
#' data <- as.data.frame(data)
#' data <- data[!(is.na(data$Y1) & is.na(data$Y2)),]
#'
#' model <- "
#' d_eta1(t) ~ eta1(t) + eta2(t)
#' d_eta2(t) ~ eta1(t) + eta2(t)
#'
#' d_eta1(t) ~~ d_eta1(t) + d_eta2(t)
#' d_eta2(t) ~~ d_eta2(t)
#'
#' eta1(t) =~ 1*Y1(t)
#' eta2(t) =~ 1*Y2(t)
#'
#' Y1(t) ~~ 0*Y1(t)
#' Y2(t) ~~ 0*Y2(t)
#'
#' Y1(t) ~ m1*1
#' Y2(t) ~ m2*1
#'
#' eta1(0) ~ 1
#' eta2(0) ~ 1
#' "
#'
#' ctsem <- lessTemplates::CTSEM(model = model,
#'                                     data = data)
#'
#' fit <- bfgs(lavaanModel = ctsem$lavaanModel,
#'             modifyModel = modifyModel(transformations = ctsem$transformation,
#'                                       transformationList = ctsem$transformationList))
#'
#' fit@parameters[,sort(fit@parameterLabels)]
#'
#' # comparison
#' AnomAuthmodel <- ctModel(LAMBDA = matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2),
#'                          Tpoints = 5, n.latent = 2, n.manifest = 2, MANIFESTVAR=diag(0, 2), TRAITVAR = NULL)
#' AnomAuthfit <- ctFit(AnomAuth, AnomAuthmodel, useOptimizer = T)
#' summary(AnomAuthfit)
#' @export
CTSEM <- function(model,
                  data){

  cat("\nSetting up a continuous time structural equation model.\n")

  if(!all(c("person", "time") %in% colnames(data)))
    stop("Could not find the columns person and time in the data set.")
  if(!is(object = data, class2 = "data.frame"))
    data <- as.data.frame(data)

  dataCTSEM <- lessTemplates:::.prepareDataCTSEM(data = data)

  arclData <- dataCTSEM$data

  # remove unnecessary white space
  syntax <- lessTemplates:::.reduceSyntax(syntax = model)
  syntax <- lessTemplates:::.removeWhitespace(syntax = syntax)
  syntax <- lessTemplates:::.makeSingleLine(syntax = syntax)

  # check if all dynamics are time dependent
  lessTemplates:::.checkDynamics(syntax = syntax)

  # find the names of all variables
  variableNames <- lessTemplates:::.getVariableNamesCTSEM(syntax = syntax)

  latents <- list(
    occasionDependent = variableNames$occasionDependent[!variableNames$occasionDependent %in% colnames(data)],
    fixed = variableNames$fixed[!variableNames$fixed %in% colnames(data)]
  )

  manifests <- list(
    occasionDependent = variableNames$occasionDependent[variableNames$occasionDependent %in% colnames(data)],
    fixed = variableNames$fixed[variableNames$fixed %in% colnames(data)]
  )

  # set up discrete time model as basis
  arclModel <- lessTemplates:::.getDiscreteBasis(syntax = syntax,
                                                 latents = latents,
                                                 manifests = manifests)

  ctMatrices <- lessTemplates:::.getCTSEMMatrices(syntax = syntax,
                                                  arclModel = arclModel)

  clpm <- lessTemplates::CLPM(model = arclModel$syntax,
                              data = arclData[,colnames(arclData) != "time"],
                              silent = TRUE,
                              meanstructure = grepl(pattern = "*1", x = arclModel$syntax))

  transformations <- lessTemplates:::.getCTSEMTransformations(arclModel = arclModel,
                                                              dataCTSEM = dataCTSEM,
                                                              ctMatrices = ctMatrices)

  fit_lavaan <- lavaan::sem(model = clpm$model,
                            data = clpm$data,
                            do.fit = FALSE,
                            missing = "ml",
                            meanstructure = TRUE)
  return(list(
    lavaanModel = fit_lavaan,
    transformation = transformations$transformation,
    transformationList = transformations$transformationList,
    internal = list(
      model = clpm$model,
      data = clpm$data,
      clpm = clpm,
      arclModel = arclModel
    )
  ))
}

.checkDynamics <- function(syntax){

  syntax_t <- gsub(pattern = "[a-zA-Z0-9]+\\*|[a-zA-Z0-9]+\\*",
                   replacement = "",
                   x = syntax)
  # only keep dynamics
  syntax_t <- syntax_t[grepl(pattern = "^d[0-9]*_",
                             x = syntax_t) & # is dynamic
                         grepl(pattern = "~",
                               x = syntax_t) & # is regression
                         !grepl(pattern = "~~",
                                x = syntax_t) # is not covariance
  ]

  # split at operators
  variables <- unlist(stringr::str_split(string = syntax_t,
                                         pattern = "\\+|~"))
  # check that all variables are time-dependent

  if(!all(grepl(pattern = "\\(t\\)$", x = variables)))
    stop("Some of your equations seem to include variables which are not time dependent: ",
         paste0(variables[grepl(pattern = "\\(t\\)$", x = variables)], collapse = ", ")
    )

}

#' .getVariableNamesCLPM
#'
#' extracts the names of the variables from the syntax
#' @param syntax string
#' @return list with names of variables
#' @keywords internal
.getVariableNamesCTSEM <- function(syntax){

  # remove all parameters
  syntax_t <- gsub(pattern = "[.a-zA-Z0-9_\\-]+\\*",
                   replacement = "",
                   x = syntax)
  # remove means
  syntax_t <- gsub(pattern = "~1",
                   replacement = "",
                   x = syntax_t)
  # split at operators
  variableNames <- unlist(stringr::str_split(string = syntax_t,
                                             pattern = "\\+|=~|~~|~"))
  # remove time indices
  isTimeDependent <- grepl(pattern = "\\(t\\)",
                           x = variableNames)
  variableNames <- gsub(pattern = "\\(t\\)",
                        replacement = "",
                        x = variableNames)
  names(isTimeDependent) <- variableNames

  return(
    list(occasionDependent = unique(variableNames[isTimeDependent]),
         fixed = unique(variableNames[!isTimeDependent]))
  )
}

.getDiscreteBasis <- function(syntax, latents, manifests){

  # keep all measurement equations
  measurementEquations <- syntax[grepl(pattern = "=~", x = syntax)]

  # replace dynamics in measurement equations
  measurementEquations <- stringr::str_replace_all(string = measurementEquations,
                                                   pattern = "\\(t\\)",
                                                   replacement = "_\\(u\\)")

  dynamicLatents <- latents$occasionDependent

  # remove all differential statements
  dynamicLatents <- unique(stringr::str_remove(string = dynamicLatents,
                                               pattern = "^d_"))

  # create names used in arcl model for dynamics:
  arclDynamics <- dynamicLatents

  arcls <- matrix(paste0("arcl_",
                         rep(1:length(arclDynamics), length(arclDynamics)),
                         "_",
                         rep(1:length(arclDynamics), each = length(arclDynamics)),
                         "_(u)"),
                  length(arclDynamics),
                  length(arclDynamics),
                  dimnames = list(paste0(arclDynamics, "_(u)"),
                                  paste0(arclDynamics, "_(u-1)")))

  # set up all discrete time arcl effects:
  arclSyntax <- c()
  for(i in 1:nrow(arcls)){

    arclSyntax <- c(arclSyntax,
                    paste0(rownames(arcls)[i], "~",
                           paste0(paste0(arcls[i,], "*", colnames(arcls)), collapse = " + ")
                    )
    )
  }

  # set up discrete time latent covariances
  covs <- matrix(paste0("cov_",
                        rep(1:length(arclDynamics), length(arclDynamics)),
                        "_",
                        rep(1:length(arclDynamics), each = length(arclDynamics)),
                        "_(u)"),
                 length(arclDynamics),
                 length(arclDynamics),
                 dimnames = list(paste0(arclDynamics, "_(u)"),
                                 paste0(arclDynamics, "_(u)")))

  for(i in 1:nrow(covs)){
    for(j in i:nrow(covs)){
      covs[i,j] <- covs[j,i]
    }
  }
  # set up all discrete time covariances:
  for(i in 1:nrow(covs)){

    arclSyntax <- c(arclSyntax,
                    paste0(rownames(covs)[i], "~~",
                           paste0(paste0(covs[i,], "*", colnames(covs)), collapse = " + ")
                    )
    )
  }


  # set up all discrete time arcl measurements
  for(i in 1:nrow(arcls)){
    tmpName <- stringr::str_replace(string = rownames(arcls)[i],
                                    pattern = "\\(",
                                    replacement = "\\\\(")
    tmpName <- stringr::str_replace(string = tmpName,
                                    pattern = "\\)",
                                    replacement = "\\\\)")
    measures <- grepl(pattern = paste0("^", tmpName),
                      x = measurementEquations)
    if(any(measures)){
      arclSyntax <- c(arclSyntax,
                      measurementEquations[measures]
      )
      measurementEquations <- measurementEquations[!measures]
    }else{

      arclSyntax <- c(
        arclSyntax,
        paste0(rownames(arcls)[i], "=~",
               paste0(paste0("0*", c(manifests$occasionDependent, manifests$fixed)), collapse = " + ")
        )
      )

    }
  }

  # add remaining measurements:
  arclSyntax <- c(arclSyntax,
                  measurementEquations)

  arclSyntax <- paste0(arclSyntax,
                       collapse = "\n")

  # add remaining manifest covariances
  measurementCovs <- syntax[grepl(pattern = "~~",
                                  x = syntax) &
                              !grepl(pattern = paste0(paste0("^", latents$occasionDependent), collapse = "|"),
                                     x = syntax)
  ] # remove latent dynamic covs; these are captured in the diffusion

  # replace dynamics in measurement equations
  measurementCovs <- stringr::str_replace_all(string = measurementCovs,
                                              pattern = "\\(t\\)",
                                              replacement = "_\\(u\\)")
  # add
  arclSyntax <- c(arclSyntax,
                  measurementCovs)

  # add manifest intercepts
  intercepts <- grepl(pattern = "[\\*]*1$", x = syntax) &
    grepl(pattern = "~", x = syntax) &
    !grepl(pattern = "~~", x = syntax)

  measurementInts <- syntax[intercepts &
                              !grepl(pattern = paste0(paste0("^", latents$occasionDependent), collapse = "|"),
                                     x = syntax)]
  measurementInts <- stringr::str_replace_all(string = measurementInts,
                                              pattern = "\\(t\\)",
                                              replacement = "_\\(u\\)")
  arclSyntax <- c(arclSyntax,
                  measurementInts)

  # check for initial means of latent variables
  initialMeans <- syntax[intercepts &
                           grepl(pattern = paste0(paste0("^", latents$occasionDependent, "\\(0\\)"), collapse = "|"),
                                 x = syntax)]
  initialMeans <- stringr::str_replace_all(string = initialMeans,
                                           pattern = "\\(0\\)",
                                           replacement = "_\\(1\\)")
  arclSyntax <- c(arclSyntax,
                  initialMeans)

  # combine
  arclSyntax <- paste0(arclSyntax,
                       collapse = "\n")

  return(list(syntax = arclSyntax,
              arcl = arcls,
              covs = covs,
              dynamicLatents = dynamicLatents
  ))
}

.prepareDataCTSEM <- function(data){
  times <- data$time
  times <- unique(sort(times))
  timeIntervals <- times - c(times[1],times[1:(length(times)-1)])
  occasions <- 1:length(times)

  data$occasion <- NA

  for(i in 1:length(times)){
    sel <- which(data$time == times[i])
    data$occasion[sel] <- occasions[i]
  }

  return(
    list(
      data = data,
      timetable = data.frame(time = times,
                             occasion = occasions,
                             timeInterval = timeIntervals)
    )
  )
}

.getCTSEMMatrices <- function(syntax, arclModel){

  dependent <- paste0("d_", arclModel$dynamicLatents)
  predictor <- arclModel$dynamicLatents

  # Set up the matrices; the individual values will be replaced later on
  DRIFT <- matrix("0.0",
                  nrow = length(predictor),
                  ncol = length(predictor),
                  dimnames = list(dependent, predictor)
  )

  for(i in 1:nrow(DRIFT)){
    predicted_in <- grepl(pattern = paste0("^", rownames(DRIFT)[i], "\\(t\\)~[0-9a-zA-Z.\\-]+"), x = syntax)
    if(!any(predicted_in)) next

    rhs <- stringr::str_remove(syntax[predicted_in],
                               paste0("^", rownames(DRIFT)[i], "\\(t\\)~"))

    for(j in 1:ncol(DRIFT)){
      # check if variable j predicts the change in i:
      if(grepl(pattern = paste0(colnames(DRIFT)[j], "\\(t\\)"), x = rhs)){
        element <- stringr::str_extract(string = rhs,
                                        pattern = paste0("[.0-9a-zA-Z\\*\\-]*",colnames(DRIFT)[j], "\\(t\\)"))
        if(grepl(pattern = "\\*", x = element)){
          DRIFT[i,j] <- stringr::str_extract(string = element,
                                             pattern = "^[\\-]*[.0-9a-zA-Z]+")
        }else{
          DRIFT[i,j] <- paste0("drift_",rownames(DRIFT)[i],"_", colnames(DRIFT)[j])
        }
      }
    }
  }

  logDiagDIFFUSION <- matrix("0.0",
                             nrow = length(dependent),
                             ncol = length(dependent),
                             dimnames = list(dependent, dependent)
  )

  for(i in 1:nrow(logDiagDIFFUSION)){
    predicted_in <- grepl(pattern = paste0("^", rownames(logDiagDIFFUSION)[i], "\\(t\\)~~[0-9a-zA-Z.\\-]+"), x = syntax)
    if(!any(predicted_in)) next

    rhs <- stringr::str_remove(syntax[predicted_in],
                               paste0("^", rownames(logDiagDIFFUSION)[i], "\\(t\\)~~"))

    for(j in 1:ncol(logDiagDIFFUSION)){
      # check if variable j predicts the change in i:
      if(grepl(pattern = paste0(colnames(logDiagDIFFUSION)[j], "\\(t\\)"), x = rhs)){
        element <- stringr::str_extract(string = rhs,
                                        pattern = paste0("[.0-9a-zA-Z\\*\\-]*",colnames(logDiagDIFFUSION)[j], "\\(t\\)"))
        if(grepl(pattern = "\\*", x = element)){
          if((i == j) & grepl(pattern = "^[.0-9\\-]+", x = element)){
            # in case of diagonal values, the parameter will be the log of the actual
            # value. This prevents negative variances.
            logDiagDIFFUSION[i,j] <- log(as.numeric(stringr::str_extract(string = element,
                                                                         pattern = "^[\\-]*[.0-9a-zA-Z]+")))
          }else{
            logDiagDIFFUSION[i,j] <- stringr::str_extract(string = element,
                                                          pattern = "^[\\-]*[.0-9a-zA-Z]+")
          }
        }else{
          if(i == j){
            # in case of diagonal values, the parameter will be the log of the actual
            # value. This prevents negative variances.
            logDiagDIFFUSION[i,j] <- paste0("log_diffusion_",rownames(logDiagDIFFUSION)[i],"_", colnames(logDiagDIFFUSION)[j])
          }else{
            logDiagDIFFUSION[i,j] <- paste0("diffusion_",rownames(logDiagDIFFUSION)[i],"_", colnames(logDiagDIFFUSION)[j])
          }
        }
        logDiagDIFFUSION[j,i] <- logDiagDIFFUSION[i,j]
      }
    }
  }

  return(
    list(DRIFT = DRIFT,
         logDiagDIFFUSION = logDiagDIFFUSION)
  )
}

.getCTSEMTransformations <- function(arclModel, dataCTSEM, ctMatrices){

  transformationList <- list(
    DRIFT = matrix(0, nrow = nrow(ctMatrices$DRIFT), ncol = ncol(ctMatrices$DRIFT)),
    logDiagDIFFUSION = matrix(0, nrow = nrow(ctMatrices$logDiagDIFFUSION), ncol = ncol(ctMatrices$logDiagDIFFUSION)),
    DIFFUSION = matrix(0, nrow = nrow(ctMatrices$logDiagDIFFUSION), ncol = ncol(ctMatrices$logDiagDIFFUSION)),
    driftHash = kronecker(diag(nrow(ctMatrices$DRIFT)), diag(nrow(ctMatrices$DRIFT)))
  ) # list with additional elements to speed up the computation

  # change values if user does not want to estimate all parameters:
  # for(i in 1:nrow(transformationList$DRIFT)){
  #   for(j in 1:ncol(transformationList$DRIFT)){
  #     if(grepl(pattern = "^[0-9]+", x = ctMatrices$DRIFT[i,j])){
  #       transformationList$DRIFT[i,j] <- as.numeric(ctMatrices$DRIFT[i,j])
  #     }
  #   }
  # }
  # for(i in 1:nrow(transformationList$logDiagDIFFUSION)){
  #   for(j in 1:ncol(transformationList$logDiagDIFFUSION)){
  #     if(grepl(pattern = "^[0-9]+", x = ctMatrices$logDiagDIFFUSION[i,j])){
  #       if(i == j){
  #         transformationList$logDiagDIFFUSION[i,j] <- log(as.numeric(ctMatrices$logDiagDIFFUSION[i,j]))
  #       }else{
  #         transformationList$logDiagDIFFUSION[i,j] <- as.numeric(ctMatrices$logDiagDIFFUSION[i,j])
  #       }
  #     }
  #   }
  # }

  lengthManif <- length(unique(dataCTSEM$timetable$timeInterval[dataCTSEM$timetable$timeInterval != 0]))
  manifelements <- vector("list", length = 2*lengthManif)
  names(manifelements) <- c(paste0("ARCL_", which(unique(dataCTSEM$timetable$timeInterval) != 0)),
                            paste0("LVCOV_", which(unique(dataCTSEM$timetable$timeInterval) != 0)))
  for(i in 1:length(manifelements)){
    if(grepl(pattern = "ARCL_", x = names(manifelements[i])))
      manifelements[[i]] <- matrix(0, nrow = nrow(ctMatrices$DRIFT), ncol = ncol(ctMatrices$DRIFT))
    if(grepl(pattern = "LVCOV_", x = names(manifelements[i])))
      manifelements[[i]] <- matrix(0, nrow = nrow(ctMatrices$logDiagDIFFUSION), ncol = ncol(ctMatrices$logDiagDIFFUSION))
  }

  transformationList <- c(transformationList, manifelements)

  parameters <- c() # names of parameters
  start <- c() # starting values
  for(i in 1:nrow(ctMatrices$DRIFT)){
    for(j in 1:ncol(ctMatrices$DRIFT)){

      if(grepl(pattern = "^[a-zA-Z]", x = ctMatrices$DRIFT[i,j])){
        parameters <- c(parameters,
                        ctMatrices$DRIFT[i,j])
        if(i == j){
          # we use the same starting values as ctsemOMX
          # Driver, C. C., Oud, J. H. L., & Voelkle, M. C. (2017).
          # Continuous time structural equation modelling with R package ctsem.
          # Journal of Statistical Software, 77(5), 1–36. https://doi.org/10.18637/jss.v077.i05

          start_i <- -.45
          names(start_i) <- ctMatrices$DRIFT[i,j]

          start <- c(start,
                     start_i
          )
        }else{
          # we use the same starting values as ctsemOMX
          # Driver, C. C., Oud, J. H. L., & Voelkle, M. C. (2017).
          # Continuous time structural equation modelling with R package ctsem.
          # Journal of Statistical Software, 77(5), 1–36. https://doi.org/10.18637/jss.v077.i05
          start_i <- -.05
          names(start_i) <- ctMatrices$DRIFT[i,j]

          start <- c(start,
                     start_i
          )
        }
      }

      if(grepl(pattern = "^[a-zA-Z]", x = ctMatrices$logDiagDIFFUSION[i,j])){
        parameters <- c(parameters,
                        ctMatrices$logDiagDIFFUSION[i,j])
        if(i == j){
          start_i <- .2
          names(start_i) <- ctMatrices$logDiagDIFFUSION[i,j]

          start <- c(start,
                     start_i
          )
        }else{
          start_i <- .1
          names(start_i) <- ctMatrices$logDiagDIFFUSION[i,j]

          start <- c(start,
                     start_i
          )
        }
      }
    }

  }

  transformations <- c("\ndouble tmpvalue = 0.0;\nbool DRIFTChanged = false;\narma::mat DRIFT = transformationList[\"DRIFT\"];\n")
  for(i in 1:nrow(ctMatrices$DRIFT)){
    for(j in 1:ncol(ctMatrices$DRIFT)){
      transformations <- c(transformations,
                           paste0(
                             "tmpvalue = DRIFT(",i-1,",", j-1,");\n",
                             "if(tmpvalue != ", ctMatrices$DRIFT[i,j],"){\n",
                             " DRIFTChanged = true;\n",
                             " DRIFT(",i-1,",", j-1,") = ", ctMatrices$DRIFT[i,j], ";"),
                           "\n}
                           "

      )
    }
  }
  transformations <- c(transformations,
                       paste0("if(DRIFTChanged) {transformationList[\"DRIFT\"] = DRIFT;}"),
                       paste0("arma::mat driftHash = transformationList[\"driftHash\"];\n",
                              "if(DRIFTChanged){\n",
                              " driftHash = kron(DRIFT, arma::eye(",nrow(ctMatrices$DRIFT), ",", nrow(ctMatrices$DRIFT),"))",
                              "+ kron(arma::eye(",nrow(ctMatrices$DRIFT),",",nrow(ctMatrices$DRIFT),"), DRIFT);\n",
                              "transformationList[\"driftHash\"] = driftHash;\n",
                              "}"))

  transformations <- c(transformations,
                       "bool logDiagDIFFUSIONChanged = false;\narma::mat logDiagDIFFUSION = transformationList[\"logDiagDIFFUSION\"];\narma::mat DIFFUSION = transformationList[\"DIFFUSION\"];")

  for(i in 1:nrow(ctMatrices$logDiagDIFFUSION)){
    for(j in 1:ncol(ctMatrices$logDiagDIFFUSION)){
      transformations <- c(transformations,
                           paste0("tmpvalue = logDiagDIFFUSION(",i-1,",", j-1,");\n",
                                  "if(tmpvalue != ", ctMatrices$logDiagDIFFUSION[i,j],"){\n",
                                  " logDiagDIFFUSIONChanged = true;\n",
                                  " logDiagDIFFUSION(",i-1,",", j-1,") = ", ctMatrices$logDiagDIFFUSION[i,j], ";"),
                           "\n}\n"
      )
    }
  }
  transformations <- c(transformations,
                       paste0("if(logDiagDIFFUSIONChanged){",
                              "transformationList[\"logDiagDIFFUSION\"] = logDiagDIFFUSION;",
                              "DIFFUSION = logDiagDIFFUSION;",
                              "DIFFUSION.diag() = arma::exp(logDiagDIFFUSION.diag());",
                              "transformationList[\"DIFFUSION\"] = DIFFUSION;",
                              "}"))

  # we will iterate over all unique time intervals and define
  # all transformations required to the drift and diffusion matrices
  # for these time intervals
  for(ti in unique(dataCTSEM$timetable$timeInterval)){
    if(ti == 0) next

    transformations <- c(transformations,
                         # create DRIFT
                         paste0("arma::mat ARCL_",which(unique(dataCTSEM$timetable$timeInterval) == ti),
                                " = transformationList[\"", "ARCL_",which(unique(dataCTSEM$timetable$timeInterval) == ti), "\"];\n",
                                "if(DRIFTChanged){\n",
                                " ARCL_",which(unique(dataCTSEM$timetable$timeInterval) == ti), " = ",
                                " arma::expmat(DRIFT*", ti, ");\n",
                                " transformationList[\"", "ARCL_",which(unique(dataCTSEM$timetable$timeInterval) == ti),
                                "\"] = ARCL_",which(unique(dataCTSEM$timetable$timeInterval) == ti),";\n}\n"
                         ),
                         # create diffusion
                         paste0(
                           "arma::mat LVCOV_",which(unique(dataCTSEM$timetable$timeInterval) == ti), " = transformationList[\"", "LVCOV_",which(unique(dataCTSEM$timetable$timeInterval) == ti), "\"];\n",
                           "if(DRIFTChanged | logDiagDIFFUSIONChanged){\n",
                           " LVCOV_",which(unique(dataCTSEM$timetable$timeInterval) == ti), " = ",
                           "arma::reshape(arma::inv(driftHash) * ",
                           "(arma::expmat(driftHash*", ti, ") - arma::eye(arma::size(arma::expmat(driftHash*", ti, "))))*",
                           "arma::vectorise(DIFFUSION),",nrow(ctMatrices$DRIFT),",",nrow(ctMatrices$DRIFT),");\n",
                           "transformationList[\"", "LVCOV_",which(unique(dataCTSEM$timetable$timeInterval) == ti), "\"] = LVCOV_",which(unique(dataCTSEM$timetable$timeInterval) == ti),";\n}\n
                          "
                         )
    )

    for(oc in which(dataCTSEM$timetable$timeInterval == ti)){

      arcl_u <- matrix(stringr::str_replace(string = arclModel$arcl,
                                            pattern = "\\(u\\)",
                                            replacement = paste0("u",oc)),
                       nrow = nrow(arclModel$arcl),
                       ncol = ncol(arclModel$arcl)
      )

      parameters <- c(parameters,
                      unique(c(arcl_u)))

      for(a1 in 1:nrow(arcl_u)){
        for(a2 in 1:ncol(arcl_u)){
          transformations <- c(transformations,
                               paste0(arcl_u[a1,a2], " = ARCL_",
                                      which(unique(dataCTSEM$timetable$timeInterval) == ti),
                                      "(", a1-1, ",", a2-1, ");"
                               )
          )
        }
      }

      cov_u <- matrix(stringr::str_replace(string = arclModel$covs,
                                           pattern = "\\(u\\)",
                                           replacement = paste0("u",oc)),
                      nrow = nrow(arclModel$covs),
                      ncol = ncol(arclModel$covs)
      )

      parameters <- c(parameters,
                      unique(c(cov_u)))

      for(c1 in 1:nrow(cov_u)){
        for(c2 in 1:ncol(cov_u)){
          transformations <- c(transformations,
                               paste0(cov_u[c1,c2], " = ",
                                      # lessSEM internally uses the exponential of
                                      # variance parameters to avoid negative variances.
                                      # We therefore take the log in case of variances:
                                      ifelse(c1==c2, "log(",""),
                                      "LVCOV_",
                                      which(unique(dataCTSEM$timetable$timeInterval) == ti),
                                      "(", c1-1, ",", c2-1, ")",
                                      ifelse(c1==c2, ")",""), ";"
                               )
          )
        }
      }

    }
  }

  # combine everything:
  parameters <- unique(parameters)
  start <- start[unique(names(start))]

  parametersStr <- paste0("parameters: ", paste0(parameters, collapse = ", "))
  startStr <- paste0("start: ", paste0(paste0(names(start), " = ", start), collapse = ", "))
  transformationsStr <- paste0(transformations, collapse = "\n")

  transform <- paste0(
    c(parametersStr,
      startStr,
      transformationsStr
    ),
    collapse = "\n"
  )

  return(list(transformation = transform,
              transformationList = transformationList))
}
