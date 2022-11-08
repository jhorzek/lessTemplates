CTSEM <- function(model,
                  data){

  cat("\nSetting up a continuous time structural equation model.\n")

  if(!all(c("person", "time") %in% colnames(data)))
    stop("Could not find the columns person and time in the data set.")
  if(!is(object = data, class2 = "data.frame"))
    data <- as.data.frame(data)

  dataCTSEM <- lessTransformations:::.prepareDataCTSEM(data = data)

  arclData <- dataCTSEM$data

  # remove unnecessary white space
  syntax <- lessSEM:::.reduceSyntax(syntax = model)
  syntax <- lessTransformations:::.removeWhitespace(syntax = syntax)
  syntax <- lessTransformations:::.makeSingleLine(syntax = syntax)

  # check if all dynamics are time dependent
  lessTransformations:::.checkDynamics(syntax = syntax)

  # find statements regarding Wiener process
  wiener <- lessTransformations:::.getWienerProcessCTSEM(syntax)

  # find the names of all variables
  variableNames <- lessTransformations:::.getVariableNamesCTSEM(syntax = syntax)

  latents <- list(
    occasionDependent = variableNames$occasionDependent[!variableNames$occasionDependent %in% colnames(data)],
    fixed = variableNames$fixed[!variableNames$fixed %in% colnames(data)]
  )

  manifests <- list(
    occasionDependent = variableNames$occasionDependent[variableNames$occasionDependent %in% colnames(data)],
    fixed = variableNames$fixed[variableNames$fixed %in% colnames(data)]
  )

  # set up discrete time model as basis
  arclModel <- lessTransformations:::.getDiscreteBasis(syntax = syntax,
                                                       latents = latents,
                                                       wiener = wiener)

  ctMatrices <- lessTransformations:::.getCTSEMMatrices(syntax = syntax,
                                                        arclModel = arclModel)

  clpm <- lessTransformations::CLPM(model = arclModel$syntax,
                                    data = arclData[,colnames(arclData) != "time"],
                                    silent = TRUE)

  transformations <- lessTransformations:::.getCTSEMTransformations(arclModel = arclModel,
                                                                    dataCTSEM = dataCTSEM,
                                                                    ctMatrices = ctMatrices)
  return(list(
    model = clpm$model,
    data = clpm$data,
    transformation = transformations
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


.getWienerProcessCTSEM <- function(syntax){

  wiener <- unlist(stringr::str_extract_all(string = syntax,
                                            pattern = "d[0-9]*_W[0-9]*\\(t\\)"))
  if(length(wiener) == 0)
    stop("Could not find a statement regarding the wiener process (e.g., d_W1(t)).")

  return(wiener)
}


#' .getVariableNamesCLPM
#'
#' extracts the names of the variables from the syntax
#' @param syntax string
#' @return list with names of variables
#' @keywords internal
.getVariableNamesCTSEM <- function(syntax){

  # remove all parameters
  syntax_t <- gsub(pattern = "[a-zA-Z0-9_]+\\*",
                   replacement = "",
                   x = syntax)
  # remove means
  syntax_t <- gsub(pattern = "~1",
                   replacement = "",
                   x = syntax_t)
  # split at operators
  variableNames <- unlist(stringr::str_split(string = syntax_t,
                                             pattern = "\\+|=~|~~|~"))
  # remove Wiener process
  variableNames <- variableNames[!grepl(pattern = "d[0-9]*_W[0-9]*\\(t\\)",
                                        x = variableNames)]
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

.getDiscreteBasis <- function(syntax, latents, wiener){

  # keep all measurement equations
  measurementEquations <- syntax[grepl(pattern = "=~", x = syntax)]

  # replace dynamics in measurement equations
  measurementEquations <- stringr::str_replace_all(string = measurementEquations,
                                                   pattern = "\\(t\\)",
                                                   replacement = "_\\(u\\)")

  # find the highest order:
  highestOrder <- max(lessTransformations:::.getHighestOrderCTSEM(latents = latents,
                                                                  wiener = wiener))

  dynamicLatents <- latents$occasionDependent

  # remove all differential statements
  dynamicLatents <- unique(stringr::str_remove(string = dynamicLatents,
                                               pattern = "^d[0-9]*_"))

  # create names used in arcl model for dynamics:
  if(highestOrder == 1){
    arclDynamics <- dynamicLatents
  }else{
    arclDynamics <- c(dynamicLatents,
                      paste0("d",1:(highestOrder-1), "_", dynamicLatents)
    )
  }

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

  # add remaining covariances
  measurementCovs <- syntax[grepl(pattern = "~~", x = syntax)]

  # replace dynamics in measurement equations
  measurementCovs <- stringr::str_replace_all(string = measurementCovs,
                                              pattern = "\\(t\\)",
                                              replacement = "_\\(u\\)")
  # add
  arclSyntax <- c(arclSyntax,
                  measurementCovs)

  # add manifest intercepts
  measurementInts <- syntax[grepl(pattern = "\\*1", x = syntax)]
  measurementInts <- stringr::str_replace_all(string = measurementInts,
                                              pattern = "\\(t\\)",
                                              replacement = "_\\(u\\)")
  arclSyntax <- c(arclSyntax,
                  measurementInts)

  arclSyntax <- paste0(arclSyntax,
                       collapse = "\n")

  return(list(syntax = arclSyntax,
              arcl = arcls,
              covs = covs,
              dynamicLatents = dynamicLatents,
              highestOrder = highestOrder
  ))
}

.getHighestOrderCTSEM <- function(latents, wiener){
  ordersDynamics <- stringr::str_extract(string = c(latents$occasionDependent),
                                         pattern = "^d[0-9]*")
  ordersWiener <- stringr::str_extract(string = c(wiener),
                                       pattern = "^d[0-9]*")
  if(all(c(is.na(ordersDynamics), is.na(ordersWiener))))
    stop("Could not find any dynamics (e.g., statements using d_eta(t)).")
  ordersDynamics <- ordersDynamics[!is.na(ordersDynamics)]
  ordersWiener <- ordersWiener[!is.na(ordersWiener)]

  getOrder <- function(orders){
    if(all(orders == "d")){
      highestOrder <- 1
    }else{
      orders <- orders[orders!="d"]
      highestOrder <- max(as.integer(stringr::str_extract(string = orders,
                                                          pattern = "[0-9]+$")))
    }
    return(highestOrder)
  }

  highestOrders <- c("dynamics" = getOrder(ordersDynamics),
                     "wiener" = getOrder(ordersWiener))

  return(highestOrders)
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

  if(arclModel$highestOrder == 1){
    dependent <- paste0("d_", arclModel$dynamicLatents)
    predictor <- arclModel$dynamicLatents
  }else{
    dependent <- paste0(rep(c("d", paste0("d",2:arclModel$highestOrder)), each = length(arclModel$dynamicLatents)),
                        "_",
                        rep(arclModel$dynamicLatents, arclModel$highestOrder))
    predictor <- paste0(rep(c("", paste0("d",1:(arclModel$highestOrder-1))), each = length(arclModel$dynamicLatents)),
                        "_",
                        rep(arclModel$dynamicLatents, arclModel$highestOrder))
  }

  # Set up the matrices; the individual values will be replaced later on
  DRIFT <- matrix(paste0("drift_", rep(dependent, length(predictor)), "_", rep(predictor, each = length(predictor))),
                  nrow = length(predictor),
                  ncol = length(predictor),
                  dimnames = list(dependent, predictor)
  )

  DIFFUSION <- matrix(paste0("diffusion_", rep(predictor, length(predictor)),
                             "_", rep(predictor, each = length(predictor))),
                      nrow = length(predictor),
                      ncol = length(predictor),
                      dimnames = list(predictor, predictor)
  )

  for(i in 1:nrow(DIFFUSION)){
    for(j in i:nrow(DIFFUSION)){
      DIFFUSION[i,j] <- DIFFUSION[j,i]
    }
  }

  if(arclModel$highestOrder == 1){

    return(
      list(DRIFT = DRIFT,
           DIFFUSION = DIFFUSION)
    )

  }else{

    stop("Higher order models not yet implemented")

  }
}

.getCTSEMTransformations <- function(arclModel, dataCTSEM, ctMatrices){

  parameters <- c() # names of parameters
  start <- c() # starting values
  for(i in 1:nrow(ctMatrices$DRIFT)){
    for(j in 1:ncol(ctMatrices$DRIFT)){

      if(grepl(pattern = "^[a-zA-Z]", x = ctMatrices$DRIFT[i,j])){
        parameters <- c(parameters,
                        ctMatrices$DRIFT[i,j])
        if(i == j){
          start_i <- -.1
          names(start_i) <- ctMatrices$DRIFT[i,j]

          start <- c(start,
                     start_i
          )
        }else{
          start_i <- 0
          names(start_i) <- ctMatrices$DRIFT[i,j]

          start <- c(start,
                     start_i
          )
        }
      }

      if(grepl(pattern = "^[a-zA-Z]", x = ctMatrices$DIFFUSION[i,j])){
        parameters <- c(parameters,
                        ctMatrices$DIFFUSION[i,j])
        if(i == j){
          start_i <- .2
          names(start_i) <- ctMatrices$DIFFUSION[i,j]

          start <- c(start,
                     start_i
          )
        }else{
          start_i <- 0
          names(start_i) <- ctMatrices$DIFFUSION[i,j]

          start <- c(start,
                     start_i
          )
        }
      }
    }

  }

  transformations <- paste0("arma::mat DRIFT(", nrow(ctMatrices$DRIFT), ",", ncol(ctMatrices$DRIFT), ")")
  for(i in 1:nrow(ctMatrices$DRIFT)){
    for(j in 1:ncol(ctMatrices$DRIFT)){
      transformations <- c(transformations,
                           paste0("DRIFT(",i-1,",", j-1,") = ", ctMatrices$DRIFT[i,j])
      )
    }
  }
  transformations <- c(transformations,
                       paste0("arma::mat driftHash = kron(DRIFT, arma::eye(",nrow(ctMatrices$DRIFT), ",", nrow(ctMatrices$DRIFT),"))",
                              "+ kron(arma::eye(",nrow(ctMatrices$DRIFT),",",nrow(ctMatrices$DRIFT),"), DRIFT)")
  )

  transformations <- c(transformations,
                       paste0("arma::mat DIFFUSION(", nrow(ctMatrices$DIFFUSION), ",", ncol(ctMatrices$DIFFUSION), ")")
  )

  for(i in 1:nrow(ctMatrices$DIFFUSION)){
    for(j in 1:ncol(ctMatrices$DIFFUSION)){
      transformations <- c(transformations,
                           paste0("DIFFUSION(",i-1,",", j-1,") = ", ctMatrices$DIFFUSION[i,j])
      )
    }
  }

  # we will iterate over all unique time intervals and define
  # all transformations required to the drift and diffusion matrices
  # for these time intervals
  for(ti in unique(dataCTSEM$timetable$timeInterval)){
    if(ti == 0) next

    transformations <- c(transformations,
                         # create DRIFT
                         paste0(
                           "arma::mat ARCL_",which(unique(dataCTSEM$timetable$timeInterval) == ti), " = ",
                           " arma::expmat(DRIFT*", ti, ")"
                         ),
                         # create diffusion
                         paste0(
                           "arma::mat LVCOV_",which(unique(dataCTSEM$timetable$timeInterval) == ti), " = ",
                           "arma::reshape(arma::inv(driftHash) * ",
                           "(arma::expmat(driftHash*", ti, ") - arma::eye(arma::size(arma::expmat(driftHash*", ti, "))))*",
                           "arma::vectorise(DIFFUSION),",nrow(ctMatrices$DRIFT),",",nrow(ctMatrices$DRIFT),")
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
                                      "(", a1-1, ",", a2-1, ")"
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
                               paste0(cov_u[c1,c2], " = LVCOV_",
                                      which(unique(dataCTSEM$timetable$timeInterval) == ti),
                                      "(", c1-1, ",", c2-1, ")"
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
      transformationsStr),
    collapse = "\n"
  )

  return(transform)
}
