#' parseARCL
#'
#' @param modelSyntax model syntax similar to lavaan
#' @returns list with parameter table, names of manifest and names of latent variables
#' @keywords internal
parseARCL <- function(modelSyntax, occasions){

  modelSyntax <- reduceSyntax(modelSyntax)

  parameterTable <- makeParameterTable(modelSyntax)

  latentVariables <- parameterTable$lhs[parameterTable$operator == "=~"]
  latentVariables <- sub("\\(.*", "", latentVariables)
  latentVariables <- unique(latentVariables)

  manifestVariables <- parameterTable$rhs[parameterTable$operator == "=~"]
  manifestVariables <- sub("\\(.*", "", manifestVariables)
  manifestVariables <- unique(manifestVariables)

  # let's check if there are regressions on manifest variables
  if(any((parameterTable$lhs %in% manifestVariables) & parameterTable$operator == "~"))
    stop("Regressions on manifest variables (e.g., y1~y2+y3) are currently not allowed.")

  # let's check if there are left hand side variables which have the wrong time index
  if(any(!(grepl("(u)",parameterTable$lhs) | grepl("\\(1",parameterTable$lhs)) & parameterTable$operator != "~~"))
    stop("Some of your variables on the left hand side of the equations are not specified with (u).")

  # for higher order lags:
  parameterTable <- addTimeInformation(parameterTable)

  return(list(
    latentVariables = latentVariables,
    manifestVariables  = manifestVariables,
    parameterTable = parameterTable
  ))
}

#' reduceSyntax
#'
#' @param modelSyntax model syntax similar to lavaan
#' @returns reduced syntax
#' @keywords internal
reduceSyntax <- function(modelSyntax){

  # first, split rows and remove everything we don't need
  # split rows
  modelSyntax <- stringr::str_split(string = modelSyntax,
                                    pattern = "\\n")[[1]]
  # remove white space
  modelSyntax <- stringr::str_replace_all(string = modelSyntax,
                                          pattern = "\\s",
                                          replacement = "")
  # remove comments
  hasComment <- stringr::str_locate(modelSyntax,
                                    "#|!")
  for(i in 1:length(modelSyntax)){
    if(!is.na(hasComment[i,1])){
      if(hasComment[i,1] == 1){
        modelSyntax[i] <- ""
        next
      }
      modelSyntax[i] <- stringr::str_trunc(
        modelSyntax[i], width = hasComment[i,1]-1, ellipsis = ""
      )
    }
  }

  # remove empty
  modelSyntax <- modelSyntax[modelSyntax != ""]

  return(modelSyntax)
}


#' makeParameterTable
#'
#' create a parameter table
#' @param modelSyntax model syntax
#' @returns data.frame with parameters
makeParameterTable <- function(modelSyntax){

  parameterTable <- c()

  operators <- c(
    "=~", # loadings
    "~~", # covariances
    "~" # regressions
  )

  for(ln in 1:length(modelSyntax)){
    for(op in operators){
      # first, get the names of all latent variables
      hasOperator <- grepl(op, modelSyntax[ln])

      if(!any(hasOperator)) next

      leftRight <- stringr::str_split(modelSyntax[ln],
                                      pattern = op
      )
      lhs <- unlist(lapply(leftRight,
                           function(x) return(x[1])))

      # now, let's check the right hand side:
      rhs <- unlist(lapply(leftRight,
                           function(x) return(x[2])))
      # elements are separated by +
      rhsElements <- stringr::str_split(string = rhs,
                                        pattern = "\\+")

      for(i in 1:length(lhs)){

        # modifiers can be labels for parameters or they can
        # be fixed parameters. They are always indicated by "*"
        hasModifier <- grepl("\\*", rhsElements[[i]])

        for(j in 1:length(hasModifier)){

          if(!hasModifier[j]){
            parameterTable <- rbind(parameterTable,
                                    data.frame(
                                      lhs = lhs[i],
                                      rhs = rhsElements[[i]][j],
                                      operator = op,
                                      label = paste0("label(",
                                                     lhs[i],
                                                     op,
                                                     rhsElements[[i]][j],
                                                     ")"),
                                      value = NA,
                                      free = TRUE,
                                      lhs_tMinus = NA,
                                      rhs_tMinus = NA,
                                      timeDependent = FALSE,
                                      location = NA,
                                      row = NA,
                                      col = NA
                                    )
            )
          }else{

            elements <- stringr::str_split(string = rhsElements[[i]][j],
                                           pattern = "\\*")[[1]]
            if(stringr::str_detect(string = elements[1], pattern = "[:letter:]")){
              # if there is a letter in the left hand side, we assume that it is a label

              # Now, we also check if the parameter differs by occasion

              isTimeDependent <- grepl(pattern = "_(u)",
                                       x = elements[1],
                                       fixed = TRUE)

              if(isTimeDependent){
                parLabel <- stringr::str_split(string = elements[1],
                                               pattern = "_")[[1]][1]
                parameterTable <- rbind(parameterTable,
                                        data.frame(
                                          lhs = lhs[i],
                                          rhs = elements[2],
                                          operator = op,
                                          label = parLabel,
                                          value = NA,
                                          free = TRUE,
                                          lhs_tMinus = NA,
                                          rhs_tMinus = NA,
                                          location = NA,
                                          timeDependent = TRUE,
                                          row = NA,
                                          col = NA
                                        )
                )
              }else{
                parameterTable <- rbind(parameterTable,
                                        data.frame(
                                          lhs = lhs[i],
                                          rhs = elements[2],
                                          operator = op,
                                          label = elements[1],
                                          value = NA,
                                          free = TRUE,
                                          lhs_tMinus = NA,
                                          rhs_tMinus = NA,
                                          location = NA,
                                          timeDependent = FALSE,
                                          row = NA,
                                          col = NA
                                        )
                )
              }

            }else if(stringr::str_count(string = elements[1],
                                        pattern = "[:digit:]")){
              # if it's all numbers, we assume that it is a fixed value
              parameterTable <- rbind(parameterTable,
                                      data.frame(
                                        lhs = lhs[i],
                                        rhs = elements[2],
                                        operator = op,
                                        label = paste0("label(",
                                                       lhs[i], op, rhsElements[[i]][j],
                                                       ")"),
                                        value = as.numeric(elements[1]),
                                        free = FALSE,
                                        lhs_tMinus = NA,
                                        rhs_tMinus = NA,
                                        location = NA,
                                        timeDependent = FALSE,
                                        row = NA,
                                        col = NA
                                      )

              )

            }else{
              stop(paste0("Unknown modifier ", elements[1]))
            }

          }

        }
      }
      # if operator was found: no need to iterate over the other ones
      break
    }# end operators
  }# end lines

  return(parameterTable)
}

#' addTimeInformation
#'
#' checks for higher order lags
#' @param parameterTable parameter table
#' @return updated parameter table
#' @keywords internal
addTimeInformation <- function(parameterTable){

  # time operators are always embraced in [] and specified with
  # u-k
  # The maximal time operator will always be on u and should appear on
  # the left hand side
  elements <- stringr::str_extract(string = parameterTable$lhs,
                                   pattern = "(?<=\\().*(?=\\))")
  if(all(is.na(elements)))
    stop("The model does not have a variable for time u (e.g., eta(u)) specified on the left hand side of an equation.")
  elements[elements == "u"] <- "u-0"
  elements <- as.integer(stringr::str_remove_all(string = elements,
                                                 pattern = "u-")
  )
  parameterTable$lhs_tMinus <- elements

  if(!any(parameterTable$lhs_tMinus == 0))
    stop("The model does not have a variable for time u (e.g., eta(u)) specified on the left hand side of an equation.")

  # the minimal time operator will be on the right hand side
  # Let's see what it is:
  elements <- stringr::str_extract(string = parameterTable$rhs,
                                   pattern = "(?<=\\().*(?=\\))")

  elements[parameterTable$rhs == "1"] <- "u"

  if(all(is.na(elements)))
    stop("The model does not have a variable for time u-k (e.g., eta(u-1)) specified on the right hand side of an equation.")
  elements[elements == "u"] <- "u-0"
  elements <- as.integer(stringr::str_remove_all(string = elements,
                                                 pattern = "u-")
  )
  parameterTable$rhs_tMinus <- elements

  if(any(
    parameterTable$lhs == parameterTable$rhs &
    parameterTable$operator == "~")
  ) stop("Self-loops are currently not allowed (e.g., eta(u) ~ eta(u).")

  if(any(
    parameterTable$lhs_tMinus == parameterTable$rhs_tMinus &
    parameterTable$operator == "~" &
    parameterTable$rhs != "1"
  )
  ) stop("Contemporaneous effects (e.g., eta(u) ~ eta2(u)) are currently not allowed.")

  return(parameterTable)

}

initializeMatrices <- function(parserResult){

  # first, let's find out what the maximal number of time lags is
  maxLag <- max(parserResult$parameterTable$rhs_tMinus,
                parserResult$parameterTable$lhs_tMinus,
                na.rm = TRUE)

  # check if maxLag is in (co) variances; if so, we must add another lag
  if(any(parserResult$parameterTable$operator[c(which(parserResult$parameterTable$rhs_tMinus == maxLag),
                                                which(parserResult$parameterTable$lhs_tMinus == maxLag)
  )] == "~~"))
    maxLag <- maxLag + 1


  nLatent <- length(parserResult$latentVariables)
  nManifest <- length(parserResult$manifestVariables)

  # Starting with the latent equations.
  ## Autoregressive and cross-lagged effects
  ARCL <- setARCL(parserResult = parserResult,
                  nLatent = nLatent,
                  maxLag = maxLag)
  parserResult <- ARCL$parserResult
  ARCL <- ARCL$ARCL

  # latent means
  latentMeans <- setLatentMeans(ARCL, parserResult)
  parserResult <- latentMeans$parserResult
  latentMeans <- latentMeans$latentMeans

  ## residual covariances
  latentVariableCov <- setLatentVariableCov(ARCL = ARCL,
                                            parserResult = parserResult)
  parserResult <- latentVariableCov$parserResult
  latentVariableCov <- latentVariableCov$latentVariableCov

  # Manifest Equations
  Loadings <- setLoadings(parserResult = parserResult,
                          ARCL = ARCL,
                          nLatent = nLatent,
                          nManifest = nManifest,
                          maxLag = maxLag)
  parserResult <- Loadings$parserResult
  Loadings <- Loadings$Loadings

  manifestMeans <- setManifestMeans(Loadings, parserResult)
  parserResult <- manifestMeans$parserResult
  manifestMeans <- manifestMeans$manifestMeans

  manifestVariableCov <- setManifestVariableCov(Loadings = Loadings,
                                                parserResult = parserResult)
  parserResult <- manifestVariableCov$parserResult
  manifestVariableCov <- manifestVariableCov$manifestVariableCov

  # check identification
  for(lv in paste0(parserResult$latentVariables, "(u)")){
    if(any(!parserResult$parameterTable$free[
      (parserResult$parameterTable$lhs == lv) & (parserResult$parameterTable$operator == "=~")
    ])
    ) next
    warning(paste0(lv, " may not be identified. Check if one of the loadings was fixed to 1."))
  }


  return(list(
    manifestVariables = parserResult$manifestVariables,
    latentVariables = parserResult$latentVariables,
    parameterTable = parserResult$parameterTable,
    model = list(
      ARCL = ARCL,
      latentMeans = latentMeans,
      latentVariableCov = latentVariableCov,
      Loadings = Loadings,
      manifestMeans = manifestMeans,
      manifestVariableCov = manifestVariableCov
    )
  ))
}

setARCL <- function(parserResult, nLatent, maxLag ){
  # We need an nLatent*nLags x nLatent*nLags matrix
  ARCL <- matrix("0",
                 nrow = nLatent*(maxLag),
                 ncol = nLatent*(maxLag))
  rownames(ARCL) <- paste0(rep(parserResult$latentVariables, each = maxLag),
                           rep(paste0("(u-", 0:(maxLag-1)), nLatent), ")")

  for(v in parserResult$latentVariables)  rownames(ARCL)[rownames(ARCL) == paste0(v, "(u-0)")] <- paste0(
    v, "(u)")

  colnames(ARCL) <- paste0(rep(parserResult$latentVariables, each = maxLag),
                           rep(paste0("(u-", seq_len(maxLag)), nLatent), ")")

  for(ro in 1:nrow(parserResult$parameterTable)){
    # insert values
    if(
      (parserResult$parameterTable$operator[ro] == "~") &&
      (parserResult$parameterTable$lhs[ro] %in% rownames(ARCL)) &&
      (parserResult$parameterTable$rhs[ro] %in% colnames(ARCL))
    ) {
      parserResult$parameterTable$location[ro] <- "ARCL"
      parserResult$parameterTable$row[ro] <- which(rownames(ARCL) == parserResult$parameterTable$lhs[ro])
      parserResult$parameterTable$col[ro] <- which(colnames(ARCL) == parserResult$parameterTable$rhs[ro])

      if(!parserResult$parameterTable$free[ro]){
        ARCL[
          parserResult$parameterTable$lhs[ro],
          parserResult$parameterTable$rhs[ro]
        ] <- parserResult$parameterTable$value[ro]
      }else{
        ARCL[
          parserResult$parameterTable$lhs[ro],
          parserResult$parameterTable$rhs[ro]
        ] <-
          parserResult$parameterTable$label[ro]
      }
    }
  }

  # for higher order models:
  for(ro in rownames(ARCL)){
    if(any(colnames(ARCL) == ro))  ARCL[ro,ro] <- "1"
  }

  return(list(ARCL = ARCL,
              parserResult = parserResult))
}

setLatentMeans <- function(ARCL, parserResult){

  latentMeans <- matrix(0, nrow = nrow(ARCL),
                        dimnames = list(rownames(ARCL),
                                        "1"))
  for(ro in 1:nrow(parserResult$parameterTable)){

    if(parserResult$parameterTable$rhs[ro] != "1") next
    if(!parserResult$parameterTable$lhs[ro] %in% paste0(parserResult$latentVariables, "(u)")) next
    if(parserResult$parameterTable$free[ro]){
      latentMeans[parserResult$parameterTable$lhs[ro],
                  parserResult$parameterTable$rhs[ro]] <- parserResult$parameterTable$label[ro]
    }else{
      latentMeans[parserResult$parameterTable$lhs[ro],
                  parserResult$parameterTable$rhs[ro]] <- parserResult$parameterTable$value[ro]
    }
    parserResult$parameterTable$location[ro] <- "latentMeans"
    parserResult$parameterTable$row[ro] <- which(rownames(latentMeans) == parserResult$parameterTable$lhs[ro])
    parserResult$parameterTable$col[ro] <- 1
  }

  return(list(
    "latentMeans" = latentMeans,
    "parserResult" = parserResult
  ))

}

setLatentVariableCov <- function(ARCL, parserResult){

  latentVariableCov <- matrix(0,
                              nrow = nrow(ARCL),
                              ncol = ncol(ARCL),
                              dimnames = list(rownames(ARCL),
                                              rownames(ARCL)))

  for(ro in 1:nrow(parserResult$parameterTable)){
    # insert values
    if (
      (parserResult$parameterTable$operator[ro] == "~~") &&
      (parserResult$parameterTable$lhs[ro] %in% rownames(latentVariableCov)) &&
      (parserResult$parameterTable$rhs[ro] %in% colnames(latentVariableCov))
    ) {

      parserResult$parameterTable$location[ro] <- "latentVariableCov"
      parserResult$parameterTable$row[ro] <- which(rownames(latentVariableCov) == parserResult$parameterTable$lhs[ro])
      parserResult$parameterTable$col[ro] <- which(colnames(latentVariableCov) == parserResult$parameterTable$rhs[ro])

      if(parserResult$parameterTable$free[ro]){
        latentVariableCov[
          parserResult$parameterTable$lhs[ro],
          parserResult$parameterTable$rhs[ro]
        ] <-
          parserResult$parameterTable$label[ro]
      }else{
        latentVariableCov[
          parserResult$parameterTable$lhs[ro],
          parserResult$parameterTable$rhs[ro]
        ] <-
          parserResult$parameterTable$value[ro]
      }
    }
  }

  # check if all variances of the first order are in the parameter table
  lvAtT <- paste0(parserResult$latentVariables, "(u)")
  if(anyNA(diag(latentVariableCov[lvAtT,lvAtT,drop = FALSE]))){

    # add those to the parameter list
    for(v in which(is.na(diag(latentVariableCov[lvAtT,lvAtT,drop = FALSE])))){
      latentVariableCov[v,v] <- .1

      parserResult$parameterTable <- rbind(
        parserResult$parameterTable,
        data.frame(
          lhs = rownames(latentVariableCov)[v],
          rhs = colnames(latentVariableCov)[v],
          operator = "~~",
          label = paste0("label(", rownames(latentVariableCov)[v], "~~", colnames(latentVariableCov)[v], ")"),
          value = .1,
          free = TRUE,
          lhs_tMinus = 0,
          rhs_tMinus = 0,
          location = "latentVariableCov",
          timeDependent = FALSE,
          row = v,
          col = v
        )
      )
    }

  }

  return(list(latentVariableCov = latentVariableCov,
              parserResult = parserResult))
}


setLoadings <- function(parserResult, ARCL, nLatent, nManifest, maxLag){
  # We need an nManifest x nLatent*maxLag matrix
  Loadings <- matrix(0,
                     nrow = nManifest,
                     ncol = nLatent*maxLag)
  rownames(Loadings) <- paste0(parserResult$manifestVariables, "(u)")

  colnames(Loadings) <- rownames(ARCL)

  for(ro in 1:nrow(parserResult$parameterTable)){
    if(parserResult$parameterTable$operator[ro] != "=~") next
    from <- parserResult$parameterTable$lhs[ro]
    to <- parserResult$parameterTable$rhs[ro]

    if(!from %in% colnames(Loadings))
      stop(paste0("Could not set ", from , " in Loadings."))
    if(!to %in% rownames(Loadings))
      stop(paste0("Could not set ", to , " in Loadings."))

    if(parserResult$parameterTable$free[ro]){
      Loadings[to,from] <- parserResult$parameterTable$label[ro]
    }else{
      Loadings[to,from] <- parserResult$parameterTable$value[ro]
    }
    parserResult$parameterTable$location[ro] <- "Loadings"
    parserResult$parameterTable$row[ro] <- which(rownames(Loadings) == parserResult$parameterTable$rhs[ro])
    parserResult$parameterTable$col[ro] <- which(colnames(Loadings) == parserResult$parameterTable$lhs[ro])
  }

  return(list(Loadings = Loadings,
              parserResult = parserResult))

}

setManifestMeans <- function(Loadings, parserResult){

  manifestMeans <- matrix(0, nrow = nrow(Loadings),
                          dimnames = list(rownames(Loadings),
                                          "1"))
  for(ro in 1:nrow(parserResult$parameterTable)){

    if(parserResult$parameterTable$rhs[ro] != "1") next
    if(!parserResult$parameterTable$lhs[ro] %in% paste0(parserResult$manifestVariables, "(u)")) next
    if(parserResult$parameterTable$free[ro]){
      manifestMeans[parserResult$parameterTable$lhs[ro],
                    parserResult$parameterTable$rhs[ro]] <- parserResult$parameterTable$label[ro]
    }else{
      manifestMeans[parserResult$parameterTable$lhs[ro],
                    parserResult$parameterTable$rhs[ro]] <- parserResult$parameterTable$value[ro]
      }
    parserResult$parameterTable$location[ro] <- "manifestMeans"
    parserResult$parameterTable$row[ro] <- which(rownames(manifestMeans) == parserResult$parameterTable$lhs[ro])
    parserResult$parameterTable$col[ro] <- 1
  }

  if(!any(parserResult$parameterTable$rhs == "1")){
    # add manifest means
    parserResult$parameterTable <- rbind(
      parserResult$parameterTable,
      data.frame(
        lhs = rownames(manifestMeans),
        rhs = "1",
        operator = "~",
        label = paste0("label(",rownames(manifestMeans), "~1)"),
        value = 0,
        free = TRUE,
        lhs_tMinus = 0,
        rhs_tMinus = 0,
        location = "manifestMeans",
        person = 0,
        time = 0,
        row = 1:nrow(manifestMeans),
        col = 1
      )
    )

  }

  return(list(
    "manifestMeans" = manifestMeans,
    "parserResult" = parserResult
  ))

}

setManifestVariableCov <- function(Loadings, parserResult){

  manifestVariableCov <- matrix(0,
                                nrow = nrow(Loadings),
                                ncol = nrow(Loadings),
                                dimnames = list(rownames(Loadings),
                                                rownames(Loadings)))
  diag(manifestVariableCov) <- NA

  for(ro in 1:nrow(parserResult$parameterTable)){
    # insert values
    if (
      (parserResult$parameterTable$operator[ro] == "~~") &&
      (parserResult$parameterTable$lhs[ro] %in% rownames(manifestVariableCov)) &&
      (parserResult$parameterTable$rhs[ro] %in% colnames(manifestVariableCov))
    ) {

      parserResult$parameterTable$location[ro] <- "manifestVariableCov"
      parserResult$parameterTable$row[ro] <- which(rownames(manifestVariableCov) == parserResult$parameterTable$lhs[ro])
      parserResult$parameterTable$col[ro] <- which(colnames(manifestVariableCov) == parserResult$parameterTable$rhs[ro])

      if(parserResult$parameterTable$free[ro]) {
        manifestVariableCov[
        parserResult$parameterTable$lhs[ro],
        parserResult$parameterTable$rhs[ro]
      ] <-
        parserResult$parameterTable$label[ro]
      }else{
        manifestVariableCov[
          parserResult$parameterTable$lhs[ro],
          parserResult$parameterTable$rhs[ro]
        ] <-
          parserResult$parameterTable$value[ro]
      }
    }
  }

  return(list(manifestVariableCov = manifestVariableCov,
              parserResult = parserResult))

}

#' makeARCL
#'
#' Creates and autoregressive cross-lagged panel model for use with lavaan or lessSEM.
#'
#' @param data data set in long format. The following columns are required: (1) "id"
#' - a column which indicates which person a specific row of the data set belongs to.
#' (2) "occasion" - a column indicating which measurement occasion a row belongs to.
#' All other columns should contain the data required for the model. Optionally, a
#' column called "lag" can be added which specifies if parameters for this lag should
#' be estimated lag-specific (see below)
#' @param latentIntercepts nLatent x 1 matrix with latent intercepts. If a parameter
#' is fixed, set it to the numeric value (e.g., "0"). If a parameter is estimated, give it
#' a name (e.g., "latentIntercept1").
#' @param ARCL nLatent x nLatent matrix with autoregressive and cross-lagged parameters
#' @param latentResiduals nLatent x nLatent matrix with covariances of latent residuals
#' @param loadings an nManifest x nLatent matrix specifying the loadings.
#' @param manifestIntercepts nManifest x 1 matrix with manifest intercepts
#' @param measurementError nManifest x nManifest matrix with measurement errors
#' @param lagSpecific vector indicating which of the matrices above differs by lag.
#' Requires the specification of the lag column.
#' @return list with lavaan model syntax and data set for lavaan
#' @keywords internal
makeARCL <- function(data,
                     latentIntercepts,
                     ARCL,
                     latentResiduals,
                     loadings,
                     manifestIntercepts,
                     measurementError,
                     lagSpecific = NULL){

  isLagSpecific <- data.frame(latentIntercepts = FALSE,
                              ARCL = FALSE,
                              latentResiduals = FALSE,
                              loadings = FALSE,
                              manifestIntercepts = FALSE,
                              measurementError = FALSE
  )

  if(!is.null(lagSpecific)){
    if(! "lag" %in% colnames(data))
      stop("Using lag specific parameters requires lag to be a column of the data set.")

    if(any(!lagSpecific %in% colnames(isLagSpecific)))
      stop("Unknown element in lagSpecific. Expected one of the following: ",
           paste0(colnames(isLagSpecific), collapse = ", "))

    isLagSpecific[,lagSpecific] <- TRUE

  }

  if(! ("id" %in% colnames(data) & "occasion" %in% colnames(data))){
    stop("data must have columns with names id and occasion")
  }

  if(!is.character(latentIntercepts)){
    latentIntercepts <- apply(latentIntercepts, 2, as.character)
  }
  if(!is.character(ARCL)){
    ARCL <- apply(ARCL, 2, as.character)
  }
  if(!is.character(latentResiduals)){
    latentResiduals <- apply(latentResiduals, 2, as.character)
  }

  if(!is.character(loadings)){
    loadings <- apply(loadings, 2, as.character)
  }
  if(!is.character(manifestIntercepts)){
    manifestIntercepts <- apply(manifestIntercepts, 2, as.character)
  }
  if(!is.character(measurementError)){
    measurementError <- apply(measurementError, 2, as.character)
  }

  ids <- unique(data$id)
  occasions <- unique(data$occasion)

  if("lag" %in% colnames(data)){
    lags <- c()
    for(occasion in occasions){
      lag <- unique(data$lag[occasions == occasion])
      if(length(lag) != 1)
        stop("Found multiple lags for occasion ", occasion)
      lags <- c(lags, lag)
    }
  }else{
    lags <- NA
  }

  observations <- data[,!(colnames(data) %in% c("id", "occasion", "lag")), drop = FALSE]

  manifestNames <- colnames(observations)
  latentNames <- paste0("eta",1:ncol(loadings))

  if(ncol(observations) != nrow(loadings))
    stop("number of observed variables does not match the rows of the loadings matrix")


  makeLagSpecific <- function(x, lag){
    isParameter <- grepl(pattern = "[A-Za-z]+[A-Za-z0-9]*", x = x)
    if(any(isParameter))
      x[isParameter] <- paste0(x[isParameter],"_", lag)

    return(x)
  }

  latentIntercept <- vector("list", length(occasions))

  for(occasion in 1:length(latentIntercept)){

    if(occasion == 1 & any(latentIntercepts != "0") ){
      latentIntercepts_1 <- latentIntercepts
      for(ro in 1:nrow(latentIntercepts_1)){

        latentIntercepts_1[ro,1] <- paste0("latentIntercept_", latentNames[ro], "_init")

      }

      latentIntercept[[occasion]] <- latentIntercepts_1

    }else if(!isLagSpecific$latentIntercepts){
      latentIntercept[[occasion]] <- latentIntercepts
    }else{
      latentIntercept[[occasion]] <- apply(latentIntercepts, 2, makeLagSpecific, lag = lags[occasion])
    }
  }

  ARCLs <- vector("list", length(occasions))

  for(occasion in 1:length(ARCLs)){
    if(occasion == 1) next

    if(!isLagSpecific$ARCL){
      ARCLs[[occasion]] <- ARCL
    }else{
      ARCLs[[occasion]] <- apply(ARCL, 2, makeLagSpecific, lag = lags[occasion])
    }

  }

  latentResidual <- vector("list", length(occasions))

  for(occasion in 1:length(latentResidual)){

    if(occasion == 1){
      latentResiduals_1 <- latentResiduals
      for(ro in 1:nrow(latentResiduals_1)){
        for(co in ro:ncol(latentResiduals_1)){
          latentResiduals_1[ro,co] <- paste0("latentVariance_", latentNames[ro], "_", latentNames[co], "_init")
          latentResiduals_1[co,ro] <- paste0("latentVariance_", latentNames[ro], "_", latentNames[co], "_init")
        }
      }
      diag(latentResiduals_1) <- paste0("latentVariance_", latentNames, "_init")

      latentResidual[[occasion]] <- latentResiduals_1

    }else if(!isLagSpecific$latentResidual){
      latentResidual[[occasion]] <- latentResiduals
    }else{
      latentResidual[[occasion]] <- apply(latentResiduals, 2, makeLagSpecific, lag = lags[occasion])
    }
  }


  Loadings <- vector("list", length(occasions))

  for(occasion in 1:length(Loadings)){

    if(!isLagSpecific$loadings){
      Loadings[[occasion]] <- loadings
    }else{
      Loadings[[occasion]] <- apply(loadings, 2, makeLagSpecific, lag = lags[occasion])
    }
  }

  manifestIntercept <- vector("list", length(occasions))

  for(occasion in 1:length(manifestIntercept)){

    if(!isLagSpecific$manifestIntercepts){
      manifestIntercept[[occasion]] <- manifestIntercepts
    }else{
      manifestIntercept[[occasion]] <- apply(manifestIntercepts, 2, makeLagSpecific, lag = lags[occasion])
    }
  }

  MeasurementError <- vector("list", length(occasions))

  for(occasion in 1:length(MeasurementError)){

    if(!isLagSpecific$measurementError){
      MeasurementError[[occasion]] <- measurementError
    }else{
      MeasurementError[[occasion]] <- apply(measurementError, 2, makeLagSpecific, lag = lags[occasion])
    }
  }


  # create lavaan model
  modelSyntax <- ""

  makeRegressions <- function(latentIntercept, ARCL, latentNames, latentNames_OccMinus1){
    a <- c()
    for (latent in latentNames){
      a <- c(a, paste0(latent, " ~ ", latentIntercept[which(latentNames==latent),], "*1 + ",
                       paste0(paste0(ARCL[which(latentNames==latent),], "*", latentNames_OccMinus1), collapse = " + ")))
    }
    return(paste0(c("# regressions",a), collapse = "\n"))
  }

  makeLatentCovariances <- function(latentResidual, latentNames){
    l <- c()
    islower <- lower.tri(latentResidual, diag = TRUE)
    for (latent in latentNames){
      l <- c(l,
             paste0(latent, " ~~ ", paste0(paste0(latentResidual[which(latentNames==latent),], "*", latentNames)[islower[which(latentNames==latent),]],
                                           collapse = " + ")))
    }
    return(paste0(c("# latent residuals",l), collapse = "\n"))
  }

  makeMeasurements <- function(manifestIntercept, loadings, latentNames, manifestNames){
    l <- c()
    m <- c()
    for (latent in latentNames){
      l <- c(l, paste0(latent, " =~ ",paste0(paste0(loadings[,which(latentNames==latent)], "*", manifestNames),
                                             collapse = " + "))
      )
    }

    for(manifest in manifestNames){

      m <- c(m, paste0(manifest, " ~ ", manifestIntercept[which(manifestNames==manifest),],
                       "*1"))

    }

    return(
      paste0(
        paste0(c("# loadings", l), collapse = "\n"),
        "\n# intercepts\n",
        paste0(m, collapse = "\n")
      )
    )
  }

  makeMeasurementError <- function(measurementError, manifestNames){
    m <- c()
    islower <- lower.tri(measurementError, diag = TRUE)
    for (manifest in manifestNames){
      m <- c(m,
             paste0(manifest, " ~~ ",
                    paste0(paste0(measurementError[which(manifestNames==manifest),], "*", manifestNames)[islower[which(manifestNames==manifest),]],
                           collapse = " + ")))
    }
    return(paste0(c("# manifest resiudals",m), collapse = "\n"))
  }

  for(occ in 1:length(occasions)){

    manifestNamesOcc <- paste0(manifestNames, "_", occ)
    latentNamesOcc <- paste0(latentNames, "_", occ)

    modelSyntax <- paste0(modelSyntax,
                          "\n\n",
                          "#### Measurement Occasion ", occ, " ####\n")

    modelSyntax <- paste0(modelSyntax,
                          "\n",
                          makeMeasurements(manifestIntercept = manifestIntercept[[occ]],
                                           loadings = Loadings[[occ]],
                                           latentNames = latentNamesOcc,
                                           manifestNames = manifestNamesOcc)
    )

    modelSyntax <- paste0(modelSyntax,
                          "\n",
                          makeMeasurementError(measurementError = MeasurementError[[occ]],
                                               manifestNames = manifestNamesOcc)
    )

    if(occ != 1){
      modelSyntax <- paste0(modelSyntax,
                            "\n",
                            makeRegressions(
                              latentIntercept = latentIntercept[[occ]],
                              ARCL = ARCLs[[occ]],
                              latentNames = latentNamesOcc,
                              latentNames_OccMinus1 = latentNamesOccMinus1)
      )
    }

    modelSyntax <- paste0(modelSyntax,
                          "\n",
                          makeLatentCovariances(latentResidual = latentResidual[[occ]],
                                                latentNames = latentNamesOcc)
    )

    latentNamesOccMinus1 <- latentNamesOcc
  }

  # finally, let's re-organize the data
  dataWide <- data %>%
    tidyr::pivot_wider(
      id_cols = id,
      names_from = occasion,
      values_from = tidyr::all_of(manifestNames)) %>%
    ungroup() %>%
    select(!c(id))

  return(list("model" = modelSyntax,
              "data" = dataWide))
}
