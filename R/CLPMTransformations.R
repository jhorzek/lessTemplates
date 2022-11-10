#' transformCLPM
#'
#' Transform the parameters of a CLPM model. Returns a string with the
#' transformations to be passed to lessSEM.
#' @param CLPM model created with CLPM function
#' @param parameters names of the parameters which should be transformed. Note: Should be of form p_(u)
#' @param transformation which transformation should be used? Currently supported: "changepoint", or "measurementInvariance"
#' @return list with (1) transformation: string to be passed to lessSEM as transformation and (2) regularized: vector specifying which parameters should be regularized.
#' @examples
#' # The following simulation and analysis of a random intercept cross-lagged panel model
#' # is based on the syntax from Jeroen D. Mulder & Ellen L. Hamaker (2021) Three Extensions of the Random
#' # Intercept Cross-Lagged Panel Model, Structural Equation Modeling: A Multidisciplinary Journal,
#' # 28:4, 638-648, DOI: 10.1080/10705511.2020.1784738
#' #
#' # See https://jeroendmulder.github.io/RI-CLPM/lavaan.html
#'
#' library(lessTransformations)
#' library(lavaan)
#' library(lessSEM)
#'
#' # Simulate Data
#' data <- simulateRICLPM()
#'
#' # Set up model
#' model <- "
#' # autoregressive and cross-lagged parameters:
#' eta1_(u) ~ a11_(u)*eta1_(u-1) + a12_(u)*eta2_(u-1)
#' eta2_(u) ~ a21_(u)*eta1_(u-1) + a22_(u)*eta2_(u-1)
#'
#' # covariances
#' eta1_(u) ~~ 0*eta2_(u) + v11*eta1_(u)
#' eta2_(u) ~~ v22*eta2_(u)
#'
#' # Add observations:
#' eta1_(u) =~ 1*y1_(u)
#' eta2_(u) =~ 1*y2_(u)
#'
#' y1_(u) ~~ 0*y1_(u)
#' y2_(u) ~~ 0*y2_(u)
#'
#' # random intercepts
#' RI_eta1 =~ 1*y1_(u)
#' RI_eta2 =~ 1*y2_(u)
#'
#' RI_eta1 ~~ vri11*RI_eta1 + vri12*RI_eta2
#' RI_eta2 ~~ vri22*RI_eta2
#' "
#'
#' # create the lavaan syntax using lessTransformations:
#' m <- lessTransformations::CLPM(model = model,
#'                                data = data,
#'                                addManifestVar = "no")
#' # fit the model:
#' fit <- sem(model = m$model,
#'            data = m$data,
#'            meanstructure = TRUE,
#'            missing = "ml")
#' # get the parameter estimates
#' coef(fit)
#'
#' # In the simulation, we assumed that the autoregressive and cross-lagged parameters
#' # all stay the same over time (e.g, a11_u1 = a11_u2 = a11_u3 = a11_u4 = a11_u5).
#' # In practice, however, we won't know that. In the following, we will test this
#' # automatically
#' transform <- lessTransformations::transformCLPM(CLPM = m,
#'                                                 parameters = c("a11_(u)", "a12_(u)", "a21_(u)", "a22_(u)"),
#'                                                 # check measurement invariance of these parameters:
#'                                                 transformation = "measurementInvariance")
#' # fit using lasso regularization:
#' lf <- lasso(lavaanModel = fit,
#'             regularized = transform$regularized,
#'             nLambdas = 50,
#'             modifyModel = modifyModel(transformations = transform$transformation))
#' # let's plot the regularized parameters:
#' plot(lf)
#' # In the best model, all regularized parameters have been set to zero:
#' coef(lf, criterion = "BIC")
#' # We can conclude that the lasso regularization suggests invariance of the
#' # autoregressive and cross-lagged parameters
#' @export
transformCLPM <- function(CLPM,
                          parameters,
                          transformation){

  if(any(!grepl(pattern = "_\\(u\\)$", x = parameters))){
    stop("The parameters must all end with _(u) (e.g., a_(u)).")
  }

  # remove _(u)
  needsRemoval <- grepl(pattern = "_\\(u\\)$", x = parameters)
  parameters[needsRemoval] <- stringr::str_remove(string = parameters[needsRemoval],
                                                  pattern = "_\\(u\\)$")

  if(transformation == "changepoint"){

    return(.transformCLPMChangePoint(CLPM = CLPM,
                                     parameters = parameters,
                                     transformation = transformation))

  }else if(transformation == "measurementInvariance"){

    cat(crayon::red("Note:"), "The model will not be identified without regularization.")
    return(.transformCLPMMeasurementInvariance(CLPM = CLPM,
                                               parameters = parameters,
                                               transformation = transformation))

  }

  stop("Could not find a valid transformation. Possible are 'changepoint', or 'measurementInvariance'")
}

.transformCLPMChangePoint <- function(CLPM,
                                      parameters,
                                      transformation){

  transformations <- c()
  transformationParameters <- c()
  regularized <- c()

  for(p in 1:length(parameters)){
    # check location of parameter
    searchFor <- paste0(parameters[p],"_u[0-9]+")
    if(any(grepl(pattern = searchFor, x = CLPM$internal$RAM@A))){
      location <- "A"
    }else if(any(grepl(pattern = searchFor, x = CLPM$internal$RAM@S))){
      location <- "S"
    }else if(any(grepl(pattern = searchFor, x = CLPM$internal$RAM@M))){
      location <- "M"
    }else{
      stop("Could not find location of ", parameters[p], ".")
    }

    loc <- slot(CLPM$internal$RAM, name = location)

    # get all instances of this parameter
    instances <- sort(loc[grepl(pattern = searchFor, x = loc)])
    if(length(instances) == 1)
      stop("Cannot create differences for", instances, " because there is only one occasion for this instance.")

    for(i in 2:length(instances)){
      transformations <- c(transformations,
                           paste0(instances[i], " = ", instances[i-1], " + Delta_", instances[i])
      )
      transformationParameters <- c(transformationParameters,
                                    instances[i],
                                    instances[i-1],
                                    paste0("Delta_", instances[i])

      )
      regularized <- c(regularized, paste0("Delta_", instances[i]))
    }
  }

  transformationParameters <- unique(transformationParameters)

  transformationString <- paste0(
    "parameters: ", paste0(transformationParameters, collapse = ", "),
    "\n\n#Transformations:\n",
    paste0(transformations, collapse = "\n")
  )

  regularized <- unique(regularized)

  return(
    list("transformation" = transformationString,
         "regularized" = regularized)
  )
}



.transformCLPMMeasurementInvariance <- function(CLPM,
                                                parameters,
                                                transformation){

  transformations <- c()
  transformationParameters <- c()
  regularized <- c()

  for(p in 1:length(parameters)){
    # check location of parameter
    searchFor <- paste0(parameters[p],"_u[0-9]+")
    if(any(grepl(pattern = searchFor, x = CLPM$internal$RAM@A))){
      location <- "A"
    }else if(any(grepl(pattern = searchFor, x = CLPM$internal$RAM@S))){
      location <- "S"
    }else if(any(grepl(pattern = searchFor, x = CLPM$internal$RAM@M))){
      location <- "M"
    }else{
      stop("Could not find location of", parameters[p], ".")
    }

    loc <- slot(CLPM$internal$RAM, name = location)

    # get all instances of this parameter
    instances <- sort(loc[grepl(pattern = searchFor, x = loc)])
    if(length(instances) == 1)
      stop("Cannot create differences for", instances, " because there is only one occasion for this instance.")

    for(i in 1:length(instances)){
      transformations <- c(transformations,
                           paste0(instances[i], " = ", parameters[p], " + Delta_", instances[i])
      )
      transformationParameters <- c(transformationParameters,
                                    parameters[p],
                                    instances[i],
                                    paste0("Delta_", instances[i])

      )
      regularized <- c(regularized, paste0("Delta_", instances[i]))
    }
  }

  transformationParameters <- unique(transformationParameters)

  transformationString <- paste0(
    "parameters: ", paste0(transformationParameters, collapse = ", "),
    "\n\n#Transformations:\n",
    paste0(transformations, collapse = "\n")
  )

  regularized <- unique(regularized)

  return(
    list("transformation" = transformationString,
         "regularized" = regularized)
  )
}
