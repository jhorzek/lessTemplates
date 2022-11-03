#' transformCLPM
#'
#' Transform the parameters of a CLPM model. Returns a string with the
#' transformations to be passed to lessSEM.
#' @param CLPM model created with CLPM function
#' @param parameters names of the parameters which should be transformed. Note: Should be of form p, not p_u or p_(u)
#' @param transformation which transformation should be used? Currently supported: "changepoint"
#' @return list with (1) transformation: string to be passed to lessSEM as transformation and (2) regularized: vector specifying which parameters should be regularized.
transformCLPM <- function(CLPM,
                          parameters,
                          transformation){

  if(any(grepl(pattern = "_u[0-9]*$|_\\([u0-9\\-]+\\)$", x = parameters))){
    isIncorrect <- grepl(pattern = "_u[0-9]*$|_\\([u0-9\\-]+\\)$", x = parameters)
    fixed <- parameters
    fixed[isIncorrect] <- stringr::str_remove(string = fixed[isIncorrect],
                                              pattern = "_u[0-9]*$|_\\([u0-9\\-]+\\)$")
    warning("The parameters argument may be misspecified.\nReceived: ", paste0(parameters[isIncorrect], collapse = ", "),
            "\nFixed: ", paste0(fixed, collapse = ", "), "\nPlease make sure that the fixed version is correct.")
    parameters <- fixed
  }

  if(transformation == "changepoint"){
    return(.transformCLPMChangePoint(CLPM = CLPM,
                                     parameters = parameters,
                                     transformation = transformation))
  }

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
      stop("Could not find location of", parameters[p], ".")
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
