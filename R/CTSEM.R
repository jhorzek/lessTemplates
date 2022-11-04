CTSEM <- function(model,
                  data){

  cat("\nSetting up a continuous time structural equation model.\n")

  if(!all(c("person", "time") %in% colnames(data)))
    stop("Could not find the columns person and time in the data set.")
  if(!is(object = data, class2 = "data.frame"))
    data <- as.data.frame(data)

  # remove unnecessary white space
  syntax <- lessSEM:::.reduceSyntax(syntax = model)
  syntax <- lessTransformations:::.removeWhitespace(syntax = syntax)
  syntax <- lessTransformations:::.makeSingleLine(syntax = syntax)

  # find statements regarding Wiener process
  wiener <- .getWienerProcessCTSEM(syntax)

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
  .getDiscreteBasis()


  RAM <- lessTransformations:::.CLPM(model = model, data = data)

}


.getWienerProcessCTSEM <- function(syntax = syntax){

  wiener <- unlist(stringr::str_extract_all(string = syntax,
                                            pattern = "d[0-9]*_W\\(t\\)"))
  if(length(wiener) == 0)
    stop("Could not find a statement regarding the wiener process (e.g., d_W(t)).")

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
  syntax_t <- gsub(pattern = "[a-zA-Z0-9]+\\*|[a-zA-Z0-9]+\\*",
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
  variableNames <- variableNames[!grepl(pattern = "d[0-9]*_W\\(t\\)",
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

.getDiscreteBasis <- function(latents, wiener){

  # keep all measurement equations
  measurementEquations <- syntax[grepl(pattern = "=~", x = syntax)]

  # find the highest order:
  highestOrder <- lessTransformations:::.getHighestOrderCTSEM(latents = latents,
                                                              wiener = wiener)

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

  getOrder <- function(ord){
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
