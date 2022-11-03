CTSEM <- function(model,
                  data){

  cat("\nSetting up a continuous time structural equation model.\n")

  if(!all(c("person", "time") %in% colnames(data)))
    stop("Could not find the columns person and time in the data set.")
  if(!is(object = data, class2 = "data.frame"))
    data <- as.data.frame(data)

  nOccasions <- length(unique(data$occasion))

  # remove unnecessary white space
  syntax <- lessSEM:::.reduceSyntax(syntax = model)
  syntax <- lessTransformations:::.removeWhitespace(syntax = syntax)
  syntax <- lessTransformations:::.makeSingleLine(syntax = syntax)

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



  RAM <- lessTransformations:::.CLPM(model = model, data = data)

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

  # remove time indices
  isTimeDependent <- grepl(pattern = "\\(t\\)",
                               x = variableNames)
  variableNames <- gsub(pattern = "\\(t\\)",
                        replacement = "",
                        x = variableNames)
  variableNames <- gsub(pattern = "^d_",
                        replacement = "",
                        x = variableNames)
  names(isTimeDependent) <- variableNames

  return(
    list(occasionDependent = unique(variableNames[isTimeDependent]),
         fixed = unique(variableNames[!isTimeDependent]))
  )
}
