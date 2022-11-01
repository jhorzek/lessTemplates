#' .CLPM
#'
#' creates a cross-laged panel model in RAM notation from a syntax similar to that
#' of lavaan.
#'
#' @param model string specifying the model syntax
#' @param data longitudinal data in long format. Must have the columns "person"
#' and "occasion"
#' @return Model in RAM notation
#' @keywords internal
.CLPM <- function(model, data){

  if(!all(c("person", "occasion") %in% colnames(data)))
    stop("Could not find the columns person and occasion in the data set.")
  if(!is(object = data, class2 = "data.frame"))
    data <- as.data.frame(data)

  nOccasions <- length(unique(data$occasion))

  message("Found ", nOccasions, " occasions in the data set.")

  # remove unnecessary white space
  syntax <- lessSEM:::.reduceSyntax(syntax = model)
  syntax <- lessTransformations:::.removeWhitespace(syntax = syntax)
  syntax <- lessTransformations:::.makeSingleLine(syntax = syntax)

  # find the names of all variables
  variableNames <- lessTransformations:::.getVariableNames(syntax = syntax)

  latents <- list(
    occasionDependent = variableNames$occasionDependent[!variableNames$occasionDependent %in% colnames(data)],
    fixed = variableNames$fixed[!variableNames$fixed %in% colnames(data)]
  )

  manifests <- list(
    occasionDependent = variableNames$occasionDependent[variableNames$occasionDependent %in% colnames(data)],
    fixed = variableNames$fixed[variableNames$fixed %in% colnames(data)]
  )

  nLatent <- data.frame("occasionDependent" = length(latents$occasionDependent),
                        "fixed" = length(latents$fixed))
  nManifest <- data.frame("occasionDependent" = length(manifests$occasionDependent),
                          "fixed" = length(manifests$fixed))
  # Now that we know the number of variables and the number of measurement occasions,
  # we can set up the matrices

  nRows <- nCols <- nManifest$occasionDependent * nOccasions + nManifest$fixed +
    nLatent$occasionDependent*nOccasions + nLatent$fixed
  varNames <- c(
    paste0(manifests$occasionDependent, "_u", rep(1:nOccasions, each = nManifest$occasionDependent)),
    manifests$fixed,
    paste0(latents$occasionDependent, "_u", rep(1:nOccasions, each = nLatent$occasionDependent)),
    latents$fixed
  )

  A <- matrix("0",
              nrow = nRows,
              ncol = nCols,
              dimnames = list(varNames,
                              varNames))

  S <- A

  M <- A[1,,drop = FALSE]

  Fmat <- matrix(0,
                 nrow = nManifest$occasionDependent * nOccasions + nManifest$fixed,
                 ncol = nCols,
                 dimnames = list(c(paste0(manifests$occasionDependent, "_u", rep(1:nOccasions, each = nManifest$occasionDependent)),
                                   manifests$fixed),
                                 varNames))
  for(i in 1:nrow(Fmat))
    Fmat[i,i] <- 1

  # Let's fill the matrices

  # get all occasion specific elements
  occasions <- data.frame(
    string = unique(unlist(stringr::str_extract_all(string = syntax,
                                                    pattern = "\\(u[-0-9]*\\)"))),
    regex = NA,
    evaluated = NA
  )
  # the following will be used to replace the (u-j) part with another string
  occasions$regex <- stringr::str_replace_all(occasions$string,
                                              pattern = "\\(",
                                              replacement = "\\\\(")
  occasions$regex <- stringr::str_replace_all(occasions$regex,
                                              pattern = "\\)",
                                              replacement = "\\\\)")

  for(u in nOccasions:1){

    # adapt syntax for this specific occasion
    syntax_u <- syntax
    occasions_u <- occasions

    for(i in 1:nrow(occasions)){

      # evaluate (u-j)
      occasions_u$evaluated[i] <- eval(parse(text = occasions_u$string[i]))

      # remove all elements with u below 1
      if(occasions_u$evaluated[i] < 1){
        syntax_u <- stringr::str_replace_all(string = syntax_u,
                                             pattern = paste0("[a-zA-Z0-9_\\(\\)]*[\\*]*[a-zA-Z0-9]+_",occasions_u$regex[i], "[\\+]*"),
                                             replacement = "")
        next
      }

      # replace (u-j) with _u-j
      syntax_u <- stringr::str_replace_all(string = syntax_u,
                                           pattern = occasions_u$regex[i],
                                           replacement = paste0("u", occasions_u$evaluated[i]))
    }

    for(i in 1:length(syntax_u)){
      splitted <- lessTransformations:::.splitEquation(equation = syntax_u[i])

      if(all(splitted$operator == "=~")){

        for(j in 1:nrow(splitted)){
          if(A[splitted$rhs[j], splitted$lhs[j]]!="0")
            stop("Redefinition of A[", splitted$rhs[j], ",", splitted$lhs[j], ".")
          A[splitted$rhs[j], splitted$lhs[j]] <- splitted$label[j]
        }

      }else if(all(splitted$operator == "~~")){

        for(j in 1:nrow(splitted)){
          if(S[splitted$rhs[j], splitted$lhs[j]]!="0")
            stop("Redefinition of S[", splitted$rhs[j], ",", splitted$lhs[j], ".")
          if(S[splitted$lhs[j], splitted$rhs[j]]!="0")
            stop("Redefinition of S[", splitted$lhs[j], ",", splitted$rhs[j], ".")
          S[splitted$rhs[j], splitted$lhs[j]] <- splitted$label[j]
          S[splitted$lhs[j], splitted$rhs[j]] <- splitted$label[j]
        }

      }else if(all(splitted$operator == "~" &  splitted$rhs == "1")){

        for(j in 1:nrow(splitted)){
          if(M[1,splitted$lhs[j]]!="0")
            stop("Redefinition of M[1,", splitted$lhs[j], ".")
          M[1,splitted$lhs[j]] <- splitted$label[j]
        }

      }else if(all(splitted$operator == "~")){

        for(j in 1:nrow(splitted)){
          if(A[splitted$lhs[j], splitted$rhs[j]]!="0")
            stop("Redefinition of A[", splitted$lhs[j], ",", splitted$rhs[j], ".")
          A[splitted$lhs[j], splitted$rhs[j]] <- splitted$label[j]
        }

      }else{

        stop("Unknown operator ", splitted$operator, ".")

      }

    }

  }

  # fill covariances of initial time points
  maxLag <- lessTransformations:::.findMaxLag(occasions = occasions)

  S <- .fillCovariances(S = S,
                        latents = latents,
                        maxLag = maxlag)

  RAM <- new("RAM")
  RAM@A <- A
  RAM@S <- S
  RAM@M <- M
  RAM@F <- Fmat
  RAM@manifest <- rownames(Fmat)
  RAM@latent <- colnames(Fmat)[!colnames(Fmat) %in% rownames(Fmat)]

  return(RAM)
}

#' .removeWhitespace
#'
#' remove any whitespace from a string
#' @param syntax string
#' @return string with removed whitespace
.removeWhitespace <- function(syntax){
  return(gsub(pattern = "\\s|\\t", replacement = "", x = syntax))
}

#' .makeSingleLine
#'
#' combine multi-line statements into one line
#' @param syntax string
#' @return string with combined multi-line-statements
.makeSingleLine <- function(syntax){

  for(i in 1:length(syntax)){
    if(grepl(pattern = "[+*]$", x = syntax[i])){
      if(i == length(syntax))
        stop("Could not parse syntax. Your syntax seems to end with a + or a *.")
      syntax[i] <- paste0(syntax[i], syntax[i+1])
      syntax[i+1] <- ""
    }
  }

  syntax <- syntax[syntax != ""]
  return(syntax)
}

#' .getVariableNames
#'
#' extracts the names of the variables from the syntax
#' @param syntax string
#' @return list with names of variables
#' @keywords internal
.getVariableNames <- function(syntax){

  # remove all parameters
  syntax_t <- gsub(pattern = "[a-zA-Z0-9]+\\*|[a-zA-Z0-9]+_\\(u[\\-]*[0-9]*\\)\\*",
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
  isOccasionDependent <- grepl(pattern = "_\\(u[\\-]*[0-9]*\\)",
                               x = variableNames)
  variableNames <- gsub(pattern = "_\\(u[\\-]*[0-9]*\\)",
                        replacement = "",
                        x = variableNames)
  names(isOccasionDependent) <- variableNames

  return(
    list(occasionDependent = unique(variableNames[isOccasionDependent]),
         fixed = unique(variableNames[!isOccasionDependent]))
  )
}


#' .splitEquation
#'
#' splits an equation in left hand side, right hand side, label, and operator
#' @param equation string with equation
#' @return data.frame with left hand side, right hand side, label, and operator
#' @keywords internal
.splitEquation <- function(equation){

  if(grepl("~~", equation)){
    operator <- "~~"
    equation <- stringr::str_replace(string = equation,
                                     pattern = "~~",
                                     replacement = "~=~")
  }else if(grepl("=~", equation)){
    operator <- "=~"
    equation <- stringr::str_replace(string = equation,
                                     pattern = "=~",
                                     replacement = "~=~")
  }else if(grepl("~", equation)){
    operator <- "~"
    equation <- stringr::str_replace(string = equation,
                                     pattern = "~",
                                     replacement = "~=~")
  }

  splitted <- unlist(stringr::str_split(string = equation, pattern = "~=~"))

  lhs <- splitted[1]
  rhs <- splitted[2]

  rhsElements <- c()
  labels <- c()
  while(rhs != ""){
    label <- NULL
    rhsElement <- NULL

    if(grepl(pattern = "^[a-zA-Z0-9_]+\\*", x = rhs)){
      # starts with an modifier
      label <- stringr::str_extract(string = rhs,
                                    pattern = "^[a-zA-Z0-9_]+\\*")
      label <- stringr::str_remove(string = label,
                                   pattern = "\\*")

      rhs <- stringr::str_remove(string = rhs,
                                 pattern = "^[a-zA-Z0-9_]+\\*")
    }

    # now that we removed the modifier, our equation starts with a variable
    rhsElement <- stringr::str_extract(string = rhs,
                                       pattern = "^[a-zA-Z0-9_]+")

    # remove that element
    rhs <- stringr::str_remove(string = rhs,
                               pattern = "^[a-zA-Z0-9_]+[\\+]*")

    if(is.null(label)){
      label <- paste0(lhs, operator, rhsElement)
    }


    labels <- c(labels, label)
    rhsElements <- c(rhsElements, rhsElement)
  }

  separated <- data.frame(
    lhs = lhs,
    operator = operator,
    label = labels,
    rhs = rhsElements
  )

  return(separated)

}

#' .fillCovariances
#'
#' fills the covariances between latent variables which are considered initial covariances
#' @param S S matrix with undirected effects
#' @param latents data.frame with names of latent variables
#' @param maxLag larges lag in the equations
#' @return updated S matrix
#' @keywords internal
.fillCovariances <- function(S, latents, maxLag){
  latentNames <- paste0(latents$occasionDependent, "_u", rep(1:maxLag, each = length(latents$occasionDependent)))
  covs <- matrix(paste0("initial_",
                        rep(1:length(latentNames), each = length(latentNames)),
                        rep(1:length(latentNames), length(latentNames))),
                 nrow = length(latentNames), ncol = length(latentNames),
                 byrow = TRUE
  )

  S[latentNames, latentNames] <- covs

  return(S)
}


#' .findMaxLag
#'
#' check the maximal time lag of the model
#' @param occasions data.frame with occasions
#' @return integer; larges time lag
#' @keywords internal
.findMaxLag <- function(occasions){
  removed_u <- stringr::str_remove_all(string = occasions$string,
                                       pattern = "\\(u[-]*|\\)")
  removed_u <- removed_u[removed_u!=""]
  return(max(as.numeric(removed_u)))
}
