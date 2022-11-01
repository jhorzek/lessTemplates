#' .reduceSyntax
#'
#' reduce user defined parameter transformation syntax to basic elements
#' @param syntax string with user defined transformations
#' @returns a cut and simplified version of the syntax
#' @keywords internal
.reduceSyntax <- function(syntax){

  # first, split rows and remove everything we don't need
  # split rows
  syntax <- stringr::str_split(string = syntax,
                               pattern = "\\n")[[1]]

  # remove comments
  hasComment <- stringr::str_locate(syntax,
                                    "#|!")
  for(i in 1:length(syntax)){
    if(!is.na(hasComment[i,1])){
      if(hasComment[i,1] == 1){
        syntax[i] <- ""
        next
      }
      syntax[i] <- stringr::str_trunc(
        syntax[i], width = hasComment[i,1]-1, ellipsis = ""
      )
    }
  }

  # remove empty elements
  syntax <- syntax[!grepl(pattern = "^\\s*$", x = syntax)]

  return(syntax)
}

#' .removeTrailingWhiteSpace
#'
#' removes trailing white space from reduces model syntax
#' @param reducedSyntax vector with syntax elements from .reduceSyntax
#' @return same vector, but without white spaces in the end
#' @keywords internal
.removeTrailingWhiteSpace <- function(reducedSyntax){
  # remove trailing white space
  for(i in 1:length(reducedSyntax)){
    if(grepl(pattern = "\\s$", x = reducedSyntax[i])){
      reducedSyntax[i] <- gsub(pattern = "\\s+$", replacement = "", x = reducedSyntax[i])
    }
  }
  return(reducedSyntax)
}


#' .mergeEquations
#'
#' combines equation elements which span multiple lines
#' @param reducedSyntax vector with syntax elements from .removeTrailingWhiteSpace
#' @return same vector, but with combined equations
#' @keywords internal
.mergeEquations <- function(reducedSyntax){
  for(i in 1:length(reducedSyntax)){
    if(grepl(pattern = "[~\\*\\+]$", x = reducedSyntax[i])){
      if(i == length(reducedSyntax))
        stop("Error while parsing your model. A line ends with + or ~")
      reducedSyntax[i] <- paste(reducedSyntax[i], reducedSyntax[i+1])
      reducedSyntax[i+1] <- NA
    }
  }
  reducedSyntax <- reducedSyntax[!is.na(reducedSyntax)]

  return(reducedSyntax)
}

##' .findMaxLag
##'
##' check the maximal time lag of the model
##' @param reducedSyntax vector with syntax elements from .mergeEquations
##' @return integer; larges time lag
##' @keywords internal
#.findMaxLag <- function(reducedSyntax){
#
#  maxLag <- 0
#  for(i in 1:length(reducedSyntax)){
#
#    if(grepl(pattern = "\\([ut]-[1-9]\\)", x = reducedSyntax[i])){
#
#      location <- gregexpr(pattern = "\\([ut]-[1-9]\\)",
#                           text = reducedSyntax[i])[[1]]
#      lengths <- attr(location, which = "match.length")
#
#      for(l in seq_len(length(location))){
#        elem <- substr(x = reducedSyntax[i],
#                       start = location[l],
#                       stop = location[l] + lengths[l] - 1)
#        lag <- as.integer(gsub(pattern = "[-\\s\\(u\\)*]", replacement = "", x = elem))
#
#        if(lag > maxLag){
#          maxLag <- lag
#        }
#      }
#
#    }
#
#  }
#
#  return(maxLag)
#}

#' .findLatents
#'
#' check the maximal time lag of the model
#' @param reducedSyntax vector with syntax elements from .mergeEquations
#' @return vector with latent names
#' @keywords internal
.findLatents <- function(reducedSyntax){

  # check for measurement model definition
  latentVars <- c()
  for(i in seq_len(length(reducedSyntax))){

    if(grepl(pattern = "=~", x = reducedSyntax[i])){
      location <- regexpr(pattern = "^[a-zA-Z0-9]+",
                          text = reducedSyntax[i])
      lengths <- attr(location, which = "match.length")
      latentVars <- c(latentVars,
                      substr(x = reducedSyntax[i],
                                         start = location,
                                         stop = location + lengths - 1))
    }

  }

  return(latentVars)

}
