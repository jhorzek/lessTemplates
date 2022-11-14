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
                                    "#")
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

  # # check if left hand side of an equation has white space -> will be
  # # a data type and a variable name
  # isDefinition <- grepl(pattern = "[a-zA-Z:]+\\s+[a-zA-Z:]+\\s*=",
  #                          x = syntax) &
  #   !grepl(pattern = "start\\s*:\\s*",
  #         x = syntax)

  return(syntax)
}
