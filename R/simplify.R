#' simplify
#'
#' simplify a model syntax created by CLPM-function
#' @param syntax model syntax from CLPM
#' @return simplified model syntax
#' @export
simplify <- function(syntax){

  # remove all things set to 0
  simpleSyntax <- stringr::str_remove_all(string = syntax,
                                          pattern = "0\\*[0-9a-zA-Z_]+[\\s]*\\+[\\s]*")
  simpleSyntax <- stringr::str_remove_all(string = simpleSyntax,
                                          pattern = "[\\+]*[\\s]*0\\*[0-9a-zA-Z_]+")
  simpleSyntax <- stringr::str_remove_all(string = simpleSyntax,
                                          pattern = "[a-zA-Z0-9]+[\\s]*[=~]+[\\s]*\\n")
  cat("Simplifying syntax.\n")
  cat(crayon::red("Note:"), "The simplified syntax is only meant to improve the readability of",
      "the model syntax. Please use the complex syntax when fitting the model with lavaan.")
  return(simpleSyntax)
}
