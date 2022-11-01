.simplify <- function(syntax){

  # remove all things set to 0
  simpleSyntax <- stringr::str_remove_all(string = syntax,
                                          pattern = "0\\*[0-9a-zA-Z_]+[\\s]*\\+[\\s]*")
  simpleSyntax <- stringr::str_remove_all(string = simpleSyntax,
                                          pattern = "[\\+]*[\\s]*0\\*[0-9a-zA-Z_]+")
  simpleSyntax <- stringr::str_remove_all(string = simpleSyntax,
                                          pattern = "[a-zA-Z0-9]+[\\s]*[=~]+[\\s]*\\n")
  return(simpleSyntax)
}
