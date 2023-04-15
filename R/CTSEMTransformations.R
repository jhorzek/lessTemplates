#' transformCTSEM
#'
#' Transform the parameters of a CTSEM model. Returns a string with the
#' transformations to be passed to lessSEM.
#' @param CTSEM model created with CTSEM function
#' @param parameters names of the parameters which should be transformed. Note: Should be of form p(t)
#' @param transformation which transformation should be used? Currently supported: "changepoint", or "measurementInvariance"
#' @return list with (1) transformation: string to be passed to lessSEM as transformation and (2) regularized: vector specifying which parameters should be regularized.
#' @examples
#' # TODO
#' @export
transformCTSEM <- function(CTSEM,
                           parameters,
                           transformation){

  if(any(!grepl(pattern = "\\(t\\)$", x = parameters))){
    stop("The parameters must all end with (t) (e.g., l1(t)).")
  }

  cat(crayon::red("Note:"), " Currently, lessTemplates only supports transformations for the measurement model!")

  # NOTE: We can build on the transformations defined for CLPM as these are
  # largely identical for CTSEM
  # redefine transformations:

  # replace (t) with _(u):
  parameters <- stringr::str_replace_all(string = parameters,
                                         pattern = "\\(t\\)",
                                         replacement = "_\\(u\\)")

  transformationCLPM <- lessTemplates::transformCLPM(CLPM = CTSEM$internal$clpm,
                                                     parameters = parameters,
                                                     transformation = transformation)

  # combine the transformations required for CTSEM and the ones defined above:

  ## extract additional parameters:
  addParametersCLPM <- stringr::str_remove_all(
    stringr::str_extract(string = transformationCLPM$transformation, pattern = "parameters:[\\sa-zA-Z0-9_,]*"),
    pattern = "\\n")

  ## extract transformation for CLPM
  addTransformationCLPM <- stringr::str_remove_all(string = transformationCLPM$transformation,
                                                   pattern = "parameters:[\\sa-zA-Z0-9_,]*")

  ## add both to the existing transformations:
  transformations <- paste0(gsub(pattern = "parameters:", x = CTSEM$transformation,
                                 replacement = paste0(addParametersCLPM, ", ")),
                            "\n\n",
                            addTransformationCLPM)

  return(
    list(transformation = transformations,
         regularized = transformationCLPM$regularized)
  )
}
