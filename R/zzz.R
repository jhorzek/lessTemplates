#' lessTemplates
#'
#' The objective of lessTemplates is to provide a simplified implementation
#' for some common models which can be used in lessSEM. Currently included
#' are:
#'
#' 1.  cross-lagged panel models (with random intercepts).
#' 2.  continuous time structural equation models.
#' 3.  SEM with definition variables.
#'
#' Additionally, lessSEM also provides some transformations for these
#' models which allow for regularizing specific model structures (e.g., to
#' test measurement invariance in cross-lagged panel models).
#'
#' More details can be found here: https://github.com/jhorzek/lessTemplates
#'
#' @docType package
#' @author Jannik Orzek <orzek@mpib-berlin.mpg.de>
#' @importFrom methods is
#' @importFrom methods new
#' @importFrom methods slot
#' @importFrom stats rnorm
#' @importFrom rlang .data
#' @import lessSEM
#' @name lessTemplates
#' @keywords internal
"_PACKAGE"
