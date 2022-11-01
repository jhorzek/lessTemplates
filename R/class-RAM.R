#' SEM in RAM notation
#' @slot A character matrix with directed effects
#' @slot S character matrix with undirected effects
#' @slot M character matrix with intercepts
#' @slot F numeric filter matrix
#' @slot latent vector with names of latent variables
#' @slot manifest vector with names of manifest variables
setClass("RAM",
         representation = representation(
           A = "matrix",
           S = "matrix",
           M = "matrix",
           F = "matrix",
           latent = "character",
           manifest = "character")
)
