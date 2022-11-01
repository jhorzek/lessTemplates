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


#' RAM2Lavaan
#'
#' @param RAM model of type RAM
#' @return string to be used as model in lavaan
RAM2Lavaan <- function(RAM){
  if(!is(RAM, "RAM")) stop("RAM must be of class RAM.")

  regressions <- c("# Regressions")
  loadings <- c("\n# Loadings")
  covariances <- c("\n# (Co-)variances")
  intercepts <- c("\n# Intercepts")

  for(l in RAM@latent){

    # loadings
    loadings <- c(
      loadings,
      paste0(l, " =~ ", paste0(paste0(RAM@A[RAM@manifest, l], "*", manifest), collapse = " + "))
    )

    # regressions
    if(all(RAM@A[l, c(RAM@manifest, RAM@latent)] == "0")) next
    regressions <- c(
      regressions,
      paste0(l, " ~ ", paste0(paste0(RAM@A[l, c(RAM@manifest, RAM@latent)], "*", c(RAM@manifest, RAM@latent)), collapse = " + "))
    )

  }

  for(m in RAM@manifest){

    # regressions
    if(all(RAM@A[l, RAM@manifest] == "0")) next
    regressions <- c(
      regressions,
      paste0(m, " ~ ", paste0(paste0(RAM@A[l, RAM@manifest], "*", RAM@manifest), collapse = " + "))
    )

  }

  for(i in c(RAM@manifest, RAM@latent)){

    # covariances
    S_i <- S[i,, drop = FALSE]
    select <- lower.tri(S, diag = TRUE)[which(rownames(S) == i),,drop = FALSE]
    S_i <- S_i[select]
    covariances <- c(
      covariances,
      paste0(i, " ~~ ", paste0(paste0(S_i, "*", c(RAM@manifest,RAM@latent)[select]), collapse = " + "))
    )

    # intercepts
    intercepts <- c(
      intercepts,
      paste0(i, " ~ ", RAM@M[,i], "*1")
    )

  }

  lavaanModel <- paste0(c(regressions,
                           loadings,
                           covariances,
                           intercepts),
                         collapse = "\n")

  return(lavaanModel)

}
