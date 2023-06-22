#' .RAM2Lavaan
#'
#' @param RAM model of type RAM
#' @param meanstructure boolean: should a mean structure be modeled?
#' @return string to be used as model in lavaan
#' @keywords internal
.RAM2Lavaan <- function(RAM, meanstructure){
  if(!is(RAM, "RAM")) stop("RAM must be of class RAM.")

  regressions <- c("# Regressions")
  loadings <- c("\n# Loadings")
  covariances <- c("\n# (Co-)variances")
  intercepts <- c("\n# Intercepts")

  # let's make the difference between parameter labels and parameter values more
  # specific:
  RAM <- .valueOrLabel(RAM)

  for(l in RAM@latent){

    # loadings
    loadings_l <- RAM@A[RAM@manifest, l,drop=FALSE]
    loadings_l <- loadings_l[loadings_l != "0",,drop = FALSE]
    loadings <- c(
      loadings,
      paste0(l, " =~ ", paste0(paste0(loadings_l[,l], "*", rownames(loadings_l)), collapse = " + "))
    )

    # regressions
    if(all(RAM@A[l, c(RAM@manifest, RAM@latent)] == "0")) next
    regressions_l <- RAM@A[l, c(RAM@manifest, RAM@latent[RAM@latent != l]),drop=FALSE]
    regressions_l <- regressions_l[,regressions_l != "0",drop = FALSE]
    regressions <- c(
      regressions,
      paste0(l, " ~ ", paste0(paste0(regressions_l[l,], "*", colnames(regressions_l)), collapse = " + "))
    )

  }

  for(m in RAM@manifest){

    # regressions
    if(all(RAM@A[m, RAM@manifest] == "0")) next
    regressions_m <- RAM@A[m, c(RAM@manifest, RAM@latent[RAM@latent != l]),drop=FALSE]
    regressions_m <- regressions_m[,regressions_m != "0",drop = FALSE]
    regressions <- c(
      regressions,
      paste0(m, " ~ ", paste0(paste0(regressions_m[m,], "*", colnames(regressions_m)), collapse = " + "))
    )

  }

  for(i in c(RAM@manifest, RAM@latent)){

    # covariances
    S_i <- RAM@S[i,, drop = FALSE]
    select <- lower.tri(RAM@S, diag = TRUE)[which(rownames(RAM@S) == i),,drop = FALSE]
    # lavaan automatically adds more latent variables if covariances between latents
    # and manifests are specified, even if they are zero. So, we will also drop all
    # zero - covariances between these variables
    if(i %in% RAM@latent){
      select <- select & !((colnames(S_i) %in% RAM@manifest) & (S_i == "0"))
    }

    S_i <- S_i[select]
    covariances <- c(
      covariances,
      paste0(i, " ~~ ", paste0(paste0(S_i, "*", c(RAM@manifest,RAM@latent)[select]), collapse = " + "))
    )

    if(meanstructure){
      # intercepts
      intercepts <- c(
        intercepts,
        paste0(i, " ~ ", RAM@M[,i], "*1")
      )
    }

  }

  lavaanModel <- paste0(c(regressions,
                          loadings,
                          covariances,
                          intercepts),
                        collapse = "\n")

  return(lavaanModel)

}

#' .valueOrLabel
#'
#' separates between values and labels in RAM notation
#' @param RAM model of type RAM
#' @return model of type RAM with distinction between value and label
#' @keywords internal
.valueOrLabel <- function(RAM){

  isLabel <- function(x)
    grepl(pattern = "^[a-zA-Z]",x)

  for(i in 1:nrow(RAM@A)){
    for(j in 1:ncol(RAM@A)){
      if(isLabel(RAM@A[i,j])){
        RAM@A[i,j] <- paste0("label('",RAM@A[i,j], "')")
      }
      if(isLabel(RAM@S[i,j])){
        RAM@S[i,j] <- paste0("label('",RAM@S[i,j], "')")
      }
    }

    if(isLabel(RAM@M[1,i])){
      RAM@M[1,i] <- paste0("label('",RAM@M[1,i], "')")
    }

  }

  return(RAM)

}
