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
  RAM <- lessTemplates:::.valueOrLabel(RAM)

  for(l in RAM@latent){

    # loadings
    loadings <- c(
      loadings,
      paste0(l, " =~ ", paste0(paste0(RAM@A[RAM@manifest, l], "*", RAM@manifest), collapse = " + "))
    )

    # regressions
    if(all(RAM@A[l, c(RAM@manifest, RAM@latent)] == "0")) next
    regressions <- c(
      regressions,
      paste0(l, " ~ ", paste0(paste0(RAM@A[l, c(RAM@manifest, RAM@latent[RAM@latent != l])], "*", c(RAM@manifest, RAM@latent[RAM@latent != l])), collapse = " + "))
    )

  }

  for(m in RAM@manifest){

    # regressions
    if(all(RAM@A[l, RAM@manifest] == "0")) next
    regressions <- c(
      regressions,
      paste0(m, " ~ ", paste0(paste0(RAM@A[m, RAM@manifest[RAM@manifest != m]], "*", RAM@manifest[RAM@manifest != m]), collapse = " + "))
    )

  }

  for(i in c(RAM@manifest, RAM@latent)){

    # covariances
    S_i <- RAM@S[i,, drop = FALSE]
    select <- lower.tri(RAM@S, diag = TRUE)[which(rownames(RAM@S) == i),,drop = FALSE]
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
