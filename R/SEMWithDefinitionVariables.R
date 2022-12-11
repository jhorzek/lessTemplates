#' SEMWithDefinitionVariables
#'
#' Allows for SEM with definition variables. This is currently not supported
#' by lavaan, but a small trick allows for setting up such models with lessSEM.
#' The basic idea is to use the multi-group implementation. The function will
#' return a vector of lavaan models. These should then be passed to any of the
#' optimizer functions of lessSEM (e.g., to bfgs() if no regularization is required).
#'
#' @param lavaanSyntax model syntax. This syntax should include the definition
#' variables as modifiers. The function will replace these with the values
#' defined in the definitionVariables
#' @param definitionVariables data.frame with definition variables. Must have
#' values for each person and must have column names. All column names must
#' occur in the lavaanSytax as modifiers
#' @param ... additional arguments passed to lavaan
#' @examples
#' # We will demonstrate the use of definition variables with a latent growth curve
#' # model. We assume that the six measurement occasions took place at subject-specific
#' # time points.
#'
#' #### Data simulation #####
#' # You can ignore the following; we just simulate some data to be used
#' # for our model.
#' ## Population parameters ##
#' intercept_mu <- 0
#' intercept_sigma <- 1
#' slope_mu <- .3
#' slope_sigma <- 1
#'
#' ## data set ##
#' N <- 50
#' intercepts <- rnorm(n = N,
#'                     mean = intercept_mu,
#'                     sd = intercept_sigma)
#' slopes <- rnorm(n = N,
#'                 mean = slope_mu,
#'                 sd = slope_sigma)
#' times <- as.data.frame(
#'   matrix(seq(0,5,1),
#'          nrow = N,
#'          ncol = 6,
#'          byrow = TRUE,
#'          dimnames = list(NULL, paste0("t", 0:5))) +
#'     cbind(0,matrix(round(runif(n = N*5, min = -.2,max = .2),2),
#'                    nrow = N,
#'                    ncol = 5,
#'                    byrow = TRUE)) # we add some jitter to make the times person-specific
#' )
#'
#' lgcData <- matrix(NA,
#'                   nrow = N,
#'                   ncol = ncol(times),
#'                   dimnames = list(NULL, paste0("x", 0:5)))
#'
#' for(i in 1:N){
#'   lgcData[i,] <- intercepts[i] + unlist(times[i,])* slopes[i] + rnorm(ncol(lgcData),0,.3)
#' }
#' lgcData <- as.data.frame(lgcData)
#'
#' #### Let's have a look at the data ####
#' # lgcData is a data.frame with the observed values at the six measurement
#' # occasions:
#' head(lgcData)
#'
#' # times is a data.frame with the time points at which the six measurements
#' # took place for each person. Note that these differ between individuals.
#' # These are going to be our definition variables.
#' head(times)
#'
#' #### Set up the base model ####
#'
#' model <- "
#' int =~ 1*x0 + 1*x1 + 1*x2 + 1*x3 + 1*x4 + 1*x5
#' slope =~ t0*x0 + t1*x1 + t2*x2 + t3*x3 + t4*x4 + t5*x5
#'
#' int ~ intMean*1
#' slope ~ slopeMean*1
#'
#' int ~~ intVar*int + intSlopeCov*slope
#' slope ~~ slopeVar*slope
#'
#' x0 ~~ v*x0; x1 ~~ v*x1; x2 ~~ v*x2; x3 ~~ v*x3; x4 ~~ v*x4; x5 ~~ v*x5
#'
#' x0 ~ 0*1; x1 ~ 0*1; x2 ~ 0*1; x3 ~ 0*1; x4 ~ 0*1; x5 ~ 0*1
#' "
#'
#' # The most important part of our model definition is that we used the
#' # modifiers t0, t1, t2, t3, t4, and t5. These have the same names as our
#' # defintion variables in the times variable. lessTemplates with
#' # replace the modifiers with the values found in the times variable.
#'
#' #### Set up the model with definition variables ####
#'
#' dvarModel <- SEMWithDefinitionVariables(lavaanSyntax = model,
#'                                         definitionVariables = times,
#'                                         data = lgcData)
#'
#' # fit with lessSEM
#' library(lessSEM)
#' dvarFit <- lessSEM::bfgs(lavaanModel = dvarModel)
#' dvarFit@parameters
#' @export
SEMWithDefinitionVariables <- function(lavaanSyntax,
                                       definitionVariables,
                                       data,
                                       ...){
  if(!is(lavaanSyntax, "character"))
    stop("lavaanSyntax must be a string specifying a lavaan model.")
  if(!is(definitionVariables, "data.frame"))
    stop("definitionVariables must be of class data.frame")
  if(is.null(colnames(definitionVariables)))
    stop("definitionVariables must have column names")
  if(nrow(definitionVariables) != nrow(data))
    stop("definitionVariables must have as many rows as there are subjects (rows) in the data set.")
  if(anyNA(definitionVariables))
    stop("NAs are not allows in definitionVariables")

  # check if all definition variables are used in the model
  for(i in colnames(definitionVariables)){
    if(!stringr::str_detect(string = lavaanSyntax,
                            pattern = paste0("[~\\+]+[\\s]*",i, "[\\s]*\\*"))){
      stop("Could not find definition variable ", i, " in the model. Checked for the ",
           "following pattern: ", paste0(i, "*"))
    }
  }

  # We don't have to specify a separate model for each person;
  # a separate model for each pattern of definition variables is enough.
  uniqueDefinitionVariables <- definitionVariables[!duplicated(definitionVariables),]
  # find out which person belongs to which unique pattern of definition variables
  pattern <- apply(definitionVariables, 1, function(x) which(apply(uniqueDefinitionVariables, 1, function(y) all(x == y))))

  models <- vector("list", nrow(uniqueDefinitionVariables))

  for(ptrn in unique(pattern)){

    it <- which(unique(pattern) == ptrn)

    syntaxPtrn <- lavaanSyntax

    for(i in colnames(definitionVariables)){
      # replace definition variables with values from data frame
      toReplace <- unlist(stringr::str_extract_all(string = syntaxPtrn,
                                                   pattern = paste0("[~\\+]+[\\s]*",i, "[\\s]*\\*")))
      replaceWith <- toReplace
      for(j in length(toReplace))
        replaceWith[j] <- stringr::str_replace_all(string = toReplace[j],
                                                   pattern = i,
                                                   replacement = paste0(uniqueDefinitionVariables[it,i]))
      for(j in 1:length(toReplace))
        syntaxPtrn <- gsub(pattern = toReplace[j],
                           replacement = replaceWith[j],
                           x = syntaxPtrn,
                           fixed = TRUE)
    }

    # Set up lavaan model. Note that we will have to pass the entire data set to
    # lavaan because lavaan does not like setting up models with N=1.
    # We will replace the data set in the next step
    models[[it]] <- lavaan::lavaan(model = syntaxPtrn,
                                 data = data,
                                 do.fit = FALSE,
                                 meanstructure = TRUE,
                                 ...)
    # replace data
    # Note: lavaan tends to change the order of the data; therefore we
    # first extract the data from the model
    internalData <- lavaan::lavInspect(models[[it]], "data")
    # replace the data set
    models[[it]]@Data@X[[1]] <- as.matrix(data[which(pattern == ptrn),
                                             colnames(internalData),
                                             drop = FALSE])
  }

  return(models)

}
