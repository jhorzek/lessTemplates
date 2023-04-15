#' simulateExample1
#'
#' Second order bivariate CLPM
#' @param seed seed used to simulate data
#' @return list with data and model
#' @export
simulateExample1 <- function(seed = 123){
  set.seed(seed)

  N <- 100
  U <- 10

  data <- c()

  for(i in 1:N){
    eta_1_um1 <- rnorm(1)
    eta_1_um2 <- rnorm(1)
    eta_2_um1 <- rnorm(1)
    eta_2_um2 <- rnorm(1)

    for(u in 1:U){
      eta_1_u <- .3*eta_1_um1 + .1*eta_1_um2 + rnorm(1,0,.2)
      eta_2_u <- .1*eta_1_um2 + .4*eta_2_um1 + rnorm(1,0,.2)

      y_1_u <- 1*eta_1_u + rnorm(1,0,.2)
      y_2_u <- .9*eta_1_u + rnorm(1,0,.2)
      y_3_u <- .9*eta_1_u + rnorm(1,0,.2)

      y_4_u <- 1*eta_2_u + rnorm(1,0,.2)
      y_5_u <- .9*eta_2_u + rnorm(1,0,.2)
      y_6_u <- .9*eta_2_u + rnorm(1,0,.2)

      data <- rbind(
        data,
        data.frame(
          person = i,
          occasion = u,
          y1 = y_1_u,
          y2 = y_2_u,
          y3 = y_3_u,
          y4 = y_4_u,
          y5 = y_5_u,
          y6 = y_6_u
        )
      )

      eta_1_um2 <- eta_1_um1
      eta_1_um1 <- eta_1_u
      eta_2_um2 <- eta_2_um1
      eta_2_um1 <- eta_2_u
    }
  }

  data_mod <- data[data$occasion %in% (U-5):U,]
  data_mod$occasion <- data_mod$occasion - (U-6)

  model <- "
eta1_(u) ~ a1*eta1_(u-1) + a2*eta1_(u-2)
eta2_(u) ~ b1*eta1_(u-2) + b2*eta2_(u-1)

eta1_(u) =~ 1*y1_(u) + l2*y2_(u) + l3*y3_(u)
eta2_(u) =~ 1*y4_(u) + l5*y5_(u) + l6*y6_(u)
"

  return(list(data = data_mod,
              model = model))

}

#' simulateRICLPM
#'
#' simulates data for a random intercept cross-lagged panel model based on the syntax from
#' Jeroen D. Mulder & Ellen L. Hamaker (2021) Three Extensions of the Random
#' Intercept Cross-Lagged Panel Model, Structural Equation Modeling: A Multidisciplinary Journal,
#' 28:4, 638-648, DOI: 10.1080/10705511.2020.1784738
#'
#' see https://jeroendmulder.github.io/RI-CLPM/lavaan.html
#' @param seed seed used to simulate data
#' @return data in long format
#' @export
simulateRICLPM <- function(seed = 123){
  set.seed(seed)

  model <- '
  # autoregressive and cross-lagged parameters:
  eta1_2 ~ .5*eta1_1 + -.2*eta2_1
  eta2_2 ~ .4*eta1_1 + .15*eta2_1

  eta1_3 ~ .5*eta1_2 + -.2*eta2_2
  eta2_3 ~ .4*eta1_2 + .15*eta2_2

  eta1_4 ~ .5*eta1_3 + -.2*eta2_3
  eta2_4 ~ .4*eta1_3 + .15*eta2_3

  eta1_5 ~ .5*eta1_4 + -.2*eta2_4
  eta2_5 ~ .4*eta1_4 + .15*eta2_4

  # initial covariances
  eta1_1 ~~ 1*eta1_1 + .5*eta2_1
  eta2_1 ~~ 1*eta2_1

  # covariances
  eta1_2 ~~ 0*eta2_2 + .25*eta1_2
  eta2_2 ~~ .25*eta2_2

  eta1_3 ~~ 0*eta2_3 + .25*eta1_3
  eta2_3 ~~ .25*eta2_3

  eta1_4 ~~ 0*eta2_4 + .25*eta1_4
  eta2_4 ~~ .25*eta2_4

  eta1_5 ~~ 0*eta2_5 + .25*eta1_5
  eta2_5 ~~ .25*eta2_5

  # Add observations:
  eta1_1 =~ 1*y1_1
  eta1_2 =~ 1*y1_2
  eta1_3 =~ 1*y1_3
  eta1_4 =~ 1*y1_4
  eta1_5 =~ 1*y1_5

  eta2_1 =~ 1*y2_1
  eta2_2 =~ 1*y2_2
  eta2_3 =~ 1*y2_3
  eta2_4 =~ 1*y2_4
  eta2_5 =~ 1*y2_5

  y1_1 ~~ 0*y1_1
  y1_2 ~~ 0*y1_2
  y1_3 ~~ 0*y1_3
  y1_4 ~~ 0*y1_4
  y1_5 ~~ 0*y1_5

  y2_1 ~~ 0*y2_1
  y2_2 ~~ 0*y2_2
  y2_3 ~~ 0*y2_3
  y2_4 ~~ 0*y2_4
  y2_5 ~~ 0*y2_5

  # random intercepts
  RI_eta1 =~ 1*y1_1 + 1*y1_2 + 1*y1_3 + 1*y1_4 + 1*y1_5
  RI_eta2 =~ 1*y2_1 + 1*y2_2 + 1*y2_3 + 1*y2_4 + 1*y2_5

  RI_eta1 ~~ .5*RI_eta1 + .1*RI_eta2
  RI_eta2 ~~ .5*RI_eta2
'

  data <- lavaan::simulateData(model = model, model.type = "sem", sample.nobs = 100)

  data_long <- tidyr::pivot_longer(data = cbind("person" = 1:nrow(data),
                                                data),
                                   cols = tidyr::starts_with(c("y1","y2")),
                                   names_to = c(".value", "occasion"),
                                   names_pattern = "(.)_(.)")

  colnames(data_long) <- c("person", "occasion", "y1", "y2")

  return(data_long)
}

#' simulateCLPM
#'
#' simulates data for a cross-lagged panel model with multiple indicators based on the syntax from
#' Jeroen D. Mulder & Ellen L. Hamaker (2021) Three Extensions of the Random
#' Intercept Cross-Lagged Panel Model, Structural Equation Modeling: A Multidisciplinary Journal,
#' 28:4, 638-648, DOI: 10.1080/10705511.2020.1784738
#'
#' see https://jeroendmulder.github.io/RI-CLPM/lavaan.html
#' @param seed seed used to simulate data
#' @return data in long format
#' @export
simulateCLPM <- function(seed = 123){
  set.seed(seed)

  model <- '
  # autoregressive and cross-lagged parameters:
  eta1_2 ~ .5*eta1_1 + -.2*eta2_1
  eta2_2 ~ .4*eta1_1 + .15*eta2_1

  eta1_3 ~ .5*eta1_2 + -.2*eta2_2
  eta2_3 ~ .4*eta1_2 + .15*eta2_2

  eta1_4 ~ .5*eta1_3 + 1*eta2_3
  eta2_4 ~ .4*eta1_3 + .15*eta2_3

  eta1_5 ~ .5*eta1_4 + 1*eta2_4
  eta2_5 ~ .4*eta1_4 + .15*eta2_4

  # initial covariances
  eta1_1 ~~ 1*eta1_1 + .5*eta2_1
  eta2_1 ~~ 1*eta2_1

  # covariances
  eta1_2 ~~ 0*eta2_2 + .25*eta1_2
  eta2_2 ~~ .25*eta2_2

  eta1_3 ~~ 0*eta2_3 + .25*eta1_3
  eta2_3 ~~ .25*eta2_3

  eta1_4 ~~ 0*eta2_4 + .25*eta1_4
  eta2_4 ~~ .25*eta2_4

  eta1_5 ~~ 0*eta2_5 + .25*eta1_5
  eta2_5 ~~ .25*eta2_5

  # Add observations:
  eta1_1 =~ 1*y1_1 + .9*y2_1 + .9*y3_1
  eta1_2 =~ 1*y1_2 + .9*y2_2 + .9*y3_2
  eta1_3 =~ 1*y1_3 + .9*y2_3 + .9*y3_3
  eta1_4 =~ 1*y1_4 + .9*y2_4 + .9*y3_4
  eta1_5 =~ 1*y1_5 + .9*y2_5 + .9*y3_5

  eta2_1 =~ 1*y4_1 + .9*y5_1 + .9*y6_1
  eta2_2 =~ 1*y4_2 + .9*y5_2 + .4*y6_2
  eta2_3 =~ 1*y4_3 + .9*y5_3 + .5*y6_3
  eta2_4 =~ 1*y4_4 + .9*y5_4 + .7*y6_4
  eta2_5 =~ 1*y4_5 + .9*y5_5 + .8*y6_5

  y1_1 ~~ .2*y1_1
  y1_2 ~~ .2*y1_2
  y1_3 ~~ .2*y1_3
  y1_4 ~~ .2*y1_4
  y1_5 ~~ .2*y1_5

  y2_1 ~~ .3*y2_1
  y2_2 ~~ .3*y2_2
  y2_3 ~~ .3*y2_3
  y2_4 ~~ .3*y2_4
  y2_5 ~~ .3*y2_5

  y3_1 ~~ .3*y3_1
  y3_2 ~~ .3*y3_2
  y3_3 ~~ .3*y3_3
  y3_4 ~~ .3*y3_4
  y3_5 ~~ .3*y3_5

  y4_1 ~~ .2*y4_1
  y4_2 ~~ .2*y4_2
  y4_3 ~~ .2*y4_3
  y4_4 ~~ .2*y4_4
  y4_5 ~~ .2*y4_5

  y5_1 ~~ .3*y5_1
  y5_2 ~~ .3*y5_2
  y5_3 ~~ .3*y5_3
  y5_4 ~~ .3*y5_4
  y5_5 ~~ .3*y5_5

  y6_1 ~~ .3*y6_1
  y6_2 ~~ .3*y6_2
  y6_3 ~~ .3*y6_3
  y6_4 ~~ .3*y6_4
  y6_5 ~~ .3*y6_5
'

  data <- lavaan::simulateData(model = model, model.type = "sem", sample.nobs = 100)

  data_long <- tidyr::pivot_longer(data = cbind("person" = 1:nrow(data),
                                                data),
                                   cols = tidyr::starts_with(c(paste0("y", 1:6))),
                                   names_to = c(".value", "occasion"),
                                   names_pattern = "(.)_(.)")

  colnames(data_long) <- c("person", "occasion", paste0("y", 1:6))

  return(data_long)
}

