#' simulateExample1
#'
#' Second order bivariate CLPM
#' @return list with data and model
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
