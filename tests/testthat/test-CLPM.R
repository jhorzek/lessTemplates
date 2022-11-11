test_that("CLPM works", {
  # The following simulation and analysis of a random intercept cross-lagged panel model
  # is based on the syntax from Jeroen D. Mulder & Ellen L. Hamaker (2021) Three Extensions of the Random
  # Intercept Cross-Lagged Panel Model, Structural Equation Modeling: A Multidisciplinary Journal,
  # 28:4, 638-648, DOI: 10.1080/10705511.2020.1784738
  #
  # See https://jeroendmulder.github.io/RI-CLPM/lavaan.html

  library(lessTemplates)
  library(lavaan)

  # Simulate Data
  data <- simulateRICLPM()

  # Set up model
  model <- "
# autoregressive and cross-lagged parameters:
eta1_(u) ~ a11_(u)*eta1_(u-1) + a12_(u)*eta2_(u-1)
eta2_(u) ~ a21_(u)*eta1_(u-1) + a22_(u)*eta2_(u-1)

# covariances
eta1_(u) ~~ 0*eta2_(u) + v11*eta1_(u)
eta2_(u) ~~ v22*eta2_(u)

# Add observations:
eta1_(u) =~ 1*y1_(u)
eta2_(u) =~ 1*y2_(u)

y1_(u) ~~ 0*y1_(u)
y2_(u) ~~ 0*y2_(u)

# random intercepts
RI_eta1 =~ 1*y1_(u)
RI_eta2 =~ 1*y2_(u)

RI_eta1 ~~ vri11*RI_eta1 + vri12*RI_eta2
RI_eta2 ~~ vri22*RI_eta2
"

  # create the lavaan syntax using lessTemplates:
  m <- lessTemplates::CLPM(model = model,
                           data = data,
                           addManifestVar = "no")
  # fit the model:
  fit <- sem(model = m$model,
             data = m$data,
             meanstructure = TRUE,
             missing = "ml")
  # get the parameter estimates
  coef(fit)

  # Manual specification adapted from https://jeroendmulder.github.io/RI-CLPM/lavaan.html
  model_comp <- "
# autoregressive and cross-lagged parameters:
eta1_u2 ~ a11_u2*eta1_u1 + a12_u2*eta2_u1
eta2_u2 ~ a21_u2*eta1_u1 + a22_u2*eta2_u1
eta1_u3 ~ a11_u3*eta1_u2 + a12_u3*eta2_u2
eta2_u3 ~ a21_u3*eta1_u2 + a22_u3*eta2_u2
eta1_u4 ~ a11_u4*eta1_u3 + a12_u4*eta2_u3
eta2_u4 ~ a21_u4*eta1_u3 + a22_u4*eta2_u3
eta1_u5 ~ a11_u5*eta1_u4 + a12_u5*eta2_u4
eta2_u5 ~ a21_u5*eta1_u4 + a22_u5*eta2_u4

# covariances
eta1_u1 ~~ initialCov_21*eta2_u1 + initialCov_11*eta1_u1
eta2_u1 ~~ initialCov_22*eta2_u1

eta1_u2 ~~ 0*eta2_u2 + v11*eta1_u2
eta2_u2 ~~ v22*eta2_u2

eta1_u3 ~~ 0*eta2_u3 + v11*eta1_u3
eta2_u3 ~~ v22*eta2_u3

eta1_u4 ~~ 0*eta2_u4 + v11*eta1_u4
eta2_u4 ~~ v22*eta2_u4

eta1_u5 ~~ 0*eta2_u5 + v11*eta1_u5
eta2_u5 ~~ v22*eta2_u5

# Add observations:
eta1_u1 =~ 1*y1_u1
eta2_u1 =~ 1*y2_u1

eta1_u2 =~ 1*y1_u2
eta2_u2 =~ 1*y2_u2

eta1_u3 =~ 1*y1_u3
eta2_u3 =~ 1*y2_u3

eta1_u4 =~ 1*y1_u4
eta2_u4 =~ 1*y2_u4

eta1_u5 =~ 1*y1_u5
eta2_u5 =~ 1*y2_u5

y1_u1 ~~ 0*y1_u1
y2_u1 ~~ 0*y2_u1

y1_u2 ~~ 0*y1_u2
y2_u2 ~~ 0*y2_u2

y1_u3 ~~ 0*y1_u3
y2_u3 ~~ 0*y2_u3

y1_u4 ~~ 0*y1_u4
y2_u4 ~~ 0*y2_u4

y1_u5 ~~ 0*y1_u5
y2_u5 ~~ 0*y2_u5

# random intercepts
RI_eta1 =~ 1*y1_u1
RI_eta2 =~ 1*y2_u1

RI_eta1 =~ 1*y1_u2
RI_eta2 =~ 1*y2_u2

RI_eta1 =~ 1*y1_u3
RI_eta2 =~ 1*y2_u3

RI_eta1 =~ 1*y1_u4
RI_eta2 =~ 1*y2_u4

RI_eta1 =~ 1*y1_u5
RI_eta2 =~ 1*y2_u5

RI_eta1 ~~ vri11*RI_eta1 + vri12*RI_eta2
RI_eta2 ~~ vri22*RI_eta2

RI_eta1 ~~ 0*eta1_u1 + 0*eta2_u1
RI_eta2 ~~ 0*eta1_u1 + 0*eta2_u1
"

  fit_comp <- sem(model = model_comp,
                  data = m$data,
                  meanstructure = TRUE,
                  missing = "ml")
  testthat::expect_equal(all(abs(coef(fit_comp)[names(coef(fit))] - coef(fit))<1e-5), TRUE)
})
