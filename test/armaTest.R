
library(bmgarch)
library(rstan)

data(stocks)
d <- stocks[1:500,]; rm(stocks)

stan_data <- list(N = nrow(d),
                  F = 3,
                  P = 0,
                  Q = 0,
                  AR_P = 1,
                  MA_P = 1,
                  AR_Q = 1,
                  MA_Q = 1,
                  y = d[,3:5],
                  prior_only = 0,
                  ahead = 8
                  )

sm <- stan_model("melsmArmaTest.stan")

out <- sampling(sm, data = stan_data, cores = 4, iter = 1000)

pars <- c("ar_location",
          "ma_location",
          "ar_scale",
          "ma_scale",
          "var_innovation_sd",
          "epsilon_cor_L",
          "nu","sigma",
          "pred")

summary(out, pars = pars)$summary

bmOut <- bmgarch(d[,3:5], parameterization = "CCC", iterations = 1000, meanstructure = "arma")

summary(bmOut)

library(brms)

## brmOut <- brm(toyota ~ arma(p = 1, q = 1), data = d, cores = 4, iter = 1000)
brmOut <- brm(mvbrmsformula(toyota ~ arma(p = 1, q = 1),
                            nissan ~ arma(p = 1, q = 1),
                            honda ~ arma(p = 1, q = 1)), data = d, cores = 4, iter = 1000)

summary(brmOut)
