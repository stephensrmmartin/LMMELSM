
library(bmgarch)
library(rstan)

# Stock data a bit too centered, and fuzzy.
## data(stocks)
## d <- stocks[1:500,]; rm(stocks)

data(package="tseries", ice.river)
d <- ice.river[1:300,c(1,2,4)]
d[,1:2] <- log(d[,1:2])
d <- scale(d)

d.train <- d[1:(nrow(d)-5),]
d.test <- d[(nrow(d)-4):nrow(d),]


stan_data <- list(N = nrow(d.train),
                  F = 3,
                  P = 0,
                  Q = 0,
                  AR_P = 1,
                  MA_P = 1,
                  AR_Q = 1,
                  MA_Q = 1,
                  y = d.train,
                  prior_only = 0,
                  ahead = 5
                  )

sm <- stan_model("melsmArmaTest.stan")

out <- sampling(sm, data = stan_data, cores = 4, iter = 1000, init = 0)

pars <- c("ar_location",
          "ma_location",
          "ar_scale",
          "ma_scale",
          ## "var_innovation_sd",
          "epsilon_cor_L",
          "nu","sigma",
          "pred")

summary(out, pars = pars)$summary

bmOut <- bmgarch(d.train, parameterization = "CCC", iterations = 1000, meanstructure = "arma")

summary(bmOut)

# Forecasts
forecast(bmOut, ahead = 5)$forecast$mean
aperm(array(summary(out, pars = "pred")$summary, c(3, 5, 10)),c(2, 3, 1))

# Plot it
ts <- 3
plot(d[,ts], type = 'l')
abline(v = 295)
# Backcast
lines(summary(out, pars = paste0("backcast[",1:295,",",ts,"]"))$summary[,"mean"],col="blue")
lines(summary(out, pars = paste0("backcast[",1:295,",",ts,"]"))$summary[,"2.5%"],col="blue",lty='dashed')
lines(summary(out, pars = paste0("backcast[",1:295,",",ts,"]"))$summary[,"97.5%"],col="blue",lty='dashed')
# Predicted
lines(296:300, summary(out, pars = paste0("pred[",1:5,",",ts,"]"))$summary[,"mean"],col="red")
lines(296:300, summary(out, pars = paste0("pred[",1:5,",",ts,"]"))$summary[,"2.5%"],col="red", lty='dashed')
lines(296:300, summary(out, pars = paste0("pred[",1:5,",",ts,"]"))$summary[,"97.5%"],col="red", lty='dashed')

# Add bmgarch
bmSamps <- fitted(bmOut, inc_samples = TRUE)$backcast$samples
bmSamps_pred <- array(rnorm(n = 2000*295*3, bmSamps$mean, sqrt(bmSamps$var)), dim(bmSamps$mean))
bmSamps_pred_sum <- apply(bmSamps_pred, c(2, 3), function(x){
    c(mean(x), quantile(x, .025), quantile(x, .975))
})
bmSamps_pred_sum <- aperm(bmSamps_pred_sum, c(2, 1, 3))

lines(bmSamps_pred_sum[,1,ts], col = "purple")
lines(bmSamps_pred_sum[,2,ts], col = "purple", lty = "dashed")
lines(bmSamps_pred_sum[,3,ts], col = "purple", lty = "dashed")

## lines(forecast(bmOut, ahead = 5)$backcast$mean[,"mean",ts], col = "green")
## lines(forecast(bmOut, ahead = 5)$backcast$mean[,"2.5%",ts], col = "green", lty = "dashed")
## lines(forecast(bmOut, ahead = 5)$backcast$mean[,"97.5%",ts], col = "green", lty = "dashed")

lines(296:300,forecast(bmOut, ahead = 5)$forecast$mean[,"mean",ts], col = "green")
lines(296:300,forecast(bmOut, ahead = 5)$forecast$mean[,"2.5%",ts], col = "green", lty = "dashed")
lines(296:300,forecast(bmOut, ahead = 5)$forecast$mean[,"97.5%",ts], col = "green", lty = "dashed")

# How does (y - mu_hat)^2 plot against innovation?
lrv <- log(abs(d.train - t(matrix(summary(out, pars = "mu_hat")$summary[,"mean"], 3, nrow(d.train)))))
inn <- t(matrix(summary(out, pars = "var_innovation")$summary[,"mean"], 3, nrow(d.train)))
lsd <- t(matrix(summary(out, pars = "logsd_hat")$summary[,"mean"], 3, nrow(d.train)))
plot(lrv, inn)
plot(lrv, lsd)


library(brms)

## brmOut <- brm(toyota ~ arma(p = 1, q = 1), data = d, cores = 4, iter = 1000)
brmOut <- brm(mvbrmsformula(flow.vat ~ arma(p = 1, q = 1),
                            flow.jok ~ arma(p = 1, q = 1),
                            temp ~ arma(p = 1, q = 1)), data = d.train, cores = 4, iter = 1000)

summary(brmOut)
