
library(magrittr)
library(bmgarch)
library(rstan)

# Stock data a bit too centered, and fuzzy.
## data(stocks)
## d <- stocks[1:500,]; rm(stocks)

# Ice river data
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

## bmOut <- bmgarch(d.train, parameterization = "CCC", iterations = 1000, meanstructure = "arma")
bmOut <- bmgarch(d.train, parameterization = "pdBEKK", iterations = 1000, meanstructure = "arma")

summary(bmOut, digits = 4)

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
# FOR Non-BEKK
## bmSamps_pred <- array(rnorm(n = 2000*295*3, bmSamps$mean, sqrt(bmSamps$var)), dim(bmSamps$mean))

# FOR BEKK
bmSamps_pred <- array(0, c(2000, 295, 3))
for(s in 1:2000) {
    for(n in 1:295) {
        bmSamps_pred[s, n,] <- mvtnorm::rmvnorm(1, bmSamps$mean[s, n,], bmSamps$var[s, n, , ])
    }
}

bmSamps_pred_sum <- apply(bmSamps_pred, c(2, 3), function(x){
    c(mean(x), quantile(x, .025), quantile(x, .975))
})
bmSamps_pred_sum <- aperm(bmSamps_pred_sum, c(2, 1, 3))

lines(bmSamps_pred_sum[,1,ts], col = "purple")
lines(bmSamps_pred_sum[,2,ts], col = "purple", lty = "dashed")
lines(bmSamps_pred_sum[,3,ts], col = "purple", lty = "dashed")

lines(296:300,forecast(bmOut, ahead = 5)$forecast$mean[,"mean",ts], col = "green")
lines(296:300,forecast(bmOut, ahead = 5)$forecast$mean[,"2.5%",ts], col = "green", lty = "dashed")
lines(296:300,forecast(bmOut, ahead = 5)$forecast$mean[,"97.5%",ts], col = "green", lty = "dashed")

# Plot variances
outVars <- as.matrix(out, pars = "logsd_hat") %>% apply(2, function(x){exp(x * 2)}) %>% colMeans %>% matrix(c(295, 3))
bmVars <- fitted(bmOut)$backcast$var[,"mean",]

par(mfrow=c(3,1))
plot(outVars[,1], type = 'l')
lines(bmVars[,1], col = "purple")
plot(outVars[,2], type = 'l')
lines(bmVars[,2], col = "purple")
plot(outVars[,3], type = 'l')
lines(bmVars[,3], col = "purple")

# How does (y - mu_hat)^2 plot against innovation?
lrv <- log(abs(d.train - t(matrix(summary(out, pars = "mu_hat")$summary[,"mean"], 3, nrow(d.train)))))
inn <- t(matrix(summary(out, pars = "var_innovation")$summary[,"mean"], 3, nrow(d.train)))
lsd <- t(matrix(summary(out, pars = "logsd_hat")$summary[,"mean"], 3, nrow(d.train)))
plot(lrv, inn)
plot(lrv, lsd)


####################
# BMGARCH Sim data #
####################

library(bmgarch)
library(rstan)

sim.bekk <- function(N,C,A,B, phi = NULL, theta = NULL) {
    if(ncol(C) != nrow(C)){
        stop("C must be symmetric, square, PD.")
    }
    if(ncol(A) != nrow(A)){
        stop("A must be square.")
    }
    if(ncol(B) != nrow(B)){
        stop("B must be square.")
    }
    nt <- ncol(C)

    y <- array(0, dim = c(N, nt))
    y[1,] <- rnorm(nt, 0, sqrt(diag(C)))

    H <- array(0, dim = c(nt, nt, N))
    H[,,1] <- C

    for(i in 2:N) {
        H[,,i] <- C + t(A) %*% (t(y[i - 1,, drop = FALSE]) %*% y[i - 1,,drop = FALSE]) %*% A + t(B) %*% H[,,i-1] %*% B
        y[i,] <- MASS::mvrnorm(1, rep(0, nt), H[,,i])
    }

    if (!is.null(phi) & !is.null(theta)) {
        ## Assume phi0 (intercept) is zero.
        if (ncol(phi) != nrow(phi)) {
            stop("phi must be square [nt, nt].")
        }
        if (ncol(theta) != nrow(theta)) {
            stop("theta must be square [nt, nt].")
        }
        if (ncol(phi) != nt) {
            stop("phi must be square [nt, nt].")
        }
        if (ncol(theta) != nt) {
            stop("theta must be square [nt, nt].")
        }
        mu <- array(0, dim = c(N, nt))
        mu[1,] <- 0
        for(i in 2:N) {
            mu[i,] <- 10 + y[i - 1, , drop = FALSE] %*% phi + (y[i - 1, ,drop = FALSE] - mu[i - 1,,drop = FALSE])%*%theta
            y[i,] <- y[i,,drop = FALSE] + mu[i,,drop = FALSE]
        }
        ## y <- mu + y
    }

    return(y)
}

# nt = 3
set.seed(13)
N <- 200
nt <- 3
C_sd <- diag(rep(2, 3))
C <- C_sd %*% rethinking::rlkjcorr(1,3, 5) %*% C_sd
A <- matrix(runif(nt^2, -.5, .5), ncol=nt)
B <- matrix(runif(nt^2, -.5, .5), ncol=nt)

# ARMA(1,1)
phi <- matrix(runif(nt^2, -.5, .5), ncol = nt)
theta <- matrix(runif(nt^2, -.5, .5), ncol = nt)
## phi <- matrix(0, ncol = nt, nrow = nt)
## theta <- matrix(0, ncol = nt, nrow = nt)
diag(phi) <- rep(.8, nt)
diag(theta) <- rep(.5, nt)

d <- sim.bekk(N, C, A, B, phi = NULL, theta =  NULL)
d.train <- d[1:190,]
d.test <- d[191:200,]


stan_data <- list(N = nrow(d.train),
                  F = ncol(d.train),
                  P = 0,
                  Q = 0,
                  AR_P = 1,
                  MA_P = 1,
                  AR_Q = 1,
                  MA_Q = 1,
                  y = d.train,
                  prior_only = 0,
                  ahead = 10
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

bmOut <- bmgarch(d.train, parameterization = "BEKK", distribution = "Gaussian", meanstructure = "arma")

summary(bmOut, digits = 4)

forecast(bmOut, ahead = 10)$forecast$mean
aperm(array(summary(out, pars = "pred")$summary, c(3, 10, 10)),c(2, 3, 1))

# Plot it

ts <- 3
plot(d[,ts], type = 'l')
abline(v = 190)
# Backcast
lines(summary(out, pars = paste0("backcast[",1:190,",",ts,"]"))$summary[,"mean"],col="blue")
lines(summary(out, pars = paste0("backcast[",1:190,",",ts,"]"))$summary[,"2.5%"],col="blue",lty='dashed')
lines(summary(out, pars = paste0("backcast[",1:190,",",ts,"]"))$summary[,"97.5%"],col="blue",lty='dashed')
# Predicted
lines(296:300, summary(out, pars = paste0("pred[",1:10,",",ts,"]"))$summary[,"mean"],col="red")
lines(296:300, summary(out, pars = paste0("pred[",1:10,",",ts,"]"))$summary[,"2.5%"],col="red", lty='dashed')
lines(296:300, summary(out, pars = paste0("pred[",1:10,",",ts,"]"))$summary[,"97.5%"],col="red", lty='dashed')

# Add bmgarch
bmSamps <- fitted(bmOut, inc_samples = TRUE)$backcast$samples

bmSamps_pred <- array(0, c(2000, 190, 3))
for(s in 1:2000) {
    for(n in 1:190) {
        bmSamps_pred[s, n,] <- mvtnorm::rmvnorm(1, bmSamps$mean[s, n,], bmSamps$var[s, n, , ])
    }
}

bmSamps_pred_sum <- apply(bmSamps_pred, c(2, 3), function(x){
    c(mean(x), quantile(x, .025), quantile(x, .975))
})
bmSamps_pred_sum <- aperm(bmSamps_pred_sum, c(2, 1, 3))

lines(bmSamps_pred_sum[,1,ts], col = "purple")
lines(bmSamps_pred_sum[,2,ts], col = "purple", lty = "dashed")
lines(bmSamps_pred_sum[,3,ts], col = "purple", lty = "dashed")

lines(296:300,forecast(bmOut, ahead = 10)$forecast$mean[,"mean",ts], col = "green")
lines(296:300,forecast(bmOut, ahead = 10)$forecast$mean[,"2.5%",ts], col = "green", lty = "dashed")
lines(296:300,forecast(bmOut, ahead = 10)$forecast$mean[,"97.5%",ts], col = "green", lty = "dashed")

# Plot variances
outVars <- as.matrix(out, pars = "logsd_hat") %>% apply(2, function(x){exp(x * 2)}) %>% colMeans %>% matrix(c(190, 3))
bmVars <- fitted(bmOut)$backcast$var[,"mean",]

par(mfrow=c(3,1))
plot(outVars[,1], type = 'l')
lines(bmVars[,1], col = "purple")
plot(outVars[,2], type = 'l')
lines(bmVars[,2], col = "purple")
plot(outVars[,3], type = 'l')
lines(bmVars[,3], col = "purple")
