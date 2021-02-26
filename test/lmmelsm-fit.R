##############################
# Univariate, Random effects #
##############################

library(dplyr)
library(LMMELSM)

set.seed(14)
d <- LMMELSM:::simulate_lmmelsm(
                   n = 50,
                   K = 10,
                   lambda = c(.7, .7, .7, .8, .9),
                   resid = rep(1, 5),
                   nu = rep(0, 5),
                   mu_beta = matrix(c(.4, -.6), ncol = 1),
                   logsd_beta = matrix(c(.4, -.6), ncol = 1),
                   P_random_ind = c(1, 2),
                   Q_random_ind = c(1, 2),
                   mu_logsd_betas_cor = diag(1, 2 + 2 + 2),
                   mu_logsd_betas_sigma = rep(.3, 2 + 2 + 2),
                   epsilon_cor = matrix(1, 1, 1)
               )

sOut <- lmmelsm(list(factor1 ~ obs_1 + obs_2 + obs_3 + obs_4 + obs_5, location ~ loc_1 + loc_2 | loc_1 + loc_2, scale ~ sca_1 + sca_2 | sca_1 + sca_2), subject, d$df, iter = 500)

d$params$mu_logsd_betas_re[,1:2]
ranef(sOut)$location[,"Mean", drop = FALSE]

####################
# Multivariate, RE #
####################

library(LMMELSM)

set.seed(13)
n <- 50
K <- 10
F <- 2
J <- 8
lambda <- matrix(c(.8,.8,.8,.8,0,0,0,0,
                    0,0,0,0,.8,.8,.8,.8),byrow=TRUE,nrow = F)
resid <- rep(1, J)
nu <- rep(0, J)
mu_beta <- matrix(c(.4,.5,.6,.7), ncol = F)
logsd_beta <- matrix(c(.4,.5,.6,.7), ncol = F)
P_random_ind <- 1
Q_random_ind <- 2
mu_logsd_betas_cor <- diag(1, 2 * F + F + F)
mu_logsd_betas_sigma <- rep(.3, 2 * F + F + F)
epsilon_cor <- diag(1,2,2)

d <- LMMELSM:::simulate_lmmelsm(
        n = n,
        K = K,
        lambda = lambda,
        resid = resid,
        nu = nu,
        mu_beta = mu_beta,
        logsd_beta = logsd_beta,
        P_random_ind = P_random_ind,
        Q_random_ind = Q_random_ind,
        mu_logsd_betas_cor = mu_logsd_betas_cor,
        mu_logsd_betas_sigma = mu_logsd_betas_sigma,
        epsilon_cor = epsilon_cor
        )

fit <- lmmelsm(list(fac1 ~ obs_1 + obs_2 + obs_3 + obs_4,
                    fac2 ~ obs_5 + obs_6 + obs_7 + obs_8,
                    location ~ loc_1 + loc_2 | loc_1,
                    scale ~ sca_1 + sca_2 | sca_2),
               group = subject, data = d$df, iter = 200)

summary(fit)

ranefs <- ranef(fit)
coefs <- coef(fit)

ranefs$location_slope %>%
    filter(subject == 1)

coefs$location_slope %>%
    filter(subject == 1)
