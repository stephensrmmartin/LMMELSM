##############################
# Univariate, Random effects #
##############################

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
