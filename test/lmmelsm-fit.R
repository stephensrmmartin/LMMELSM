library(LMMELSM)

##############
# Univariate #
##############
set.seed(14)
d <- LMMELSM:::simulate.uni.fe(
                     n = 20,
                     K = 200,
                     lambda = rep(.8, 10),
                     resid = rep(sqrt(1 - .8^2), 10),
                     nu = rep(0, 10),
                     mu_beta = c(.4, .8),
                     logsd_beta = c(.4, .8),
                     mu_logsd_cor = matrix(c(1, .5, .5, 1), 2, 2),
                     mu_logsd_sigma = c(.3, .3),
                     L2_pred_only = FALSE,
                     X_loc = NULL,
                     X_sca = NULL)

head(d$df)

sOut <- melsm_latent(list(myfactor ~ obs_1 + obs_2 + obs_3 + obs_4 + obs_5 + obs_6 + obs_7 + obs_8 + obs_9 + obs_10, location ~ loc_1 + loc_2, scale ~ sca_1 + sca_2), subject, d$df)
summary(sOut, pars = c("lambda", "nu", "sigma"))$summary
summary(sOut, pars = c("mu_logsd_random_sigma", "Omega_eta", "Omega_mean_logsd"))$summary
summary(sOut, pars = c("mu_beta","logsd_beta"))$summary

################
# Multivariate #
################

library(LMMELSM)
set.seed(14)
d <- LMMELSM:::simulate.multi.fe(
                   n = 20,
                   K = 200,
                   lambda = matrix(c(.7, .7, .7, .8, .9, 0, 0, 0, 0, 0,
                                     0, 0, 0, 0, 0, .7, .7, .7, .8, .9), nrow = 2, ncol = 10, byrow = TRUE),
                   resid = rep(1, 10),
                   nu = rep(0, 10),
                   mu_beta = matrix(c(.4, -.4,
                                      .6, -.6), nrow = 2, ncol = 2, byrow = TRUE),
                   logsd_beta = matrix(c(.4, -.4,
                                      .6, -.6), nrow = 2, ncol = 2, byrow = TRUE),
                   mu_logsd_cor = matrix(c(1, .6, 0, 0,
                                           .6, 1, 0, 0,
                                           0, 0, 1, .8,
                                           0, 0, .8, 1), ncol = 4, nrow = 4, byrow = TRUE),
                   mu_logsd_sigma = rep(.3, 4),
                   epsilon_cor = matrix(c(1, .5,
                                          .5, 1), nrow = 2, ncol = 2, byrow = TRUE),
                   L2_pred_only = FALSE,
                   X_loc = NULL,
                   X_sca = NULL
               )

sOut <- melsm_latent(list(factor1 ~ obs_1 + obs_2 + obs_3 + obs_4 + obs_5, factor2 ~ obs_6 + obs_7 + obs_8 + obs_9 + obs_10, location ~ loc_1 + loc_2, scale ~ sca_1 + sca_2), subject, d$df)
summary(sOut, pars = c("lambda", "nu", "sigma"))$summary
summary(sOut, pars = c("mu_logsd_random_sigma", "Omega_eta", "Omega_mean_logsd"))$summary
summary(sOut, pars = c("mu_beta","logsd_beta"))$summary


##############################
# Univariate, Random effects #
##############################

library(LMMELSM)

set.seed(14)
d <- LMMELSM:::simulate.multi.re(
                   n = 20,
                   K = 200,
                   lambda = matrix(c(.7, .7, .7, .8, .9), nrow = 1, ncol = 5, byrow = TRUE),
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

sOut <- melsm_latent(list(factor1 ~ obs_1 + obs_2 + obs_3 + obs_4 + obs_5, location ~ loc_1 + loc_2 | loc_1 + loc_2, scale ~ sca_1 + sca_2 | sca_1 + sca_2), subject, d$df)

summary(sOut, pars = c("lambda", "nu", "sigma"))$summary
summary(sOut, pars = c("mu_logsd_betas_random_sigma", "Omega_eta", "Omega_mean_logsd"))$summary
summary(sOut, pars = c("mu_beta","logsd_beta"))$summary
head(sort(summary(sOut, pars = c("mu_beta_random","logsd_beta_random"))$summary[,"Rhat"], decreasing = TRUE))


####################
# Multivariate, RE #
####################

library(LMMELSM)

set.seed(14)

d <- LMMELSM:::simulate.multi.re(
                   n = 50,
                   K = 50,
                   lambda = matrix(c(.7, .7, .7, .8, .9, 0, 0, 0, 0, 0,
                                     0, 0, 0, 0, 0, .7, .7, .7, .8, .9), nrow = 2, ncol = 10, byrow = TRUE),
                   resid = rep(1, 10),
                   nu = rep(0, 10),
                   mu_beta = matrix(c(.4, -.6), ncol = 2, nrow = 2),
                   logsd_beta = matrix(c(.4, -.6), ncol = 2, nrow = 2),
                   P_random_ind = c(1),
                   Q_random_ind = c(1),
                   mu_logsd_betas_cor = diag(1, 2 + 2 + 2*2),
                   mu_logsd_betas_sigma = rep(.3, 2 + 2 + 2*2),
                   epsilon_cor = matrix(c(1, -.4,
                                          -.4, 1), 2, 2)
               )

sOut <- melsm_latent(list(factor1 ~ obs_1 + obs_2 + obs_3 + obs_4 + obs_5,
                          factor2 ~ obs_6 + obs_7 + obs_8 + obs_9 + obs_10,
                          location ~ loc_1 + loc_2 | loc_1,
                          scale ~ sca_1 + sca_2 | sca_1), subject, d$df, iter = 1000)
