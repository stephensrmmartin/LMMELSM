library(LMMELSM)

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
