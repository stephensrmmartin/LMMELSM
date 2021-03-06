library(LMMELSM)

set.seed(13)

# Settings
n <- 50
K <- 100
J <- 12
F <- 2
P <- 2
P_random <- 1
Q <- 2
Q_random <- 1
R <- 1

fac1_loadings <- c(rep(.6, J/2), rep(0, J/2))
fac2_loadings <- c(rep(0, J/2), rep(.8, J/2))
lambda <- rbind(fac1_loadings, fac2_loadings)
resid <- rep(.8, J)
nu <- rep(0, J)

mu_beta <- matrix(c(.6, -.6,
                    -.8, .8), ncol = F, byrow = TRUE)

logsd_beta <- matrix(c(-.6, .6,
                       .8, -.8), ncol = F, byrow = TRUE)

zeta <- matrix(c(0, .4, 0, .5, rep(0, F*P_random + F*Q_random)), nrow = 1)

re_sigmas <- rep(.4, 2*F + F*P_random + F*Q_random)
re_cor <- rethinking::rlkjcorr(1, 2*F + F*P_random + F*Q_random)

group <- rep(1:K, each = n)
X_l2 <- rnorm(K)
X_l1_loc <- rnorm(n*K)
X_l1_sca <- rnorm(n*K)

X_bet <- matrix(X_l2[group], ncol = 1)
X_loc <- cbind(X_l1_loc, X_bet)
X_sca <- cbind(X_l1_sca, X_bet)

factor_cor <- matrix(c(1, -.5, -.5, 1), F, F)

d <- LMMELSM:::simulate_lmmelsm(n,
                                       K,
                                       lambda,
                                       resid,
                                       nu,
                                       mu_beta,
                                       logsd_beta,
                                       1,
                                       1,
                                       re_cor,
                                       re_sigmas,
                                       factor_cor,
                                       zeta,
                                       X_loc,
                                       X_sca,
                                       X_bet,
                                       FALSE
                                       )

sim_data <- d$df
names(sim_data)[1:J] <- c(paste0("A_", 1:(J/2)), paste0("N_", 1:(J/2)))
names(sim_data)[14:16] <- c("x1", "baseline", "x2")
sim_data[c(17:18)] <- NULL

sim_data <- sim_data[, c("subject","baseline","x1","x2",paste0("A_",1:(J/2)), paste0("N_",1:(J/2)))]
