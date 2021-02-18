
##' @title Simulate covariates without correlation.
##' @param n Number of observations per k unit.
##' @param K Number of k units (Level 2 units).
##' @param P Number of covariates.
##' @param L2_pred_only Logical (Default: FALSE). If generated predictors should be level-2 only.
##' @return Matrix.
##' @author Stephen R. Martin
##' @keywords internal
##' @importFrom MASS mvrnorm
.simulate.X <- function(n, K, P, L2_pred_only) {
    if(L2_pred_only) {
        X_L2 <- mvrnorm(K, mu = rep(0, P), Sigma = diag(1, P, P))
        X <- X_L2[rep(1:K, each = n),, drop = FALSE]
    } else {
        X <- mvrnorm(n * K, mu = rep(0, P), Sigma = diag(1, P, P))
    }
    return(X)
}

##' @title Simulate data from latent uni/multidimensional MELSM
##' @param n Integer. Number of repeated observations per group.
##' @param K Integer. Number of groups.
##' @param lambda Matrix (FxJ). Loading matrix.
##' @param resid Numeric vector (J). Residual SDs.
##' @param nu Numeric vector (J). Intercepts.
##' @param mu_beta Matrix (PxF). Location coefficient matrix.
##' @param logsd_beta Matrix (QxF). Scale coefficient matrix.
##' @param P_random_ind Integer vector (P_random). Which location predictors have random slopes.
##' @param Q_random_ind Integer vector (Q_random). Which scale predictors have random slopes.
##' @param mu_logsd_betas_cor Matrix (Symmetric, SPD; F*2 + P_random*F + Q_random*F). Correlation matrix of random effects (slopes and intercepts, for location and scale models).
##' @param mu_logsd_betas_sigma Numeric vector (Positive; F*2 + P_random*F + Q_random*F). RE SDs (intercepts on exponentiated scale, if zeta is specified).
##' @param epsilon_cor Matrix (Symmetric, SPD; F). Stochastic error term correlation between factors.
##' @param zeta Matrix (`Rx[F*2 + P_random*F + Q_random*F]`). Coefficient matrix for predicting RE SDs.
##' @param X_loc Matrix (Optional; NxP). Location design matrix.
##' @param X_sca Matrix (Optional; NxQ). Scale design matrix.
##' @param X_bet Matrix (Optional; NxR). Between-SD design matrix.
##' @param L2_pred_only Logical. Whether predictors should be group-level (TRUE) or observation level (FALSE).
##' @return List of params (list), data (list), and df (data.frame).
##' @author Stephen R. Martin
##' @keywords internal
simulate_lmmelsm <- function(n,
                             K,
                             lambda,
                             resid,
                             nu,
                             mu_beta = NULL,
                             logsd_beta = NULL,
                             P_random_ind = NULL,
                             Q_random_ind = NULL,
                             mu_logsd_betas_cor,
                             mu_logsd_betas_sigma,
                             epsilon_cor,
                             zeta = NULL,
                             X_loc = NULL,
                             X_sca = NULL,
                             X_bet = NULL,
                             L2_pred_only = FALSE
                             ) {
    # Restructure; If not matrix, then make matrix, assume unidimensional model.
    if(!is.matrix(lambda)) lambda <- matrix(lambda, nrow = 1)
    if(!is.matrix(mu_beta) & !is.null(mu_beta)) mu_beta <- matrix(mu_beta, ncol = 1)
    if(!is.matrix(logsd_beta) & !is.null(mu_beta)) logsd_beta <- matrix(logsd_beta, ncol = 1)

    # Dimensions
    J <- ncol(lambda) 
    F <- nrow(lambda) 
    P <- nrow(mu_beta) %IfNull% 0
    Q <- nrow(logsd_beta) %IfNull% 0
    R <- nrow(zeta) %IfNull% 0
    P_random <- length(P_random_ind) %IfNull% 0
    Q_random <- length(Q_random_ind) %IfNull% 0
    P_random_ind <- P_random_ind %IfNull% array(0, dim = c(P_random))
    Q_random_ind <- Q_random_ind %IfNull% array(0, dim = c(Q_random))
    N <- n * K

    # Generate groups, if not provided
    ## Removing this for now; requires more care, since n may not be balanced across groups.
    ## group <- group %IfNull% rep(1:K, each = n)
    group <- rep(1:K, each = n)

    # Predictor data
    if(P > 0) X_loc <- X_loc %IfNull% .simulate.X(n, K, P, L2_pred_only)
    if(Q > 0) X_sca <- X_sca %IfNull% .simulate.X(n, K, Q, L2_pred_only)
    if(R > 0) {
        X_bet <- X_bet %IfNull% .simulate.X(n, K, R, L2_pred_only = TRUE)
        X_bet_L2 <- X_bet[match(1:K, group), ]
    }

    # Generate latent values
    eta_mu <- matrix(0, nrow = N, ncol = F)
    eta_logsd <- matrix(0, nrow = N, ncol = F)

    ## Fixed effects
    if(P > 0) eta_mu <- eta_mu + X_loc %*% mu_beta
    if(Q > 0) eta_logsd <- eta_logsd + X_sca %*% logsd_beta

    ## Random effects
    num_rands <- 2*F + P_random*F + Q_random*F
    mu_logsd_betas_re <- matrix(0, nrow = K, ncol = num_rands)
    if(num_rands != length(mu_logsd_betas_sigma)) stop("Incorrect number of mu_logsd_betas_sigma. Should contain ", num_rands, "elements")

    #### Between variance model
    mu_logsd_betas_sigmas <- matrix(mu_logsd_betas_sigma, nrow = K, ncol = num_rands, byrow = TRUE)
    if(R > 0) {
        if(ncol(zeta) != num_rands) stop("Mismatch between zeta dimensions and total number of random effects.")
        mu_logsd_betas_sigmas <- mu_logsd_betas_sigmas * exp(X_bet_L2 %*% zeta)
    }
    #### Add RE contributions
    for(k in 1:K) {
        mu_logsd_betas_re[k,] <- mvrnorm(1, rep(0, num_rands), Sigma = diag(mu_logsd_betas_sigmas[k,], num_rands, num_rands) %*% mu_logsd_betas_cor %*% diag(mu_logsd_betas_sigmas[k,], num_rands, num_rands))
    }

    eta_mu <- eta_mu + mu_logsd_betas_re[group, 1:F]
    eta_logsd <- eta_logsd + mu_logsd_betas_re[group, (F + 1):(2 * F)]

    if(P_random > 0) {
        mu_beta_random <- array(t(mu_logsd_betas_re[, (2*F + 1):(2*F + P_random*F)]), dim = c(P_random, F, K))
        for(n in 1:N) {
            eta_mu[n,] <- eta_mu[n,] + X_loc[n, P_random_ind, drop = FALSE] %*% .array_extract(mu_beta_random, group[n])
        }
    } else {mu_beta_random <- numeric()}

    if(Q_random > 0) {
        logsd_beta_random <- array(t(mu_logsd_betas_re[, (2*F + P_random*F + 1):(2*F + P_random*F + Q_random*F)]), dim = c(Q_random, F, K))
        for(n in 1:N) {
            eta_logsd[n,] <- eta_logsd[n,] + X_sca[n, Q_random_ind, drop = FALSE] %*% .array_extract(logsd_beta_random, group[n])
        }
    } else {logsd_beta_random <- numeric()}

    # Generate Etas
    eta <- eta_mu
    for(n in 1:N) {
        eta[n,] <- eta[n,] + mvrnorm(1, rep(0, F), Sigma = diag(exp(eta_logsd[n,]), F, F) %*% epsilon_cor %*% diag(exp(eta_logsd[n,]), F, F))
    }

    # Measurement model
    Y <- matrix(nu, nrow = N, ncol = J, byrow = TRUE) + eta %*% lambda + matrix(rnorm(N*J, 0, resid), nrow = N, ncol = J, byrow = TRUE)

    # Construct output

    df <- as.data.frame(Y)
    colnames(df) <- paste0("obs_", 1:J)
    df$subject <- group

    if(P > 0) {
        colnames(X_loc) <- paste0("loc_", 1:P)
        df <- cbind(df, X_loc)
    } else {
        X_loc <- array(0, dim =c(N, 0))
    }
    if(Q > 0) {
        colnames(X_sca) <- paste0("sca_", 1:Q)
        df <- cbind(df, X_sca)
    } else {
        X_sca <- array(0, dim =c(N, 0))
    }

    if(R > 0) {
        colnames(X_bet) <- paste0("bet_", 1:R)
        df <- cbind(df, X_bet)
    } else {
        X_bet <- array(0, dim =c(N, 0))
    }

    J_f <- apply(lambda, 1, function(r) {
        sum(r != 0)
    })
    F_ind <- matrix(0, nrow = F, ncol = J)
    for(r in 1:F) {
        F_ind[r,1:(J_f[r])] <- which(lambda[r,] != 0)
    }

    out <- list(params = nlist(
                    N, J, F, K, P, Q, n, R,
                    P_random, Q_random,
                    P_random_ind, Q_random_ind,
                    lambda, resid, nu,
                    mu_beta, logsd_beta, zeta,
                    mu_logsd_betas_cor, mu_logsd_betas_sigma,
                    epsilon_cor,
                    eta, eta_logsd, mu_logsd_betas_re
                ),
                data = nlist(
                    N, J, F, K, P, Q, R,
                    P_random, Q_random,
                    P_random_ind, Q_random_ind,
                    x_loc = X_loc,
                    x_sca = X_sca,
                    x_bet = X_bet,
                    y = Y,
                    group,
                    J_f,
                    F_ind,
                    L2_pred_only = L2_pred_only,
                    prior_only = FALSE
                ),
                df = df)

    return(out)

}
