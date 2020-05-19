##' @title Simulate univariate fixed effects predictive LM-MELSM
##' @param n Integer. Number of observations per group.
##' @param K Integer. Number of groups.
##' @param lambda Vector. Loadings.
##' @param resid Vector. Residual SD.
##' @param nu Vector. Intercepts.
##' @param mu_beta Vector. \code{P} Predictive coefficients for latent location.
##' @param logsd_beta Vector. \code{Q} predictive coefficients for latent IIV.
##' @param mu_logsd_cor Matrix. 2x2 correlation matrix between location and scale.
##' @param mu_logsd_sigma Vector (Length 2). SDs for mu and logsd RE's.
##' @param L2_pred_only Logical (Default: FALSE). If generated predictors should be level-2 only.
##' @param X_loc Matrix (Optional). If supplied, will be used as the location predictors.
##' @param X_sca Matrix (Optional). If supplied, will be used as the scale predictors.
##' @return List of data and params.
##' @author Stephen R. Martin
##' @keywords internal
##' @importFrom MASS mvrnorm
simulate.uni.fe <- function(n, K, lambda, resid, nu, mu_beta = numeric(), logsd_beta = numeric(), mu_logsd_cor, mu_logsd_sigma, L2_pred_only = FALSE, X_loc = NULL, X_sca = NULL) {
    P <- length(mu_beta)
    Q <- length(logsd_beta)
    F <- 1
    N <- n * K

    if(P > 0) X_loc <- X_loc %IfNull% .simulate.X(n, K, P, L2_pred_only)
    if(Q > 0) X_sca <- X_sca %IfNull% .simulate.X(n, K, Q, L2_pred_only)
}
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
        X_L2 <- mvrnorm(K, Sigma = diag(1, P, P))
        X <- t(apply(X_L2, 1, function(x){
            rep(x, n)
        }))
    } else {
        X <- mvrnorm(n * K, Sigma = diag(1, P, P))
    }
    return(X)
}
