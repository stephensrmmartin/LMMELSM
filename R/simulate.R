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
##' @return List of data and params.
##' @author Stephen R. Martin
##' @keywords internal
simulate.uni.fe <- function(n, K, lambda, resid, nu, mu_beta, logsd_beta, mu_logsd_cor, mu_logsd_sigma) {
    # TODO: Implement this
    
}
