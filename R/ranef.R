##' Ranef method for lmmelsm objects.
##'
##' Extracts the random effects from the lmmelsm object.
##' Note that this is different from the random \emph{coefficients}.
##' E.g., if \eqn{\beta_{0i} = \beta_0 + u_{0i}}, then \code{coef} extracts \eqn{\beta_{0i}} and \code{ranef} extracts \eqn{u_{0i}}.
##' @title Extract random effects.
##' @return List of ranef summaries (random_mu_intercept, random_logsd_intercept, random_mu_coef, and random_logsd_coef), or samples (if summarize = FALSE).
##' @author Stephen R. Martin
##' @importFrom nlme ranef
##' @export ranef
##' @aliases ranef
##' @param object lmmelsm object.
##' @param prob Numeric (Default: .95). Amount of probability mass contained in the credible interval.
##' @param summarize Logical (Default: TRUE). Whether to return posterior summaries (TRUE) or MCMC samples (FALSE).
##' @param ... Not used.
##' @export
ranef.lmmelsm <- function(object, prob = .95, summarize = TRUE, ...) {
    x <- object

    IS <- x$meta$indicator_spec
    PS <- x$meta$pred_spec
    GS <- x$meta$group_spec

    # Meta-data
    pnames <- PS$pname
    fnames <- unlist(IS$fname)
    F <- IS$F
    P_random <- PS$P_random
    Q_random <- PS$Q_random
    P_random_ind <- PS$P_random_ind
    Q_random_ind <- PS$Q_random_ind
    re_total <- 2 * F + F * P_random + F * Q_random

    if(summarize) {
        out <- summary(x, prob = prob)$summary[c("random_mu_intercept",
                                                 "random_logsd_intercept",
                                                 "random_mu_coef",
                                                 "random_logsd_coef")]
        return(out)
    }

    mu_random <- as.matrix(x$fit, pars = "mu_random")
    logsd_random <- as.matrix(x$fit, pars = "logsd_random")

    mu_beta_random <- logsd_beta_random <- NA

    if(P_random > 0) {
        mu_beta_random <- as.matrix(x$fit, pars = "mu_beta_random")
    }
    if(Q_random > 0) {
        logsd_beta_random <- as.matrix(x$fit, pars = "logsd_beta_random")
    }

    out <- list(random_mu_intercept = mu_random,
                random_logsd_intercept = logsd_random,
                random_mu_coef = mu_beta_random,
                random_logsd_coef = logsd_beta_random
                )

    return(out)
}
##' Coef method for lmmelsm objects.
##'
##' Extracts all group-specific coefficients from lmmelsm object.
##' Note that this is different from \code{\link{ranef}}.
##' Whereas \code{ranef} extracts the zero-centered random effects, \code{coef} extracts the group-specific effects, defined as the sum of the fixed effect and random effect.
##' @title Extract group-specific coefficients.
##' @param object lmmelsm object.
##' @param prob Numeric (Default: .95). Amount of probability mass contained in the credible interval.
##' @param summarize Logical (Default: TRUE). Whether to return posterior summaries (TRUE) or MCMC samples (FALSE).
##' @param ... Not used.
##' @return List of summaries (if \code{summarize} is TRUE), or list of MCMC samples.
##' @author Stephen R Martin
##' @export
##' @importFrom stats coef
##' @aliases coef
##' @method coef lmmelsm
coef.lmmelsm <- function(object, prob = .95, summarize = TRUE, ...) {
    x <- object

    IS <- x$meta$indicator_spec
    PS <- x$meta$pred_spec
    GS <- x$meta$group_spec

    # Meta-data
    pnames <- PS$pname
    fnames <- unlist(IS$fname)
    F <- IS$F
    P_random <- PS$P_random
    Q_random <- PS$Q_random
    P_random_ind <- PS$P_random_ind
    Q_random_ind <- PS$Q_random_ind
    re_total <- 2 * F + F * P_random + F * Q_random
    latent <- x$meta$latent

    # Get samples
    ## Fixed effects
    ## Intercepts
    ## Mu and logsd fixef = 0 WHEN LATENT is true
    mu_coef <- as.matrix(x$fit, pars = "mu_random")
    logsd_coef <- as.matrix(x$fit, pars = "logsd_random")
    if(!latent) {
        S <- nrow(mu_coef)
        nu <- as.matrix(x$fit, pars = "nu")
        sigma <- as.matrix(x$fit, pars = "sigma") # sigma = log_sigma
        rep_cols <- rep(1:F, each = GS$K) # F=J=number of nu/sigma
        mu_coef <- mu_coef + nu[, rep_cols]
        logsd_coef <- logsd_coef + sigma[, rep_cols]
    }
    

    S <- nrow(mu_coef) # samples

    mu_beta_coef <- logsd_beta_coef <- NA
    if(P_random > 0) { # Random location coefficients
        mu_fixed_inds <- as.matrix(expand.grid(P_random_ind, 1:F))
        mu_fixed_inds <- paste0("[",mu_fixed_inds[,1],",",mu_fixed_inds[,2],"]")

        ## Fixed
        mu_beta <- as.matrix(x$fit, pars = paste0("mu_beta", mu_fixed_inds))
        ## Random
        mu_beta_random <- as.matrix(x$fit, pars = "mu_beta_random")

        ## Restructure for vectorized addition
        ### As array
        mu_beta_random_arr <- array(mu_beta_random, dim = c(S, GS$K, P_random, F))

        ### Reorder to be [Samples, predictor, factor, group]
        mu_beta_coef <- aperm(mu_beta_random_arr, c(1, 3, 4, 2))

        ### Add to this a recycled fixed effect array
        mu_beta_coef <- mu_beta_coef + array(mu_beta, c(S, P_random, F, GS$K))

        ### Repermute to be more like ranef order (Column major)
        mu_beta_coef <- aperm(mu_beta_coef, c(1, 4, 2, 3))

        ### Wind it back down to be [S, Coefs] to pass to .summarize, or to return.
        mu_beta_coef <- array(mu_beta_coef, c(S, P_random*F*GS$K))

        ## Regenerate matrix labels
        col_renames_mu <- expand.grid(1:GS$K, P_random_ind, 1:F)
        colnames(mu_beta_coef) <- paste0("mu_beta[", col_renames_mu[,1], ",", col_renames_mu[,2], ",", col_renames_mu[,3], "]")
    }

    if(Q_random > 0) { # Random scale coefficients
        logsd_fixed_inds <- as.matrix(expand.grid(Q_random_ind, 1:F))
        logsd_fixed_inds <- paste0("[",logsd_fixed_inds[,1],",",logsd_fixed_inds[,2],"]")

        ## Fixed
        logsd_beta <- as.matrix(x$fit, pars = paste0("logsd_beta", logsd_fixed_inds))
        ## Random
        logsd_beta_random <- as.matrix(x$fit, pars = "logsd_beta_random")

        ## Restructure for vectorized addition
        ### As array
        logsd_beta_random_arr <- array(logsd_beta_random, dim = c(S, GS$K, Q_random, F))

        ### Reorder to be [Samples, predictor, factor, group]
        logsd_beta_coef <- aperm(logsd_beta_random_arr, c(1, 3, 4, 2))

        ### Add to this a recycled fixed effect array
        logsd_beta_coef <- logsd_beta_coef + array(logsd_beta, c(S, Q_random, F, GS$K))

        ### Repermute to be more like ranef order (Column major)
        logsd_beta_coef <- aperm(logsd_beta_coef, c(1, 4, 2, 3))

        ### Wind it back down to be [S, Coefs] to pass to .summarize, or to return.
        logsd_beta_coef <- array(logsd_beta_coef, c(S, Q_random*F*GS$K))

        ## Regenerate matrix labels
        col_renames_logsd <- expand.grid(1:GS$K, Q_random_ind, 1:F)
        colnames(logsd_beta_coef) <- paste0("logsd_beta[", col_renames_logsd[,1], ",", col_renames_logsd[,2], ",", col_renames_logsd[,3], "]")
    }


    ### Summarize if needed.
    if(summarize) {
        # Summarize
        mu_coef <- .summarize(mu_coef, NULL, prob = prob)
        mu_coef <- .tidy_summary(mu_coef, c(GS$name, "factor"), GS$map$label, fnames)

        logsd_coef <- .summarize(logsd_coef, NULL, prob = prob)
        logsd_coef <- .tidy_summary(logsd_coef, c(GS$name, "factor"), GS$map$label, fnames)
        if(P_random > 0) {
            mu_beta_coef <- .summarize(mu_beta_coef, NULL, prob = prob)
            mu_beta_coef <- .tidy_summary(mu_beta_coef, c(GS$name, "predictor", "factor"), GS$map$label, pnames$location, fnames)
        }
        if(Q_random > 0) {
            logsd_beta_coef <- .summarize(logsd_beta_coef, NULL, prob = prob)
            logsd_beta_coef <- .tidy_summary(logsd_beta_coef, c(GS$name, "predictor", "factor"), GS$map$label, pnames$scale, fnames)
        }

        # Package up
        out <- list(mu_intercept = mu_coef,
                    logsd_intercept = logsd_coef,
                    mu_coef = mu_beta_coef,
                    logsd_coef = logsd_beta_coef
                    )
        return(out)
    } else {
        out <- list(mu_intercept = mu_coef,
                    logsd_intercept = logsd_coef,
                    mu_coef = mu_beta_coef,
                    logsd_coef = logsd_beta_coef
                    )
        return(out)
    }
}
