##' Ranef method for lmmelsm objects.
##'
##' Extracts the random effects from the lmmelsm object.
##' Note that this is different from the random \emph{coefficients}.
##' E.g., if $$\beta_{0i} = \beta_0 + u_{0i} $$, then \code{coef} extracts $\beta_{0i}$ and \code{ranef} extracts $u_{0i}$.
##' @title Extract random effects from lmmelsm objects.
##' @return List of ranef summaries, or samples (if summarize = FALSE).
##' @author Stephen R. Martin
##' @importFrom nlme ranef
##' @export ranef
##' @param object lmmelsm object.
##' @param prob Numeric (Default: .95). Amount of probability mass contained in the credible interval.
##' @param summarize Logical (Default: TRUE). Whether to return posterior summaries (TRUE) or MCMC samples (FALSE).
##' @export
ranef.lmmelsm <- function(object, prob = .95, summarize = TRUE) {
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
        out <- summary(x, prob = prob)$summary[c("mu_random","logsd_random","mu_beta_random","logsd_beta_random")]
        names(out) <- c("location","scale","location_slope","scale_slope")
        return(out)
    }

    mu_random <- as.matrix(x$fit, pars = "mu_random")
    logsd_random <- as.matrix(x$fit, pars = "logsd_random")
    mu_beta_random <- as.matrix(x$fit, pars = "mu_beta_random")
    logsd_beta_random <- as.matrix(x$fit, pars = "logsd_beta_random")

    out <- list(location = mu_random,
                scale = logsd_random,
                location_slope = mu_beta_random,
                scale_slope = logsd_beta_random
                )

    return(out)
}
##' Coef method for lmmelsm objects.
##'
##' Extracts all group-specific coefficients from lmmelsm object.
##' Note that this is different from \code{\link{ranef}}.
##' Whereas \code{ranef} extracts the zero-centered random effects, \code{coef} extracts the group-specific effects, defined as the sum of the fixed effect and random effect.
##' @title Extract random coefficients for each group.
##' @param object lmmelsm object.
##' @param prob Numeric (Default: .95). Amount of probability mass contained in the credible interval.
##' @param summarize Logical (Default: TRUE). Whether to return posterior summaries (TRUE) or MCMC samples (FALSE).
##' @return List of summaries (if \code{summarize} is TRUE), or list of MCMC samples.
##' @author Stephen R Martin
##' @export
##' @importFrom nlme coef
##' @export coef
##' @aliases coef
##' @method coef lmmelsm
coef.lmmelsm <- function(object, prob = .95, summarize = TRUE) {
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

    # Get samples
    ## Fixed effects
    ## Mu and logsd fixef = 0
    mu_coef <- as.matrix(x$fit, pars = "mu_random")
    logsd_coef <- as.matrix(x$fit, pars = "logsd_random")

    mu_fixed_inds <- as.matrix(expand.grid(P_random_ind, 1:F))
    mu_fixed_inds <- paste0("[",mu_fixed_inds[,1],",",mu_fixed_inds[,2],"]")

    logsd_fixed_inds <- as.matrix(expand.grid(Q_random_ind, 1:F))
    logsd_fixed_inds <- paste0("[",logsd_fixed_inds[,1],",",logsd_fixed_inds[,2],"]")

    mu_beta <- as.matrix(x$fit, pars = paste0("mu_beta", mu_fixed_inds))
    logsd_beta <- as.matrix(x$fit, pars = paste0("logsd_beta", logsd_fixed_inds))

    ## Random effects
    mu_beta_random <- as.matrix(x$fit, pars = "mu_beta_random")
    logsd_beta_random <- as.matrix(x$fit, pars = "logsd_beta_random")

    ## Restructure for vectorized addition
    S <- nrow(mu_coef) # samples

    ### As array
    mu_beta_random_arr <- array(mu_beta_random, dim = c(S, GS$K, P_random, F))
    logsd_beta_random_arr <- array(logsd_beta_random, dim = c(S, GS$K, Q_random, F))

    ### Reorder to be [Samples, predictor, factor, group]
    mu_beta_coef <- aperm(mu_beta_random_arr, c(1, 3, 4, 2))
    logsd_beta_coef <- aperm(logsd_beta_random_arr, c(1, 3, 4, 2))

    ### Add to this a recycled fixed effect array
    mu_beta_coef <- mu_beta_coef + array(mu_beta, c(S, P_random, F, GS$K))
    logsd_beta_coef <- logsd_beta_coef + array(logsd_beta, c(S, Q_random, F, GS$K))

    ### Wind it back down to be [S, Coefs] to pass to .summarize, or to return.
    mu_beta_coef <- array(mu_beta_coef, c(S, P_random*F*GS$K))
    logsd_beta_coef <- array(logsd_beta_coef, c(S, Q_random*F*GS$K))

    col_renames_mu <- expand.grid(P_random_ind, 1:F, 1:GS$K)
    col_renames_logsd <- expand.grid(Q_random_ind, 1:F, 1:GS$K)
    colnames(mu_beta_coef) <- paste0("mu_beta[", col_renames_mu[,3], ",", col_renames_mu[,1], ",", col_renames_mu[,2], "]")
    colnames(logsd_beta_coef) <- paste0("logsd_beta[", col_renames_logsd[,3], ",", col_renames_logsd[,1], ",", col_renames_logsd[,2], "]")

    ### Summarize if needed.
    if(summarize) {
        # Summarize
        mu_coef <- .summarize(mu_coef, NULL, prob = prob)
        logsd_coef <- .summarize(logsd_coef, NULL, prob = prob)
        mu_beta_coef <- .summarize(mu_beta_coef, NULL, prob = prob)
        logsd_beta_coef <- .summarize(logsd_beta_coef, NULL, prob = prob)
        # Tidy up
        mu_coef <- .tidy_summary(mu_coef, c(GS$name, "Factor"), GS$map$label, fnames)
        logsd_coef <- .tidy_summary(logsd_coef, c(GS$name, "Factor"), GS$map$label, fnames)
        mu_beta_coef <- .tidy_summary(mu_beta_coef, c(GS$name, "Predictor", "Factor"), GS$map$label, pnames$location, fnames)
        logsd_beta_coef <- .tidy_summary(logsd_beta_coef, c(GS$name, "Predictor", "Factor"), GS$map$label, pnames$scale, fnames)
        # Package up [Remember to rearrange later]
        out <- list(location = mu_coef,
                    scale = logsd_coef,
                    location_slope = mu_beta_coef,
                    scale_slope = logsd_beta_coef
                    )
        return(out)
    } else {
        out <- list(location = mu_coef,
                    scale = logsd_coef,
                    location_slope = mu_beta_coef,
                    scale_slope = logsd_beta_coef
                    )
        return(out)
    }
}
