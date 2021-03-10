predict.lmmelsm <- function(object, newdata = NULL, prob = .95, summarize = TRUE, what = c("latent", "indicators"), include_error = TRUE, ...) {
    x <- object

    # If no newdata, use old data.
    newdata <- newdata %IfNull% x$data

    # Make numeric and NA group columns
    newdata <- .add_group_codings(x$meta$group_spec, newdata)

    # Group numeric
    group <- newdata[,x$meta$group_spec$name]

    # Get formulas
    latent <- object$meta$latent
    plist <- object$meta$pred_spec$plist
    mlist <- object$meta$indicator_spec$mlist

    # Get newdata predictor and indicator specifications
    pred_pred_spec <- .parse_formula.predictor(plist, newdata, seq_len(nrow(newdata)))
    if(latent) {
        pred_indicator_spec <- .parse_formula.indicators(mlist, newdata)
    } else {
        pred_indicator_spec <- .parse_formula.observed(mlist, newdata)
    }

    # Get samples
    zeta <- .extract_transform_to_list(object$fit, par = "zeta")
    random_sigma <- .extract_transform_to_list(object$fit, par = "mu_logsd_betas_random_sigma")
    random_cor <- .extract_transform_to_list(object$fit, par = "Omega_mean_logsd")

    mu_beta <- .extract_transform_to_list(object$fit, par = "mu_beta")
    logsd_beta <- .extract_transform_to_list(object$fit, par = "logsd_beta")
    factor_cor <- .extract_transform_to_list(object$fit, par = "Omega_eta")

    mu_random <- .extract_transform_to_list(object$fit, par = "mu_random")
    logsd_random <- .extract_transform_to_list(object$fit, par = "logsd_random")
    mu_beta_random <- .extract_transform_to_list(object$fit, par = "mu_beta_random")
    logsd_beta_random <- .extract_transform_to_list(object$fit, par = "logsd_beta_random")

    # Set up predictions

    pred_list <- list()
    N <- nrow(newdata)
    S <- length(mu_random)

    for(n in seq_len(N_pred)) { # Each row of newdata
        # Between-group variance
        ## pred_pred_spec$x_bet[n,, drop = FALSE] %*% 
        between <- random_sigma

    }
}

# Converts fitted random effects to uncorrelated z
# Returns list of random_z's; one list per sample.
# Does not split into different REs.
.random_to_z <- function(fit) {
    ranef_pars <- c("mu_random","logsd_random")
    if(fit$meta$pred_spec$P_random > 0) ranef_pars <- c(ranef_pars, "mu_beta_random")
    if(fit$meta$pred_spec$Q_random > 0) ranef_pars <- c(ranef_pars, "logsd_beta_random")
    K <- fit$meta$group_spec$K
    F <- fit$meta$indicator_spec$F
    N_ranefs <- F + F + F*fit$meta$pred_spec$P_random + F*fit$meta$pred_spec$Q_random

    # ranef_pars-length List of Samples-length Lists of KxRandom matrices
    ranefs_list <- lapply(ranef_pars, function(par) {.extract_transform_to_list(fit$fit, par)})
    # List of Samples-length KxRandom matrices
    ranefs <- do.call(.list_zip, ranefs_list)

    random_sigma <- .extract_transform_to_list(fit$fit, par = "mu_logsd_betas_random_sigma")
    random_cor <- .extract_transform_to_list(fit$fit, par = "Omega_mean_logsd")
    random_cor_U <- lapply(random_cor, function(R){chol(R)})

    S <- length(random_sigma)

    R <- fit$meta$pred_spec$R
    if(R > 0) { # If between model
        zeta <- .extract_transform_to_list(fit$fit, par = "zeta")
        x_bet_l2 <- fit$stan_data$x_bet[match(1:K, fit$stan_data$group),]
        between_pred <- lapply(zeta, function(s) {x_bet_l2 %*% s})
    } else { # not between model; fill with zeroes
        between_pred <- lapply(seq_len(S), function(s){matrix(0, K, 2*F)})
    }

    # For each row in ranefs...
    for(k in 1:K) {
        # For each sample
        for(s in seq_len(S)) {
            # RE' = 0 + z' * L' * re_sigma
            # z = RE' * L'^-1 * re_sigma^-1
            re_sigma_k_s <- random_sigma[[s]]
            re_sigma_k_s[1:(2*F)] <- re_sigma_k_s[1:(2*F)] * exp(between_pred[[s]][k,])
            ranefs[[s]][k,] <- ranefs[[s]][k,,drop = FALSE] %*% solve(random_cor_U[[s]]) %*% (diag(1 / re_sigma_k_s))
        }
    }
    
    return(ranefs)
}

.row_multiply_list_mats <- function(r, mats) {
    lapply(mats, function(m) {
        r %*% m
    })
}

##' @title Adds group codings for predictive df.
##' @param gs group spec
##' @param newdata data frame
##' @return data frame, with numerics where group matches, and NAs otherwise.
##' @author Stephen Martin
##' @keywords internal
.add_group_codings <- function(gs, newdata) {
    if(!(gs$name %in% colnames(newdata))) {
        newdata[, gs$name] <- NA
    }
    newdata[, gs$name] <- gs$numeric[match(newdata[, gs$name], gs$data)]
    return(newdata)
}

##' Extract model fitted variates.
##'
##' Extracts model fitted variates. When a latent MMELSM, these are the latent score expectations and log standard deviations. When an observed MMELSM, these are the observed score expectations and log standard deviations.
##' @title Extracted model fitted variates.
##' @param object lmmelsm object.
##' @inheritParams ranef.lmmelsm
##' @param ... Not used.
##' @return List of eta and eta_logsd. If summarize is \code{TRUE}, then the list contains formatted summary tables. If \code{FALSE}, then the list contains MCMC samples for all variables.
##' @author Stephen Martin
fitted.lmmelsm <- function(object, prob = .95, summarize = TRUE, ...) {
    pars <- c("eta", "eta_logsd")
    facVar <- ifelse(object$meta$latent, "factor", "variable")
    fnames <- unlist(object$meta$indicator_spec$fname)

    eta <- as.matrix(object, pars = pars[1])
    eta_logsd <- as.matrix(object, pars = pars[2])

    if(summarize) {
        eta <- .summarize(eta, prob = prob)
        eta_logsd <- .summarize(eta_logsd, prob = prob)

        eta <- .tidy_summary(eta, c("observation", facVar), seq_len(nrow(object$data)), fnames)
        eta_logsd <- .tidy_summary(eta_logsd, c("observation", facVar), seq_len(nrow(object$data)), fnames)

        eta <- .summary_rearrange(eta, c("observation", facVar))
        eta_logsd <- .summary_rearrange(eta_logsd, c("observation", facVar))

        return(nlist(eta, eta_logsd))
    } else {
        return(nlist(eta, eta_logsd))
    }
}
