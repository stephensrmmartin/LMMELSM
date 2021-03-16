##' Generates posterior predictions from fitted LMMELSM object.
##'
##' If the grouping variable is missing, or contains NAs, then new random effects are generated from the posterior random effect distribution.
##' Where the grouping variables are not missing, the posterior standardized, orthogonalized random effects are obtained from the fitted model, and used as a basis for predicted random effects.
##' Because the standardized, orthogonalized random effects are used, one can include different between-group variance values to predict new RE variances, and therefore a different random effect value.
##' That is, the random effect, conditional on new between-group variance model covariates, is equal to: \deqn{u_i = z_i U D_i}, where \eqn{D_i} is the predicted between-group random effect SD, U is the upper cholesky factorization of the random effect correlations, and \eqn{z_i} is the standardized, orthogonalized random effect for group i.
##' @title Predict method for lmmelsm objects.
##' @param object lmmelsm object.
##' @param newdata Data.frame (Default: NULL). If NULL, uses original data.frame.
##' @inheritParams ranef.lmmelsm
##' @param include_error Logical (Default: TRUE). If TRUE, then include stochastic realizations in outcome variables.
##' @param ... Not used.
##' @return List. If summarize is TRUE, then a list of outcome (eta, eta_logsd) and indicator (y) posterior predictive distribution summaries. If FALSE, an N-length list of lists of outcome and indicator MCMC samples.
##' @author Stephen Martin
##' @export
predict.lmmelsm <- function(object, newdata = NULL, prob = .95, summarize = TRUE,  include_error = TRUE, ...) {
    x <- object
    include_error <- include_error

    facVar <- ifelse(has_latent(x), "factor", "variable")

    # If no newdata, use old data.
    newdata <- newdata %IfNull% x$data

    # Make numeric and NA group columns
    newdata <- .add_group_codings(x$meta$group_spec, newdata)

    # Group numeric
    group <- newdata[,x$meta$group_spec$name]

    # Get formulas
    latent <- has_latent(x)
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
    random_sigma <- .extract_transform_to_list(object$fit, par = "mu_logsd_betas_random_sigma")
    random_cor <- .extract_transform_to_list(object$fit, par = "Omega_mean_logsd")
    random_cor_U <- lapply(random_cor, function(s){chol(s)})
    if(has_between(object)) {
        zeta <- .extract_transform_to_list(object$fit, par = "zeta")
    } else {
        zeta <- lapply(seq_len(length(random_sigma)), function(s){0})
    }

    if(has_location(x)) {
        mu_beta <- .extract_transform_to_list(object$fit, par = "mu_beta")
        if(has_random_location(x)) {
            P_random_ind <- get_P_random_ind(x)
        }
    }
    if(has_scale(x)) {
        logsd_beta <- .extract_transform_to_list(object$fit, par = "logsd_beta")
        if(has_random_scale(x)) {
            Q_random_ind <- get_Q_random_ind(x)
        }
    }
    factor_cor <- .extract_transform_to_list(object$fit, par = "Omega_eta")

    random_z <- .random_to_z(x)

    if(has_latent(x)) {
        lambda <- .extract_transform_to_list(object$fit, par = "lambda")
        sigma <- .extract_transform_to_list(object$fit, par = "sigma")
    }
    nu <- .extract_transform_to_list(object$fit, par = "nu")

    # Set up predictions

    pred_list <- list()
    N <- nrow(newdata)
    J <- get_J(x)
    S <- get_S(x)
    re_total <- get_number_re(x)
    re_indices <- get_re_indices(x)
    re_intercept_indices <- unlist(re_indices[1:2])
    P_random <- get_P_random(x)
    Q_random <- get_Q_random(x)
    P_random_ind <- get_P_random_ind(x)
    Q_random_ind <- get_Q_random_ind(x)
    F <- get_F(x)

    for(n in seq_len(N)) { # Each row of newdata
        # Between-group variance
        pred_bet <- random_sigma
        if(has_between(x)) {
            for(s in seq_len(S)) {
                pred_bet[[s]][re_intercept_indices] <- t(pred_bet[[s]][re_intercept_indices]) * exp(pred_pred_spec$x_bet[n,,drop=FALSE] %*% zeta[[s]])
            }
        }

        # Random effects
        if(is.na(group[n])) { # If grouping variable missing; random normal variates
            pred_random_z <- lapply(seq_len(S), function(s) {
                matrix(rnorm(re_total), nrow = 1)
            })
        } else { # Use samples from group[n]
            pred_random_z <- lapply(random_z, function(s) {
                s[group[n], , drop = FALSE]
            })
        }
        pred_random <- lapply(seq_len(S), function(s) {
            pred_random_z[[s]] %*% random_cor_U[[s]] %*% diag(pred_bet[[s]])
        })

        # Eta, log SD(Eta)
        ## Intercepts
        eta <- lapply(pred_random, function(s) {
            s[, re_indices$mu_random, drop = FALSE]
        })
        eta_logsd <- lapply(pred_random, function(s) {
            s[, re_indices$logsd_random, drop = FALSE]
        })

        # Location model
        if(has_location(x)) {
            eta <- lapply(seq_len(S), function(s) {
                eta[[s]] + pred_pred_spec$x_loc[n,, drop = FALSE] %*% mu_beta[[s]]
            })
            if(has_random_location(x)) {
                eta <- lapply(seq_len(S), function(s) {
                    random_coef_matrix <- matrix(pred_random[[s]][re_indices$mu_beta_random], nrow = P_random, ncol = F)
                    eta[[s]] + pred_pred_spec$x_loc[n,P_random_ind, drop = FALSE] %*% random_coef_matrix
                })
            }
        }

        # Scale model
        if(has_scale(x)) {
            eta_logsd <- lapply(seq_len(S), function(s) {
                eta_logsd[[s]] + pred_pred_spec$x_sca[n,, drop = FALSE] %*% logsd_beta[[s]]
            })
            if(has_random_scale(x)) {
                eta_logsd <- lapply(seq_len(S), function(s) {
                    random_coef_matrix <- matrix(pred_random[[s]][re_indices$logsd_beta_random], nrow = Q_random, ncol = F)
                    eta_logsd[[s]] + pred_pred_spec$x_sca[n,Q_random_ind, drop = FALSE] %*% random_coef_matrix
                })
            }
        }

        # Add eta realizations
        if(include_error) {
            eta <- lapply(seq_len(S), function(s) {
                D <- diag(exp(eta_logsd[[s]][,]), F, F)
                cov_matrix <- D %*% factor_cor[[s]] %*% D
                eta[[s]] + mvrnorm(1, rep(0, F), Sigma = cov_matrix)
            })
        }
        # Return if predicting latents
        pred_list[[n]] <- nlist(eta, eta_logsd)
        
        # Indicators, if latent; if not latent, then y and eta are the same
        y <- NA
        if(has_latent(x)) {
            y <- lapply(seq_len(S), function(s) {
                eta[[s]] %*% lambda[[s]] + t(nu[[s]])
            })
            if(include_error) {
                y <- lapply(seq_len(S), function(s) {
                    y[[s]] + rnorm(J, 0, sigma[[s]])
                })
            }
        } else if(!has_latent(x)) {
            y <- pred_list[[n]][["eta"]]
        }
        pred_list[[n]][["y"]] <- y
    } # For each n in 1:N

    # Restructure
    for(n in seq_len(N)) {
        # Each list element to S-rows matrix
        pred_list[[n]] <- lapply(pred_list[[n]], function(l) {
            do.call(rbind, l)
        })

        # Rename
        eta_names <- paste0("eta[",n, ",", seq_len(F),"]")
        eta_logsd_names <- paste0("eta_logsd[",n, ",", seq_len(F),"]")
        y_names <- paste0("y[",n,",", seq_len(get_J(x)), "]")
        colnames(pred_list[[n]][["eta"]]) <- eta_names
        colnames(pred_list[[n]][["eta_logsd"]]) <- eta_logsd_names
        colnames(pred_list[[n]][["y"]]) <- y_names

        if(summarize) {
            pred_list[[n]][c("eta","eta_logsd")] <- lapply(pred_list[[n]][c("eta","eta_logsd")], .summarize, prob = prob)
            pred_list[[n]][c("eta","eta_logsd")] <- lapply(pred_list[[n]][c("eta","eta_logsd")], .tidy_summary, labs = c("observation", facVar), seq_len(N), get_factor_names(x))
            pred_list[[n]][c("eta","eta_logsd")] <- lapply(pred_list[[n]][c("eta","eta_logsd")], .summary_rearrange, cols = c("observation", facVar))

            pred_list[[n]][["y"]] <- .summarize(pred_list[[n]][["y"]], prob = prob)
            pred_list[[n]][["y"]] <- .tidy_summary(pred_list[[n]][["y"]], labs = c("observation", "indicator"), seq_len(N), get_indicator_names(x))
            pred_list[[n]][["y"]] <- .summary_rearrange(pred_list[[n]][["y"]], cols = c("observation", "indicator"))
        }
    }

    if(summarize) {
        # End goal: List of each type ([eta, eta_logsd, y]); rbinded across N
        out <- do.call(.list_zip, c(pred_list, f = rbind))
        names(out) <- c("eta","eta_logsd","y")
        return(out)
    }
    # Else...
    pred_list
}

.lambda_matrix <- function(F, J, vec) {
    matrix(vec, F, J)
}

# Converts fitted random effects to uncorrelated z
# Returns list of random_z's; one list per sample.
# Does not split into different REs.
.random_to_z <- function(fit) {
    ranef_pars <- c("mu_random","logsd_random")
    if(has_random_location(fit)) ranef_pars <- c(ranef_pars, "mu_beta_random")
    if(has_random_scale(fit)) ranef_pars <- c(ranef_pars, "logsd_beta_random")
    K <- get_K(fit)
    F <- get_F(fit)
    N_ranefs <- get_number_re(fit)

    # ranef_pars-length List of Samples-length Lists of KxRandom matrices
    ranefs_list <- lapply(ranef_pars, function(par) {.extract_transform_to_list(fit$fit, par)})
    names(ranefs_list) <- ranef_pars
    # Flatten out random coef matrices -- My method does too good a job. Can't cbind the array.
    if(has_random_location(fit)) {
        ranefs_list[["mu_beta_random"]] <- lapply(ranefs_list[["mu_beta_random"]], function(s) {
            t(apply(s, 1, function(r) {r}))
        })
    }
    if(has_random_scale(fit)) {
        ranefs_list[["logsd_beta_random"]] <- lapply(ranefs_list[["logsd_beta_random"]], function(s) {
            t(apply(s, 1, function(r) {r}))
        })
    }
    # List of Samples-length KxRandom matrices/arrays
    ranefs <- do.call(.list_zip, ranefs_list)

    random_sigma <- .extract_transform_to_list(fit$fit, par = "mu_logsd_betas_random_sigma")
    random_cor <- .extract_transform_to_list(fit$fit, par = "Omega_mean_logsd")
    random_cor_U <- lapply(random_cor, function(R){chol(R)})

    S <- get_S(fit)

    R <- get_R(fit)
    if(has_between(fit)) { # If between model
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
##' @method fitted lmmelsm
##' @export
fitted.lmmelsm <- function(object, prob = .95, summarize = TRUE, ...) {
    pars <- c("eta", "eta_logsd")
    facVar <- ifelse(has_latent(object), "factor", "variable")
    fnames <- get_factor_names(object)

    eta <- as.matrix(object$fit, pars = pars[1])
    eta_logsd <- as.matrix(object$fit, pars = pars[2])

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
