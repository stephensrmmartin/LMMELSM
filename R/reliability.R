icc <- function(x, ...) {
    UseMethod(icc, x)
}

icc.lmmelsm <- function(x, prob = .95, ...) {
    # Unconditional ICC(1): V_i(u^mu_i) / (V_i(u^mu_i) + V_i(epsilon_i))
    # Unconditional ICC(2): V_i(u^mu_i) / (V_i(u^mu_i) + V_i(epsilon_i) / n_(k_i))
    # Conditional ICC(1): V_i(u_0^mu_i) / (V_i(u_0^mu_i) + V_i(epsilon_i)); u_0^mu_i = location intercept; location covariates are present.
    # Adjusted, conditional ICC(1): V_i(RE_projections) / (V_i(RE_projections) + V_i(epsilon_i))
    # V_i(RE_projections) = V_i(Z_i u_i) = Z_i V_i(u_i) Z_i'
    # V_i(RE_projections) = [According to paper] E_i(...) = mean(Z_i V_i(u_i) Z_i'); I think we just did Z_i V_i(u_i) Z_i'; not a mean over Z_i.

    unconditional <- !has_location(x)
    if(!unconditional | has_between(x)) {
        stop("Only unconditional models with constant between-group variances are currently supported.")
    }

    F <- get_F(x)
    S <- get_S(x)
    gn <- x$meta$group_spec$numeric
    re_indices <- get_re_indices(x)

    if(has_between(x)) {
        x_bet <- x$stan_data$x_bet
        x_bet_l2 <- x_bet[match(seq_len(get_K(x), gn))]
    }

    between <- .extract_transform_to_list(x$fit, "mu_logsd_betas_random_sigma")
    between <- lapply(seq_len(S), function(s) {
        matrix(between[[s]], byrow = TRUE, nrow = get_K(x), ncol = get_number_re(x))
    })
    if(has_between(x)) {
        zeta <- .extract_transform_to_list(x, "zeta")
        between <- lapply(seq_len(S), function(s) {
            out <- between[[s]]
            out[[s]][, seq_len(2*F), drop = FALSE] <- out[[s]][, seq_len(2*F), drop = FALSE] * exp(x_bet_l2 %*% zeta[[s]])
        })
    }

    eta_logsd <- .extract_transform_to_list(x$fit, "eta_logsd")

    icc_unconditional <- lapply(seq_len(S), function(s) {
        between_vars <- (between[[s]][,re_indices$mu_random, drop = FALSE])^2 # Vars of RE Location Intercepts
        between_vars <- between_vars[gn,, drop = FALSE]
        iccs <- between_vars / (between_vars + exp(2 * eta_logsd[[s]]))
        iccs <- matrix(iccs, nrow = 1)
        cns <- paste0("icc[", seq_len(get_N(x)), ",", rep(seq_len(F), each = get_N(x)), "]")
        colnames(iccs) <- cns
        iccs
    })

    icc_unconditional <- do.call(rbind, icc_unconditional)

    icc_unconditional <- .summarize(icc_unconditional, prob = prob)
    facVar <- ifelse(has_latent(x), "factor", "variable")
    icc_unconditional <- .tidy_summary(icc_unconditional, labs = c("observation", facVar), seq_len(get_N(x)), get_factor_names(x))

    icc_unconditional
}

omega <- function(x, ...) {
    UseMethod(omega, x)
}

omega.lmmelsm <- function(x) {
    
}
