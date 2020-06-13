
##' @title loo method for LMMELSM objects.
##' @param x lmmelsm object.
##' @param ... Not used.
##' @return loo object.
##' @author Stephen R. Martin
##' @importFrom loo loo
##' @export loo
##' @export
loo.lmmelsm <- function(x, type = c("observation", "group"), ...) {
    type <- match.arg(type)

    ll <- .log_liks(x)

    ll_fun <- switch(type,
                     observation = .ll_obs,
                     group = .ll_group)

    ll <- ll_fun(ll, x)

    ## r_eff <- loo::relative_eff(exp(ll))
    ## looOut <- loo::loo(ll, r_eff = r_eff)
    looOut <- loo::loo(ll)

    return(looOut)
}

.ll_obs <- function(ll, ...) {

    out <- apply(ll, 3, function(x) { # Each sample
        rowSums(x) # Observation-LL; 
    }) # N x S
    out <- t(out) # Make S x N

    return(out)
    
}

.ll_group <- function(ll, x) {
    ll_obs <- .ll_obs(ll)
    ll_obs <- t(ll_obs) # Make N x S

    gn <- x$meta$group_spec$numeric
    N <- nrow(ll_obs)
    N_seq <- seq_len(N)
    ## out <- tapply(ll_obs, gn, colSums)
    out <- tapply(N_seq, gn, function(x) {
        colSums(ll_obs[x,])
    }, simplify = FALSE) # K-length list; ordered by 1:K
    out <- do.call(rbind, out) # N x S
    out <- t(out)  # S x N

    return(out)
}

.log_liks <- function(x) {
    J <- x$meta$indicator_spec$J
    N <- x$meta$indicator_spec$N
    y <- x$stan_data$y
    # Extract and transform param samples.
    etas <- .extract_transform(x$fit, par = "eta")
    lambdas <- .extract_transform(x$fit, par = "lambda")
    nus <- .extract_transform(x$fit, par = "nu")
    sigmas <- .extract_transform(x$fit, par = "sigma")

    S <- dim(etas)[3]

    # Predictions
    preds <- array(0, dim = c(N, J, S))
    for(s in seq_len(S)) {
        preds[,, s] <- .ones(N) %*% t(nus[, s]) + etas[,, s] %*% lambdas[,, s]
    }

    # Likelihoods
    log_lik <- array(0, dim = c(N, J, S))
    for(n in seq_len(N)) {
        log_lik[n, ,] <- dnorm(y[n, ], preds[n, , ], sigmas, log = TRUE)
    }

    return(log_lik)
}
